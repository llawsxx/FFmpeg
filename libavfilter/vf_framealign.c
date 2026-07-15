/*
 * Copyright (c) 2026
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 */

/**
 * @file
 * Align two video streams by comparing windows of frame fingerprints.
 */

#include <inttypes.h>
#include <math.h>
#include <string.h>

#include "libavutil/avstring.h"
#include "libavutil/error.h"
#include "libavutil/fifo.h"
#include "libavutil/imgutils.h"
#include "libavutil/mem.h"
#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"

#include "avfilter.h"
#include "filters.h"
#include "formats.h"
#include "video.h"

enum Metric {
    METRIC_SSIM,
    METRIC_PSNR,
};

enum FailAction {
    FAIL_ERROR,
    FAIL_WARN,
    FAIL_PASS,
};

enum AlignState {
    STATE_SEARCHING,
    STATE_LOCKED,
    STATE_PASS,
};

enum RefUVMode {
    REF_UV_NONE,
    REF_UV_TOP,
    REF_UV_BOTTOM,
    REF_UV_ALL,
};

typedef struct FrameItem {
    AVFrame *frame;
    uint16_t *fingerprint;
} FrameItem;

typedef struct FrameAlignContext {
    const AVClass *class;

    int buffer_size;
    int window;
    int locked_window;
    int refine;
    int sample_w;
    int sample_h;
    int metric;
    double threshold;
    int fail_action;
    int max_searches;
    int ref_uv;

    AVFifo *queue[2];
    int input_status[2];
    int64_t input_status_pts[2];
    int request_turn;
    int state;
    int searches;
    int configured;
    int depth;
    int chroma_h;
    int chroma_linesize;
    int fingerprint_size;
    uint16_t *line;
    double *score_matrix;
    size_t score_matrix_elems;
    int64_t last_pts;
    double last_score;
    int64_t aligned_pairs;
} FrameAlignContext;

#define OFFSET(x) offsetof(FrameAlignContext, x)
#define FLAGS (AV_OPT_FLAG_VIDEO_PARAM | AV_OPT_FLAG_FILTERING_PARAM)

static const AVOption framealign_options[] = {
    { "buffer", "number of frames to buffer before a full search",
      OFFSET(buffer_size), AV_OPT_TYPE_INT, { .i64 = 100 }, 2, 10000, FLAGS },
    { "window", "number of consecutive frames used for each score",
      OFFSET(window), AV_OPT_TYPE_INT, { .i64 = 25 }, 1, 10000, FLAGS },
    { "locked_window", "number of consecutive frames checked after alignment",
      OFFSET(locked_window), AV_OPT_TYPE_INT, { .i64 = 10 }, 1, 10000, FLAGS },
    { "refine", "maximum forward slip searched after alignment",
      OFFSET(refine), AV_OPT_TYPE_INT, { .i64 = 5 }, 0, 1000, FLAGS },
    { "sample_w", "fingerprint width",
      OFFSET(sample_w), AV_OPT_TYPE_INT, { .i64 = 32 }, 4, 512, FLAGS },
    { "sample_h", "fingerprint height",
      OFFSET(sample_h), AV_OPT_TYPE_INT, { .i64 = 18 }, 4, 512, FLAGS },
    { "metric", "similarity metric", OFFSET(metric), AV_OPT_TYPE_INT,
      { .i64 = METRIC_SSIM }, METRIC_SSIM, METRIC_PSNR, FLAGS, .unit = "metric" },
        { "ssim", NULL, 0, AV_OPT_TYPE_CONST, { .i64 = METRIC_SSIM },
          0, 0, FLAGS, .unit = "metric" },
        { "psnr", NULL, 0, AV_OPT_TYPE_CONST, { .i64 = METRIC_PSNR },
          0, 0, FLAGS, .unit = "metric" },
    { "threshold", "minimum average window score required for alignment",
      OFFSET(threshold), AV_OPT_TYPE_DOUBLE, { .dbl = 0.90 }, -1.0, 100.0, FLAGS },
    { "fail_action", "action after alignment repeatedly fails",
      OFFSET(fail_action), AV_OPT_TYPE_INT, { .i64 = FAIL_ERROR },
      FAIL_ERROR, FAIL_PASS, FLAGS, .unit = "fail_action" },
        { "error", NULL, 0, AV_OPT_TYPE_CONST, { .i64 = FAIL_ERROR },
          0, 0, FLAGS, .unit = "fail_action" },
        { "warn", NULL, 0, AV_OPT_TYPE_CONST, { .i64 = FAIL_WARN },
          0, 0, FLAGS, .unit = "fail_action" },
        { "pass", NULL, 0, AV_OPT_TYPE_CONST, { .i64 = FAIL_PASS },
          0, 0, FLAGS, .unit = "fail_action" },
    { "max_searches", "failed full searches before fail_action, 0 is unlimited",
      OFFSET(max_searches), AV_OPT_TYPE_INT, { .i64 = 10 }, 0, INT_MAX, FLAGS },
    { "ref_uv", "replace main chroma with fields or all chroma from reference",
      OFFSET(ref_uv), AV_OPT_TYPE_INT, { .i64 = REF_UV_NONE },
      REF_UV_NONE, REF_UV_ALL, FLAGS, .unit = "ref_uv" },
        { "none", "keep all main chroma", 0, AV_OPT_TYPE_CONST,
          { .i64 = REF_UV_NONE }, 0, 0, FLAGS, .unit = "ref_uv" },
        { "top", "take top-field chroma from reference", 0, AV_OPT_TYPE_CONST,
          { .i64 = REF_UV_TOP }, 0, 0, FLAGS, .unit = "ref_uv" },
        { "bottom", "take bottom-field chroma from reference", 0, AV_OPT_TYPE_CONST,
          { .i64 = REF_UV_BOTTOM }, 0, 0, FLAGS, .unit = "ref_uv" },
        { "all", "take all chroma from reference", 0, AV_OPT_TYPE_CONST,
          { .i64 = REF_UV_ALL }, 0, 0, FLAGS, .unit = "ref_uv" },
    { NULL }
};

AVFILTER_DEFINE_CLASS(framealign);

static size_t queue_size(const FrameAlignContext *s, int input)
{
    return av_fifo_can_read(s->queue[input]);
}

static int queue_peek(const FrameAlignContext *s, int input, size_t offset,
                      FrameItem *item)
{
    return av_fifo_peek(s->queue[input], item, 1, offset);
}

static void free_item(FrameItem *item)
{
    av_frame_free(&item->frame);
    av_freep(&item->fingerprint);
}

static void queue_drop(FrameAlignContext *s, int input, size_t count)
{
    while (count--) {
        FrameItem item;

        if (av_fifo_read(s->queue[input], &item, 1) < 0)
            return;
        free_item(&item);
    }
}

static int make_fingerprint(FrameAlignContext *s, const AVFrame *frame,
                            uint16_t **fingerprint)
{
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(frame->format);
    const unsigned max = (1U << s->depth) - 1;
    uint16_t *dst;

    dst = av_malloc_array(s->fingerprint_size, sizeof(*dst));
    if (!dst)
        return AVERROR(ENOMEM);

    for (int oy = 0; oy < s->sample_h; oy++) {
        const int y = ((2 * oy + 1) * frame->height) / (2 * s->sample_h);

        av_read_image_line(s->line, (const uint8_t **)frame->data,
                           frame->linesize, desc, 0, y, 0, frame->width, 0);
        for (int ox = 0; ox < s->sample_w; ox++) {
            const int x = ((2 * ox + 1) * frame->width) / (2 * s->sample_w);
            dst[oy * s->sample_w + ox] = (s->line[x] * 255U + max / 2) / max;
        }
    }

    *fingerprint = dst;
    return 0;
}

static int queue_frame(FrameAlignContext *s, int input, AVFrame *frame)
{
    FrameItem item = { .frame = frame };
    int ret;

    ret = make_fingerprint(s, frame, &item.fingerprint);
    if (ret < 0) {
        av_frame_free(&frame);
        return ret;
    }

    ret = av_fifo_write(s->queue[input], &item, 1);
    if (ret < 0)
        free_item(&item);
    return ret;
}

static double fingerprint_ssim(const uint16_t *a, const uint16_t *b, int n)
{
    const double c1 = 2.55 * 2.55;
    const double c2 = 7.65 * 7.65;
    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    double mx, my, vx, vy, cov, score;

    for (int i = 0; i < n; i++) {
        sx  += a[i];
        sy  += b[i];
        sxx += a[i] * a[i];
        syy += b[i] * b[i];
        sxy += a[i] * b[i];
    }

    mx = sx / n;
    my = sy / n;
    vx = FFMAX(0.0, sxx / n - mx * mx);
    vy = FFMAX(0.0, syy / n - my * my);
    cov = sxy / n - mx * my;
    score = ((2 * mx * my + c1) * (2 * cov + c2)) /
            ((mx * mx + my * my + c1) * (vx + vy + c2));

    return av_clipd(score, -1.0, 1.0);
}

static double fingerprint_psnr(const uint16_t *a, const uint16_t *b, int n)
{
    double mse = 0.0;

    for (int i = 0; i < n; i++) {
        const double d = a[i] - b[i];
        mse += d * d;
    }
    mse /= n;

    return mse == 0.0 ? 100.0 : FFMIN(100.0, 10.0 * log10(255.0 * 255.0 / mse));
}

static double frame_score(const FrameAlignContext *s,
                          const FrameItem *a, const FrameItem *b)
{
    if (s->metric == METRIC_PSNR)
        return fingerprint_psnr(a->fingerprint, b->fingerprint,
                                s->fingerprint_size);
    return fingerprint_ssim(a->fingerprint, b->fingerprint,
                            s->fingerprint_size);
}

static int ensure_score_matrix(FrameAlignContext *s, size_t elems)
{
    double *matrix;

    if (elems <= s->score_matrix_elems)
        return 0;
    matrix = av_realloc_array(s->score_matrix, elems, sizeof(*matrix));
    if (!matrix)
        return AVERROR(ENOMEM);
    s->score_matrix = matrix;
    s->score_matrix_elems = elems;
    return 0;
}

static int find_alignment(FrameAlignContext *s, int count0, int count1,
                          int search_window, int max_start0, int max_start1,
                          int *best0, int *best1, double *best_score)
{
    const int starts0 = FFMIN(max_start0 + 1, count0 - search_window + 1);
    const int starts1 = FFMIN(max_start1 + 1, count1 - search_window + 1);
    const int used0 = starts0 + search_window - 1;
    const int used1 = starts1 + search_window - 1;
    int ret;

    if (starts0 <= 0 || starts1 <= 0)
        return AVERROR(EAGAIN);
    ret = ensure_score_matrix(s, (size_t)used0 * used1);
    if (ret < 0)
        return ret;

    for (int i = 0; i < used0; i++) {
        FrameItem a;

        queue_peek(s, 0, i, &a);
        for (int j = 0; j < used1; j++) {
            FrameItem b;

            queue_peek(s, 1, j, &b);
            s->score_matrix[i * used1 + j] = frame_score(s, &a, &b);
        }
    }

    *best_score = -INFINITY;
    *best0 = *best1 = 0;
    for (int i = 0; i < starts0; i++) {
        for (int j = 0; j < starts1; j++) {
            double score = 0.0;

            for (int k = 0; k < search_window; k++)
                score += s->score_matrix[(i + k) * used1 + j + k];
            score /= search_window;
            if (score > *best_score) {
                *best_score = score;
                *best0 = i;
                *best1 = j;
            }
        }
    }

    return 0;
}

static int alignment_failed(AVFilterContext *ctx, double score)
{
    FrameAlignContext *s = ctx->priv;

    if (s->fail_action == FAIL_ERROR) {
        av_log(ctx, AV_LOG_ERROR,
               "Could not align inputs: best %s score %.6f is below threshold %.6f\n",
               s->metric == METRIC_SSIM ? "SSIM" : "PSNR", score, s->threshold);
        return AVERROR_INVALIDDATA;
    }

    if (s->fail_action == FAIL_WARN)
        av_log(ctx, AV_LOG_WARNING,
               "Could not align inputs: best %s score %.6f is below threshold %.6f; passing paired frames through\n",
               s->metric == METRIC_SSIM ? "SSIM" : "PSNR", score, s->threshold);
    s->state = STATE_PASS;
    return 0;
}

static int full_search(AVFilterContext *ctx)
{
    FrameAlignContext *s = ctx->priv;
    const int count0 = FFMIN((int)queue_size(s, 0), s->buffer_size);
    const int count1 = FFMIN((int)queue_size(s, 1), s->buffer_size);
    int start0, start1, ret;
    double score;

    ret = find_alignment(s, count0, count1, s->window,
                         count0 - s->window, count1 - s->window,
                         &start0, &start1, &score);
    if (ret < 0)
        return ret;

    s->last_score = score;
    if (score >= s->threshold) {
        queue_drop(s, 0, start0);
        queue_drop(s, 1, start1);
        s->state = STATE_LOCKED;
        s->searches = 0;
        av_log(ctx, AV_LOG_INFO,
               "Aligned inputs at offsets %d and %d, average %s score %.6f\n",
               start0, start1, s->metric == METRIC_SSIM ? "SSIM" : "PSNR", score);
        return 0;
    }

    s->searches++;
    av_log(ctx, AV_LOG_DEBUG,
           "Alignment search %d did not reach threshold: best score %.6f\n",
           s->searches, score);
    if ((s->max_searches && s->searches >= s->max_searches) ||
        (s->input_status[0] && s->input_status[1]))
        return alignment_failed(ctx, score);

    queue_drop(s, 0, FFMAX(1, count0 - s->window + 1));
    queue_drop(s, 1, FFMAX(1, count1 - s->window + 1));
    return 0;
}

static void set_metadata(AVFrame *frame, double score, int aligned)
{
    char value[64];

    snprintf(value, sizeof(value), "%.6f", score);
    av_dict_set(&frame->metadata, "lavfi.framealign.score", value, 0);
    av_dict_set_int(&frame->metadata, "lavfi.framealign.aligned", aligned, 0);
}

static int replace_main_uv(FrameAlignContext *s, AVFrame *main,
                           const AVFrame *reference)
{
    int ret;

    if (s->ref_uv == REF_UV_NONE)
        return 0;

    ret = av_frame_make_writable(main);
    if (ret < 0)
        return ret;

    for (int plane = 1; plane <= 2; plane++) {
        for (int y = 0; y < s->chroma_h; y++) {
            if (s->ref_uv == REF_UV_TOP && (y & 1))
                continue;
            if (s->ref_uv == REF_UV_BOTTOM && !(y & 1))
                continue;
            memcpy(main->data[plane] + (ptrdiff_t)y * main->linesize[plane],
                   reference->data[plane] + (ptrdiff_t)y * reference->linesize[plane],
                   s->chroma_linesize);
        }
    }

    return 0;
}

static int output_pair(AVFilterContext *ctx, double score, int aligned)
{
    FrameAlignContext *s = ctx->priv;
    FrameItem item[2];
    int ret;

    ret = av_fifo_read(s->queue[0], &item[0], 1);
    if (ret < 0)
        return ret;
    ret = av_fifo_read(s->queue[1], &item[1], 1);
    if (ret < 0) {
        free_item(&item[0]);
        return ret;
    }

    av_freep(&item[0].fingerprint);
    av_freep(&item[1].fingerprint);
    ret = replace_main_uv(s, item[0].frame, item[1].frame);
    if (ret < 0) {
        av_frame_free(&item[0].frame);
        av_frame_free(&item[1].frame);
        return ret;
    }
    item[1].frame->pts = item[0].frame->pts;
    item[0].frame->time_base = ctx->outputs[0]->time_base;
    set_metadata(item[0].frame, score, aligned);
    s->last_pts = item[0].frame->pts;
    s->aligned_pairs += aligned;

    if (ctx->nb_outputs == 1) {
        av_frame_free(&item[1].frame);
        return ff_filter_frame(ctx->outputs[0], item[0].frame);
    }

    item[1].frame->time_base = ctx->outputs[1]->time_base;
    set_metadata(item[1].frame, score, aligned);
    ret = ff_filter_frame(ctx->outputs[0], item[0].frame);
    if (ret < 0) {
        av_frame_free(&item[1].frame);
        return ret;
    }
    return ff_filter_frame(ctx->outputs[1], item[1].frame);
}

static int try_locked_output(AVFilterContext *ctx)
{
    FrameAlignContext *s = ctx->priv;
    const int count0 = queue_size(s, 0);
    const int count1 = queue_size(s, 1);
    int start0 = 0, start1 = 0, ret;
    double score;

    if (!count0 || !count1)
        return AVERROR(EAGAIN);

    if (count0 < s->locked_window || count1 < s->locked_window) {
        if (s->input_status[0] || s->input_status[1])
            return output_pair(ctx, s->last_score, 1);
        return AVERROR(EAGAIN);
    }

    ret = find_alignment(s, count0, count1, s->locked_window, 0, 0,
                         &start0, &start1, &score);
    if (ret < 0)
        return ret;
    if (score >= s->threshold) {
        s->last_score = score;
        return output_pair(ctx, score, 1);
    }

    if (count0 >= s->locked_window + s->refine &&
        count1 >= s->locked_window + s->refine && s->refine > 0) {
        ret = find_alignment(s, count0, count1, s->locked_window,
                             s->refine, s->refine,
                             &start0, &start1, &score);
        if (ret < 0)
            return ret;
        if (score >= s->threshold) {
            queue_drop(s, 0, start0);
            queue_drop(s, 1, start1);
            s->last_score = score;
            av_log(ctx, AV_LOG_INFO,
                   "Corrected alignment by dropping %d and %d frames, score %.6f\n",
                   start0, start1, score);
            return output_pair(ctx, score, 1);
        }
    }

    av_log(ctx, AV_LOG_WARNING,
           "Alignment lost at output pair %"PRId64" (score %.6f); restarting full search\n",
           s->aligned_pairs, score);
    s->state = STATE_SEARCHING;
    s->searches = 0;
    s->last_score = score;
    return 0;
}

static int inputs_ready_for_search(const FrameAlignContext *s)
{
    const int count0 = queue_size(s, 0);
    const int count1 = queue_size(s, 1);

    return count0 >= s->window && count1 >= s->window &&
           (count0 >= s->buffer_size || s->input_status[0]) &&
           (count1 >= s->buffer_size || s->input_status[1]);
}

static int finish_if_done(AVFilterContext *ctx)
{
    FrameAlignContext *s = ctx->priv;

    if ((s->input_status[0] && !queue_size(s, 0)) ||
        (s->input_status[1] && !queue_size(s, 1))) {
        for (int i = 0; i < ctx->nb_outputs; i++)
            ff_outlink_set_status(ctx->outputs[i], AVERROR_EOF, s->last_pts);
        queue_drop(s, 0, queue_size(s, 0));
        queue_drop(s, 1, queue_size(s, 1));
        return 1;
    }
    return 0;
}

static int select_input(const FrameAlignContext *s, int target)
{
    const size_t count0 = queue_size(s, 0);
    const size_t count1 = queue_size(s, 1);
    const int need0 = !s->input_status[0] && count0 < target;
    const int need1 = !s->input_status[1] && count1 < target;

    if (need0 && need1) {
        if (count0 < count1)
            return 0;
        if (count1 < count0)
            return 1;
        return s->request_turn;
    }
    if (need0)
        return 0;
    if (need1)
        return 1;
    return -1;
}

static int activate(AVFilterContext *ctx)
{
    FrameAlignContext *s = ctx->priv;
    int wanted = 0;
    AVFrame *frame = NULL;
    int input, target, ret, status;
    int64_t pts;

    for (int i = 0; i < ctx->nb_outputs; i++) {
        wanted |= ff_outlink_frame_wanted(ctx->outputs[i]);
        FF_FILTER_FORWARD_STATUS_BACK_ALL(ctx->outputs[i], ctx);
    }

    if (s->state == STATE_SEARCHING && inputs_ready_for_search(s)) {
        ret = full_search(ctx);
        if (ret < 0)
            return ret;
    }

    if (s->state == STATE_SEARCHING &&
        ((s->input_status[0] && queue_size(s, 0) < s->window) ||
         (s->input_status[1] && queue_size(s, 1) < s->window))) {
        ret = alignment_failed(ctx, s->last_score);
        if (ret < 0)
            return ret;
    }

    if (wanted && s->state == STATE_LOCKED) {
        ret = try_locked_output(ctx);
        if (ret < 0 && ret != AVERROR(EAGAIN))
            return ret;
        if (ret != AVERROR(EAGAIN) && s->state == STATE_LOCKED)
            return ret;
    } else if (wanted && s->state == STATE_PASS &&
               queue_size(s, 0) && queue_size(s, 1)) {
        return output_pair(ctx, s->last_score, 0);
    }

    if (finish_if_done(ctx))
        return 0;

    if (s->state != STATE_SEARCHING && !wanted)
        return FFERROR_NOT_READY;

    target = s->state == STATE_SEARCHING ? s->buffer_size :
             s->state == STATE_LOCKED ? s->locked_window + s->refine : 1;
    input = select_input(s, target);
    if (input < 0)
        return FFERROR_NOT_READY;

    /* First consume a frame which is already queued on the selected link. */
    ret = ff_inlink_consume_frame(ctx->inputs[input], &frame);
    if (ret < 0)
        return ret;
    if (ret > 0) {
        ret = queue_frame(s, input, frame);
        if (ret < 0)
            return ret;
        s->request_turn = 1 - input;
        ff_filter_set_ready(ctx, 100);
        return 0;
    }

    /* No queued frame: acknowledge EOF/error before requesting upstream. */
    if (ff_inlink_acknowledge_status(ctx->inputs[input], &status, &pts)) {
        s->input_status[input] = status;
        s->input_status_pts[input] = pts;
        ff_filter_set_ready(ctx, 100);
        return 0;
    }

    /* consume_frame() returned no frame and the input is still live. Only
     * now ask this one upstream link for another frame. */
    ff_inlink_request_frame(ctx->inputs[input]);
    s->request_turn = 1 - input;
    return 0;
}

static int config_output(AVFilterLink *outlink)
{
    AVFilterContext *ctx = outlink->src;
    FrameAlignContext *s = ctx->priv;
    AVFilterLink *main = ctx->inputs[0];
    FilterLink *ol = ff_filter_link(outlink);
    FilterLink *ml = ff_filter_link(main);
    const AVPixFmtDescriptor *desc;

    if (ctx->inputs[0]->w != ctx->inputs[1]->w ||
        ctx->inputs[0]->h != ctx->inputs[1]->h) {
        av_log(ctx, AV_LOG_ERROR, "Input dimensions must match\n");
        return AVERROR(EINVAL);
    }
    if (ctx->inputs[0]->format != ctx->inputs[1]->format) {
        av_log(ctx, AV_LOG_ERROR, "Input pixel formats must match\n");
        return AVERROR(EINVAL);
    }

    outlink->w = main->w;
    outlink->h = main->h;
    outlink->time_base = main->time_base;
    outlink->sample_aspect_ratio = main->sample_aspect_ratio;
    ol->frame_rate = ml->frame_rate;

    if (s->configured)
        return 0;
    desc = av_pix_fmt_desc_get(main->format);
    if (!desc || !desc->nb_components)
        return AVERROR(EINVAL);
    if (s->ref_uv != REF_UV_NONE &&
        (desc->nb_components < 3 ||
         desc->flags & AV_PIX_FMT_FLAG_RGB ||
         desc->comp[0].plane != 0 ||
         desc->comp[1].plane != 1 ||
         desc->comp[2].plane != 2)) {
        av_log(ctx, AV_LOG_ERROR,
               "ref_uv requires a planar YUV pixel format with separate U and V planes\n");
        return AVERROR(EINVAL);
    }
    s->depth = desc->comp[0].depth;
    if (s->ref_uv != REF_UV_NONE) {
        s->chroma_h = AV_CEIL_RSHIFT(main->h, desc->log2_chroma_h);
        s->chroma_linesize = av_image_get_linesize(main->format, main->w, 1);
        if (s->chroma_linesize < 0)
            return s->chroma_linesize;
    }
    s->line = av_malloc_array(main->w, sizeof(*s->line));
    if (!s->line)
        return AVERROR(ENOMEM);
    s->configured = 1;
    return 0;
}

static int query_formats(const AVFilterContext *ctx,
                         AVFilterFormatsConfig **cfg_in,
                         AVFilterFormatsConfig **cfg_out)
{
    return ff_set_common_formats2(ctx, cfg_in, cfg_out,
                                  ff_formats_pixdesc_filter(0,
                                      AV_PIX_FMT_FLAG_BITSTREAM |
                                      AV_PIX_FMT_FLAG_HWACCEL |
                                      AV_PIX_FMT_FLAG_PAL));
}

static av_cold int init(AVFilterContext *ctx)
{
    FrameAlignContext *s = ctx->priv;
    AVFilterPad pad = {
        .type = AVMEDIA_TYPE_VIDEO,
        .config_props = config_output,
    };
    int ret;

    if (s->window > s->buffer_size) {
        av_log(ctx, AV_LOG_ERROR, "window must not be larger than buffer\n");
        return AVERROR(EINVAL);
    }
    if (s->metric == METRIC_PSNR && s->threshold <= 1.0)
        av_log(ctx, AV_LOG_WARNING,
               "PSNR thresholds are in dB; threshold %.3f is probably too low\n",
               s->threshold);

    pad.name = "main";
    ret = ff_append_outpad(ctx, &pad);
    if (ret < 0)
        return ret;
    if (s->ref_uv == REF_UV_NONE) {
        pad.name = "reference";
        ret = ff_append_outpad(ctx, &pad);
        if (ret < 0)
            return ret;
    }

    s->fingerprint_size = s->sample_w * s->sample_h;
    s->last_pts = AV_NOPTS_VALUE;
    s->last_score = -INFINITY;
    s->state = STATE_SEARCHING;
    for (int i = 0; i < 2; i++) {
        s->queue[i] = av_fifo_alloc2(s->buffer_size + s->window + s->refine,
                                     sizeof(FrameItem), AV_FIFO_FLAG_AUTO_GROW);
        if (!s->queue[i])
            return AVERROR(ENOMEM);
    }
    return 0;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    FrameAlignContext *s = ctx->priv;

    for (int i = 0; i < 2; i++) {
        if (s->queue[i])
            queue_drop(s, i, queue_size(s, i));
        av_fifo_freep2(&s->queue[i]);
    }
    av_freep(&s->line);
    av_freep(&s->score_matrix);
}

static const AVFilterPad framealign_inputs[] = {
    {
        .name = "main",
        .type = AVMEDIA_TYPE_VIDEO,
    }, {
        .name = "reference",
        .type = AVMEDIA_TYPE_VIDEO,
    },
};

const FFFilter ff_vf_framealign = {
    .p.name        = "framealign",
    .p.description = NULL_IF_CONFIG_SMALL("Align two video streams by frame similarity."),
    .p.priv_class  = &framealign_class,
    .p.outputs     = NULL,
    .p.flags       = AVFILTER_FLAG_DYNAMIC_OUTPUTS,
    .priv_size     = sizeof(FrameAlignContext),
    .init          = init,
    .uninit        = uninit,
    .activate      = activate,
    FILTER_INPUTS(framealign_inputs),
    FILTER_QUERY_FUNC2(query_formats),
};
