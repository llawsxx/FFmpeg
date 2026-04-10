#include "libavutil/avassert.h"
#include "libavutil/opt.h"
#include "libavutil/frame.h"
#include "libavutil/imgutils.h"
#include "libavutil/pixdesc.h"
#include "libavutil/mem.h"
#include "libavutil/rational.h"
#include "libavfilter/formats.h"
#include "avfilter.h"
#include "filters.h"
#include "video.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include "rife/rife_wrapper.h"

typedef struct RIFEContext {
    const AVClass* class;
    AVRational dest_frame_rate;         ///< output frames per second
    AVRational srce_time_base;          ///< timebase of source
    AVRational dest_time_base;          ///< timebase of destination

    char* model_path;
    int gpu_id;
    int gpu_thread;
    int tta;
    int tta_temporal_mode;
    int uhd;
    int factor;
    int full_interp;
    double max_ts_jump;
    double scene_threshold;
    double score;
    double prev_mafd;

    RIFEWrapperHandle rife;

    int rife_version;
    int interpolation_mode;
    int padding;

    AVFrame* work;

    AVFrame* f0;                        ///< last frame
    AVFrame* f1;                        ///< current frame
    int64_t pts0;                       ///< last frame pts in dest_time_base
    int64_t pts1;                       ///< current frame pts in dest_time_base
    int f0_used;
    int f1_used;
    double pts_diff_max;
    int64_t delta;                      ///< pts1 to pts0 delta
    int64_t max_delta;                  ///< max jump delta
    int64_t start_pts;                  ///< pts of the first output frame
    int64_t n;
    int cur_interp_i;
    int flush;
} RIFEContext;

#define OFFSET(x) offsetof(RIFEContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM | AV_OPT_FLAG_VIDEO_PARAM

static const AVOption rife_options[] = {
    { "model_path", "path to model directory", OFFSET(model_path), AV_OPT_TYPE_STRING, {.str = NULL }, 0, 0, FLAGS },
    { "version", "model version", OFFSET(rife_version), AV_OPT_TYPE_INT, {.i64 = 0 }, 0, 1, FLAGS, .unit = "version"},
        {"v2",   "RIFE_v2",   0, AV_OPT_TYPE_CONST, {.i64 = 0}, 0, 0, FLAGS, .unit = "version"},
        {"v4",   "RIFE_v4",  0, AV_OPT_TYPE_CONST,  {.i64 = 1}, 0, 0, FLAGS, .unit = "version"},
    { "mode", "interpolation mode", OFFSET(interpolation_mode), AV_OPT_TYPE_INT, {.i64 = 1 }, 0, 1, FLAGS, .unit = "mode"},
        {"fps",    "target fps mode",   0, AV_OPT_TYPE_CONST, {.i64 = 0}, 0, 0, FLAGS, .unit = "mode"},
        {"factor", "factor mode (default 2x)",  0, AV_OPT_TYPE_CONST,  {.i64 = 1}, 0, 0, FLAGS, .unit = "mode"},
    { "fps",  "required output frames per second rate", OFFSET(dest_frame_rate), AV_OPT_TYPE_VIDEO_RATE, {.str = "60"}, 0, INT_MAX, FLAGS},
    { "factor", "factor numerator", OFFSET(factor), AV_OPT_TYPE_INT, {.i64 = 2 }, 2, INT_MAX, FLAGS },
    { "gpu_id", "gpu device id", OFFSET(gpu_id), AV_OPT_TYPE_INT, {.i64 = 0 }, -1, INT_MAX, FLAGS },
    { "gpu_thread", "gpu thread count", OFFSET(gpu_thread), AV_OPT_TYPE_INT, {.i64 = 2 }, 1, INT_MAX, FLAGS },
    { "tta", "tta mode", OFFSET(tta), AV_OPT_TYPE_BOOL, {.i64 = 0 }, 0, 1, FLAGS },
    { "tta_temporal_mode", "tta temporal mode", OFFSET(tta_temporal_mode), AV_OPT_TYPE_BOOL, {.i64 = 0 }, 0, 1, FLAGS },
    { "uhd", "uhd mode", OFFSET(uhd), AV_OPT_TYPE_BOOL, {.i64 = 0 }, 0, 1, FLAGS },
    { "max_ts_jump", "timestamp jump threshold in seconds", OFFSET(max_ts_jump), AV_OPT_TYPE_DOUBLE, {.dbl = 1.0 }, 0.0, 3600.0, FLAGS },
    { "scene_threshold", "scene change threshold", OFFSET(scene_threshold), AV_OPT_TYPE_DOUBLE, {.dbl = 8.2 }, 0, 100, FLAGS },
    { "full_interp", "enable full interpolation mode", OFFSET(full_interp), AV_OPT_TYPE_BOOL, {.i64 = 0 }, 0, 1, FLAGS },
    { "padding", "image padding", OFFSET(padding), AV_OPT_TYPE_INT, {.i64 = 32 }, 32, 128, FLAGS },
    { NULL }
};

AVFILTER_DEFINE_CLASS(rife);

static void rife_planes(const AVFrame* f, const float** r, const float** g, const float** b, ptrdiff_t* stride)
{
    *g = (const float*)f->data[0];
    *b = (const float*)f->data[1];
    *r = (const float*)f->data[2];
    *stride = f->linesize[0] / (ptrdiff_t)sizeof(float);
}

static void rife_planes_w(AVFrame* f, float** r, float** g, float** b, ptrdiff_t* stride)
{
    *g = (float*)f->data[0];
    *b = (float*)f->data[1];
    *r = (float*)f->data[2];
    *stride = f->linesize[0] / (ptrdiff_t)sizeof(float);
}

static int rife_process_pair(RIFEContext* s, const AVFrame* src0, const AVFrame* src1, AVFrame* dst, float t)
{
    const float* src0R, * src0G, * src0B;
    const float* src1R, * src1G, * src1B;
    float* dstR, * dstG, * dstB;
    ptrdiff_t stride0, stride1, strided;

    rife_planes(src0, &src0R, &src0G, &src0B, &stride0);
    rife_planes(src1, &src1R, &src1G, &src1B, &stride1);
    rife_planes_w(dst, &dstR, &dstG, &dstB, &strided);

    return rife_wrapper_process(
        s->rife,
        src0R, src0G, src0B,
        src1R, src1G, src1B,
        dstR, dstG, dstB,
        src0->width, src0->height, stride0, stride1 , strided, t
    );
}

static double rife_scene_score(RIFEContext* s, const AVFrame* a, const AVFrame* b)
{
    const float* ar, * ag, * ab, * br, * bg, * bb;
    ptrdiff_t sa, sb;
    rife_planes(a, &ar, &ag, &ab, &sa);
    rife_planes(b, &br, &bg, &bb, &sb);

    double sum = 0.0;
    int cnt = 0;
    int step = 8;
    double mafd, diff, ret;
    for (int y = 0; y < a->height; y += step) {
        for (int x = 0; x < a->width; x += step) {
            int idx_a = y * sa + x;
            int idx_b = y * sb + x;
            double d = (fabs(ar[idx_a] - br[idx_b]) +
                fabs(ag[idx_a] - bg[idx_b]) +
                fabs(ab[idx_a] - bb[idx_b])) / 3.0;
            sum += d;
            cnt++;
        }
    }

    if (cnt <= 0)
        return 0.0;

    mafd = sum * 100 / cnt;
    diff = fabs(mafd - s->prev_mafd);
    ret = av_clipf(FFMIN(mafd, diff), 0, 100.);
    s->prev_mafd = mafd;

    return ret;
}

static int rife_send_interp(AVFilterContext* ctx, const AVFrame* a, const AVFrame* b, double t)
{
    AVFilterLink* outlink = ctx->outputs[0];
    RIFEContext* s = ctx->priv;
    AVFrame* mid = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    int ret;

    if (!mid)
        return AVERROR(ENOMEM);

    ret = rife_process_pair(s, a, b, mid, (float)t);
    if (ret != 0) {
        av_frame_free(&mid);
        return AVERROR_EXTERNAL;
    }
    av_frame_copy_props(mid, a);
    s->work = mid;
    return 0;
}

static int64_t compute_work_pts(RIFEContext* s, int64_t n) {
    return s->start_pts + av_rescale_q(n, av_inv_q(s->dest_frame_rate), s->dest_time_base);
}


static int process_work_frame_factor_mode(AVFilterContext* ctx)
{
    RIFEContext* s = ctx->priv;
    int ret;
    int64_t work_pts;

    if (!s->f1)
        return 0;
    if (!s->f0 && !s->flush)
        return 0;

    if (!s->f0) {
        av_assert1(s->flush);
        s->work = s->f1;
        work_pts = s->pts1;
        s->f1 = NULL;
    }
    else {
        if (s->cur_interp_i == 0) {
            s->work = av_frame_clone(s->f0);
            work_pts = s->pts0;
            s->cur_interp_i++;
        }
        else{
            int icnt = s->factor - 1;
            
            if (s->cur_interp_i <= icnt) {
                int use_scene_detection = (s->scene_threshold > 0.0);
                int is_scene_change = 0;
                if (use_scene_detection) {
                    if (s->score < 0.0) {
                        s->score = rife_scene_score(s, s->f0, s->f1);
                    }
                    is_scene_change = (s->score >= s->scene_threshold);

                    if (is_scene_change) {
                        av_log(ctx, AV_LOG_INFO, "Scene change detected: score=%.2f, threshold=%.2f\n",
                            s->score, s->scene_threshold);
                    }
                }

                double t = (double)s->cur_interp_i / (double)(icnt + 1);
                if (use_scene_detection && is_scene_change) {
                    s->work = av_frame_clone(s->f0);
                }
                else {
                    ret = rife_send_interp(ctx, s->f0, s->f1, t);
                    if (ret < 0)
                        return ret;
                }
                work_pts = (int64_t)(s->pts0 + s->delta * t);
                s->cur_interp_i++;
            }
            else {
                if (s->flush) {
                    s->work = s->f1;
                    work_pts = s->pts1;
                    s->f1 = NULL;
                }
                else {
                    return 0;
                }
            }

        }

    }

    if (!s->work)
        return AVERROR(ENOMEM);

    s->work->pts = work_pts;
    s->n++;

    return 1;
}


static int process_work_frame(AVFilterContext* ctx)
{
    RIFEContext* s = ctx->priv;
    int64_t work_pts;
    double interpolate;
    int ret;

    if (!s->f1)
        return 0;
    if (!s->f0 && !s->flush)
        return 0;

    work_pts = compute_work_pts(s, s->n);

    if (work_pts >= s->pts1 && !s->flush)
        return 0;

    if (!s->f0) {
        av_assert1(s->flush);
        s->work = s->f1;
        s->f1 = NULL;
    }
    else {
        if (work_pts >= s->pts1 + s->delta && s->flush)
            return 0;

        interpolate = (double)(work_pts - s->pts0) / s->delta;
        double work_pts_diff_pts0 = (double)FFABS(work_pts - s->pts0);
        double work_pts_diff_pts1 = (double)FFABS(work_pts - s->pts1);
        if (interpolate <= 0.0 || (!s->full_interp && !s->f0_used && work_pts_diff_pts0 < s->pts_diff_max)) {
            s->work = av_frame_clone(s->f0);
            s->f0_used = 1;
        }
        else if (interpolate >= 1.0 || (!s->full_interp && work_pts_diff_pts1 < s->pts_diff_max)) {
            s->work = av_frame_clone(s->f1);
            s->f1_used = 1;
        }
        else {
            int use_scene_detection = (s->scene_threshold > 0.0);
            int is_scene_change = 0;

            if (use_scene_detection) {
                if (s->score < 0.0) {
                    s->score = rife_scene_score(s, s->f0, s->f1);
                }
                is_scene_change = (s->score >= s->scene_threshold);

                if (is_scene_change) {
                    av_log(ctx, AV_LOG_INFO, "Scene change detected: score=%.2f, threshold=%.2f\n",
                        s->score, s->scene_threshold);
                }
            }

            if (use_scene_detection && is_scene_change) {
                s->work = av_frame_clone(s->f0);
            }
            else {
                ret = rife_send_interp(ctx, s->f0, s->f1, interpolate);
                if (ret < 0)
                    return ret;
            }
        }
    }

    if (!s->work)
        return AVERROR(ENOMEM);

    s->work->pts = work_pts;
    s->n++;

    return 1;
}

static int activate(AVFilterContext* ctx)
{
    int ret, status;
    AVFilterLink* inlink = ctx->inputs[0];
    AVFilterLink* outlink = ctx->outputs[0];
    RIFEContext* s = ctx->priv;
    AVFrame* inpicref;
    int64_t pts;

    FF_FILTER_FORWARD_STATUS_BACK(outlink, inlink);

retry:
    if(s->interpolation_mode == 0)
        ret = process_work_frame(ctx);
    else {
        ret = process_work_frame_factor_mode(ctx);
    }
    if (ret < 0)
        return ret;
    else if (ret == 1)
        return ff_filter_frame(outlink, s->work);

    ret = ff_inlink_consume_frame(inlink, &inpicref);
    if (ret < 0)
        return ret;

    if (inpicref) {
        if (inpicref->flags & AV_FRAME_FLAG_INTERLACED)
            av_log(ctx, AV_LOG_WARNING, "Interlaced frame found - the output will not be correct.\n");

        if (inpicref->pts == AV_NOPTS_VALUE) {
            av_log(ctx, AV_LOG_WARNING, "Ignoring frame without PTS.\n");
            av_frame_free(&inpicref);
        }
    }

    if (inpicref) {
        pts = av_rescale_q(inpicref->pts, s->srce_time_base, s->dest_time_base);

        if (s->f1 && pts == s->pts1) {
            av_log(ctx, AV_LOG_WARNING, "Ignoring frame with same PTS.\n");
            av_frame_free(&inpicref);
        }
    }

    if (inpicref) {
        av_frame_free(&s->f0);
        s->f0 = s->f1;
        s->pts0 = s->pts1;
        s->f1 = inpicref;
        s->pts1 = pts;
        s->f0_used = s->f1_used;
        s->f1_used = 0;
        s->score = -1.0;
        s->cur_interp_i = 0;
        if (s->pts0 != AV_NOPTS_VALUE)
            s->delta = s->pts1 - s->pts0;

        if (s->delta < 0 || s->delta > s->max_delta) {
            if(s->delta < 0)
                av_log(ctx, AV_LOG_WARNING, "PTS discontinuity\n");
            else
                av_log(ctx, AV_LOG_WARNING, "PTS difference is too large, consider increase max_ts_jump\n");
            s->start_pts = s->pts1;
            s->n = 0;
            if (s->f0) {
                if(s->f0_used){
                    av_frame_free(&s->f0);
                }else{
                    s->work = s->f0;
                    s->work->pts = s->pts0;
                    s->f0 = NULL;
                    return ff_filter_frame(outlink, s->work);
                }
            }
        }

        if (s->start_pts == AV_NOPTS_VALUE)
            s->start_pts = s->pts1;

        goto retry;
    }

    if (ff_inlink_acknowledge_status(inlink, &status, &pts)) {
        if (!s->flush) {
            s->flush = 1;
            goto retry;
        }
        ff_outlink_set_status(outlink, status, pts);
        return 0;
    }

    FF_FILTER_FORWARD_WANTED(outlink, inlink);

    return FFERROR_NOT_READY;
}


static int config_output(AVFilterLink* outlink)
{
    AVFilterContext* ctx = outlink->src;
    AVFilterLink* inlink = ctx->inputs[0];
    FilterLink* il = ff_filter_link(inlink);
    FilterLink* ol = ff_filter_link(outlink);
    RIFEContext* s = ctx->priv;
    int exact = 1;

    s->srce_time_base = inlink->time_base;

    if (s->interpolation_mode == 0) {
        // make sure timebase is small enough to hold the framerate
        exact = av_reduce(&s->dest_time_base.num, &s->dest_time_base.den,
            av_gcd((int64_t)s->srce_time_base.num * s->dest_frame_rate.num,
                (int64_t)s->srce_time_base.den * s->dest_frame_rate.den),
            (int64_t)s->srce_time_base.den * s->dest_frame_rate.num, INT_MAX);
        if (!exact) {
            av_log(ctx, AV_LOG_WARNING, "Timebase conversion is not exact\n");
        }
        ol->frame_rate = s->dest_frame_rate;
        outlink->time_base = s->dest_time_base;
    }
    else {
        s->dest_frame_rate = ol->frame_rate = av_mul_q(il->frame_rate, (AVRational) { s->factor, 1 });
        s->dest_time_base = outlink->time_base = av_mul_q(inlink->time_base, (AVRational) { 1, s->factor });
    }

    s->pts_diff_max = (double)av_rescale_q(1, av_inv_q(s->dest_frame_rate), s->dest_time_base) * 0.6;
    s->max_delta = s->max_ts_jump * av_q2d(av_inv_q(s->dest_time_base));

    av_log(ctx, AV_LOG_INFO,
        "time base:%u/%u -> %u/%u exact:%d\n",
        s->srce_time_base.num, s->srce_time_base.den,
        s->dest_time_base.num, s->dest_time_base.den, exact);

    av_log(ctx, AV_LOG_INFO, "fps:%u/%u max_ts_jump:%lf scene_threshold:%lf\n",
        s->dest_frame_rate.num, s->dest_frame_rate.den, s->max_ts_jump, 
        s->scene_threshold);

    return 0;
}

static av_cold int rife_init(AVFilterContext* ctx)
{
    RIFEContext* s = ctx->priv;
    int ret = 0;

    s->f0 = s->f1 = NULL;
    s->start_pts = AV_NOPTS_VALUE;
    s->pts0 = AV_NOPTS_VALUE;
    s->pts1 = AV_NOPTS_VALUE;
    s->prev_mafd = 0.0;
    int rife_v2 = 1;
    int rife_v4 = 0;
    if (s->rife_version == 0) {
        rife_v2 = 1;
        rife_v4 = 0;
    }
    else if (s->rife_version == 1) {
        rife_v2 = 0;
        rife_v4 = 1;
    }

    if (s->model_path == NULL) {
        av_log(ctx, AV_LOG_ERROR,
            "model path is empty\n");
        return AVERROR(EINVAL);
    }

    if (rife_v2) {
        if(s->interpolation_mode == 0){
            av_log(ctx, AV_LOG_ERROR,
                "v2 version does not support fps mode\n");
            return AVERROR(EINVAL);
        }
        if(s->interpolation_mode == 1 && s->factor != 2){
            av_log(ctx, AV_LOG_ERROR,
                "v2 version only support factor 2\n");
            return AVERROR(EINVAL);
        }
    }


    s->rife = rife_wrapper_create(
        s->gpu_id, s->tta, s->tta_temporal_mode, s->uhd, s->gpu_thread,
        rife_v2, rife_v4, s->padding
    );
    if (!s->rife)
        return AVERROR_EXTERNAL;

#ifdef _WIN32
    {
        int len = MultiByteToWideChar(CP_ACP, 0, s->model_path, -1, NULL, 0);
        if (len <= 0)
            return AVERROR(EINVAL);
        {
            wchar_t* wpath = (wchar_t*)av_mallocz(len * sizeof(wchar_t));
            if (!wpath)
                return AVERROR(ENOMEM);
            MultiByteToWideChar(CP_ACP, 0, s->model_path, -1, wpath, len);
            ret = rife_wrapper_load(s->rife, wpath);
            av_free(wpath);
        }
    }
#else
    ret = rife_wrapper_load(s->rife, s->model_path);
#endif

    if (ret != 0)
        return AVERROR_EXTERNAL;

    return 0;
}

static av_cold void rife_uninit(AVFilterContext* ctx)
{
    RIFEContext* s = ctx->priv;
    av_frame_free(&s->f0);
    av_frame_free(&s->f1);
    if (s->rife)
        rife_wrapper_destroy(s->rife);
    av_freep(&s->model_path);
}

static const AVFilterPad rife_inputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO
    },
};

static const AVFilterPad rife_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
        .config_props = config_output
    },
};

static const enum AVPixelFormat pix_fmts[] = {
    AV_PIX_FMT_GBRPF32LE,
    AV_PIX_FMT_NONE
};

const FFFilter ff_vf_rife = {
    .p.name = "rife",
    .p.description = NULL_IF_CONFIG_SMALL("RIFE frame interpolation filter"),
    .p.priv_class = &rife_class,
    .p.flags = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC,
    .priv_size = sizeof(RIFEContext),
    .init = rife_init,
    .uninit = rife_uninit,
    FILTER_INPUTS(rife_inputs),
    FILTER_OUTPUTS(rife_outputs),
    FILTER_PIXFMTS_ARRAY(pix_fmts),
    .activate = activate
};
