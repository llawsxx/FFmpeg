#include "libavutil/avassert.h"
#include "libavutil/opt.h"
#include "libavutil/frame.h"
#include "avfilter.h"
#include "filters.h"
#include "video.h"

#include <stdint.h>

typedef struct FPSSmoothContext {
    const AVClass* class;
    AVRational dest_frame_rate;         ///< output frames per second
    AVRational srce_time_base;          ///< timebase of source
    AVRational dest_time_base;          ///< timebase of destination

    int factor;
    double max_ts_jump;

    int interpolation_mode;

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
} FPSSmoothContext;

#define OFFSET(x) offsetof(FPSSmoothContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM | AV_OPT_FLAG_VIDEO_PARAM

static const AVOption fps_smooth_options[] = {
    { "fps",  "required output frames per second rate", OFFSET(dest_frame_rate), AV_OPT_TYPE_VIDEO_RATE, {.str = "25"}, 0, INT_MAX, FLAGS},
    { "factor", "factor numerator", OFFSET(factor), AV_OPT_TYPE_INT, {.i64 = 2 }, 2, INT_MAX, FLAGS },
    { "mode", "interpolation mode", OFFSET(interpolation_mode), AV_OPT_TYPE_INT, {.i64 = 0 }, 0, 1, FLAGS, .unit = "mode"},
        {"fps",    "target fps mode",   0, AV_OPT_TYPE_CONST, {.i64 = 0}, 0, 0, FLAGS, .unit = "mode"},
        {"factor", "factor mode (default 2x)",  0, AV_OPT_TYPE_CONST,  {.i64 = 1}, 0, 0, FLAGS, .unit = "mode"},
    { "max_ts_jump", "timestamp jump threshold in seconds", OFFSET(max_ts_jump), AV_OPT_TYPE_DOUBLE, {.dbl = 1.0 }, 0.0, 3600.0, FLAGS },
    { NULL }
};

AVFILTER_DEFINE_CLASS(fps_smooth);


static int fps_smooth_send_interp(AVFilterContext* ctx, const AVFrame* a, const AVFrame* b, double t)
{
    AVFilterLink* outlink = ctx->outputs[0];
    FPSSmoothContext* s = ctx->priv;
    AVFrame* mid = NULL;
    AVFrame* src = NULL;
    int ret;

    src = t <= 0.5 ? a : b;

    mid = av_frame_clone(src);
    if (!mid) {
        return AVERROR(ENOMEM);
    }

    av_frame_copy_props(mid, src);
    s->work = mid;
    return 0;
}

static int64_t compute_work_pts(FPSSmoothContext* s, int64_t n) {
    return s->start_pts + av_rescale_q(n, av_inv_q(s->dest_frame_rate), s->dest_time_base);
}


static int process_work_frame_factor_mode(AVFilterContext* ctx)
{
    FPSSmoothContext* s = ctx->priv;
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
                double t = (double)s->cur_interp_i / (double)(icnt + 1);
                ret = fps_smooth_send_interp(ctx, s->f0, s->f1, t);
                if (ret < 0)
                    return ret;
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
    FPSSmoothContext* s = ctx->priv;
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
        if (interpolate <= 0.0 || (!s->f0_used && work_pts_diff_pts0 < s->pts_diff_max)) {
            s->work = av_frame_clone(s->f0);
            s->f0_used = 1;
        }
        else if (interpolate >= 1.0 || (work_pts_diff_pts1 < s->pts_diff_max)) {
            s->work = av_frame_clone(s->f1);
            s->f1_used = 1;
        }
        else {
            ret = fps_smooth_send_interp(ctx, s->f0, s->f1, interpolate);
            if (ret < 0)
                return ret;
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
    FPSSmoothContext* s = ctx->priv;
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
    FPSSmoothContext* s = ctx->priv;
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

    av_log(ctx, AV_LOG_INFO, "fps:%u/%u max_ts_jump:%lf\n",
        s->dest_frame_rate.num, s->dest_frame_rate.den, s->max_ts_jump);

    return 0;
}

static av_cold int fps_smooth_init(AVFilterContext* ctx)
{
    FPSSmoothContext* s = ctx->priv;
    int ret = 0;

    s->f0 = s->f1 = NULL;
    s->start_pts = AV_NOPTS_VALUE;
    s->pts0 = AV_NOPTS_VALUE;
    s->pts1 = AV_NOPTS_VALUE;
    return 0;
}

static av_cold void fps_smooth_uninit(AVFilterContext* ctx)
{
    FPSSmoothContext* s = ctx->priv;
    av_frame_free(&s->f0);
    av_frame_free(&s->f1);
}

static const AVFilterPad fps_smooth_inputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO
    },
};

static const AVFilterPad fps_smooth_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
        .config_props = config_output
    },
};

const FFFilter ff_vf_fps_smooth = {
    .p.name = "fps_smooth",
    .p.description = NULL_IF_CONFIG_SMALL("Framerate conversion filter for unstable timestamp video."),
    .p.priv_class = &fps_smooth_class,
    .p.flags = AVFILTER_FLAG_METADATA_ONLY,
    .priv_size = sizeof(FPSSmoothContext),
    .init = fps_smooth_init,
    .uninit = fps_smooth_uninit,
    FILTER_INPUTS(fps_smooth_inputs),
    FILTER_OUTPUTS(fps_smooth_outputs),
    .activate = activate
};
