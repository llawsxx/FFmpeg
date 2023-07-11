#include "libavutil/common.h"
#include "libavutil/imgutils.h"
#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"
#include "avfilter.h"
#include "formats.h"
#include "internal.h"
#include "video.h"


typedef struct ADEINTContext {
    const AVClass *class;
    int parity;           ///< frame field parity
    int deint;
    int mode;
    double lthres1;
    double lthres2;
    double cthres1;
    double cthres2;
    int ela;
    
    int internal_lthres1;
    int internal_lthres2;
    int internal_cthres1;
    int internal_cthres2;
    
    int nb_threads;
    int eof;
    int (*deint_func)(AVFrame* prev, AVFrame* cur, AVFrame* next, AVFrame* out, int plane, int width, int height, int jobnr, int nb_jobs, int is_second, int field_order, int thres1, int thres2, int enable_ela);
    AVFrame *prev, *cur, *next;  ///< previous, current, next frames
    const AVPixFmtDescriptor *format_desc;
} ADEINTContext;

#define OFFSET(x) offsetof(ADEINTContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_RUNTIME_PARAM
#define CONST(name, help, val, unit) { name, help, 0, AV_OPT_TYPE_CONST, {.i64=val}, 0, 0, FLAGS, unit }

static const AVOption adeint_options[] = {
    { "parity", "specify the assumed picture field parity", OFFSET(parity), AV_OPT_TYPE_INT, {.i64=-1}, -1, 1, FLAGS, "parity" },
    CONST("tff",  "assume top field first",     0, "parity"),
    CONST("bff",  "assume bottom field first",  1, "parity"),
    CONST("auto", "auto detect parity",        -1, "parity"),
    { "deint",  "specify which frames to deinterlace", OFFSET(deint), AV_OPT_TYPE_INT, {.i64=0}, 0, 1, FLAGS, "deint" },
    CONST("all",        "deinterlace all frames",                       0, "deint"),
    CONST("interlaced", "only deinterlace frames marked as interlaced", 1, "deint"),
    { "mode",   "specify the interlacing mode", OFFSET(mode), AV_OPT_TYPE_INT, {.i64=1}, 0, 1, FLAGS, "mode"},
    CONST("frame", "send one frame for each frame", 0, "mode"),
    CONST("field", "send one frame for each field", 1, "mode"),
    { "lthres1",  "specify luma motion threshold starts to using intra deint", OFFSET(lthres1), AV_OPT_TYPE_DOUBLE, {.dbl=0.03}, 0, 1, FLAGS, "lthres1" },
    { "lthres2",  "specify luma motion threshold completely using intra deint", OFFSET(lthres2), AV_OPT_TYPE_DOUBLE, {.dbl=0.05}, 0, 1, FLAGS, "lthres2" },
    { "cthres1",  "specify chroma motion threshold starts to using intra deint", OFFSET(cthres1), AV_OPT_TYPE_DOUBLE, {.dbl=0.02}, 0, 1, FLAGS, "cthres1" },
    { "cthres2",  "specify chroma motion threshold completely using intra deint", OFFSET(cthres2), AV_OPT_TYPE_DOUBLE, {.dbl=0.03}, 0, 1, FLAGS, "cthres1" },
    { "ela",  "enable 5 tap edge line average", OFFSET(ela), AV_OPT_TYPE_INT, {.i64=1}, 0, 1, FLAGS, "ela" },

    { NULL }
};

AVFILTER_DEFINE_CLASS(adeint);

static const enum AVPixelFormat pix_fmts[] = {
    AV_PIX_FMT_YUV410P, AV_PIX_FMT_YUV411P,
    AV_PIX_FMT_YUV420P, AV_PIX_FMT_YUV422P,
    AV_PIX_FMT_YUV440P, AV_PIX_FMT_YUV444P,
    AV_PIX_FMT_YUVJ444P, AV_PIX_FMT_YUVJ440P,
    AV_PIX_FMT_YUVJ422P, AV_PIX_FMT_YUVJ420P,
    AV_PIX_FMT_YUVJ411P,
    AV_PIX_FMT_YUVA420P, AV_PIX_FMT_YUVA422P, AV_PIX_FMT_YUVA444P,
    AV_PIX_FMT_GBRP, AV_PIX_FMT_GBRAP,
    AV_PIX_FMT_GRAY8,
    AV_PIX_FMT_GRAY9, AV_PIX_FMT_GRAY10, AV_PIX_FMT_GRAY12, AV_PIX_FMT_GRAY14, AV_PIX_FMT_GRAY16,
    AV_PIX_FMT_YUV420P9, AV_PIX_FMT_YUV422P9, AV_PIX_FMT_YUV444P9,
    AV_PIX_FMT_YUV420P10, AV_PIX_FMT_YUV422P10, AV_PIX_FMT_YUV444P10,
    AV_PIX_FMT_YUV440P10,
    AV_PIX_FMT_YUV420P12, AV_PIX_FMT_YUV422P12, AV_PIX_FMT_YUV444P12,
    AV_PIX_FMT_YUV440P12,
    AV_PIX_FMT_YUV420P14, AV_PIX_FMT_YUV422P14, AV_PIX_FMT_YUV444P14,
    AV_PIX_FMT_YUV420P16, AV_PIX_FMT_YUV422P16, AV_PIX_FMT_YUV444P16,
    AV_PIX_FMT_GBRP9, AV_PIX_FMT_GBRP10, AV_PIX_FMT_GBRP12, AV_PIX_FMT_GBRP14, AV_PIX_FMT_GBRP16,
    AV_PIX_FMT_YUVA444P9, AV_PIX_FMT_YUVA444P10, AV_PIX_FMT_YUVA444P12, AV_PIX_FMT_YUVA444P16,
    AV_PIX_FMT_YUVA422P9, AV_PIX_FMT_YUVA422P10, AV_PIX_FMT_YUVA422P12, AV_PIX_FMT_YUVA422P16,
    AV_PIX_FMT_YUVA420P9, AV_PIX_FMT_YUVA420P10, AV_PIX_FMT_YUVA420P16,
    AV_PIX_FMT_GBRAP10,   AV_PIX_FMT_GBRAP12,    AV_PIX_FMT_GBRAP16,
    AV_PIX_FMT_NONE
};

#define ELA_TAP 5
#define B_IS_MEDIAN(a,b,c) (((a)<=(b) && (b)<=(c))||((a)>=(b) && (b)>=(c))) 

#define ELA_INTRA_DEINT(outname) \
if (enable_ela && x >= (ELA_TAP - 1) / 2 && x < width - (ELA_TAP - 1) / 2) {\
for (int i = 0; i < sizeof(scores) / sizeof(scores[0]); i++) {\
        int j = i - (ELA_TAP - 1) / 2;\
        scores[i] = FFABS((cur_p[-cur_linesize + (j)]) - cur_p[+cur_linesize - (j)]);\
};\
    spatial_score = spatial_max;\
for (int i = 0; i < sizeof(scores) / sizeof(scores[0]); i++) {\
    if (spatial_score > scores[i]) {\
        spatial_score = scores[i];\
        spatial_min_index = i;\
    }\
}\
{\
    int j = spatial_min_index - (ELA_TAP - 1) / 2;\
    int ela_point = ((cur_p[-cur_linesize + (j)]) + cur_p[+cur_linesize - (j)]) >> 1;\
    if (B_IS_MEDIAN(cur_p[-cur_linesize], ela_point, cur_p[+cur_linesize]))\
        outname = ela_point;\
    else\
        outname = ((cur_p[-cur_linesize]) + cur_p[+cur_linesize]) >> 1;\
}\
                }\
    else {\
        outname = (cur_p[-cur_linesize] + cur_p[+cur_linesize]) >> 1;\
}

#define GET_X_POINT_START(type) \
type prev_p, * cur_p, * next_p,*adj_p, * out_p;\
for (int x = 0; x < width; x++) {\
    prev_p = &((type)prev_line)[x];\
    cur_p = &((type)cur_line)[x];\
    next_p = &((type)next_line)[x];\
    adj_p = &((type)adj_line)[x];\
    out_p = &((type)out_line)[x];

#define GET_X_POINT_END(mul) \
}\
    y_out += 2;\
    prev_line += prev_linesize * 2 * mul;\
    cur_line += cur_linesize * 2 * mul;\
    next_line += next_linesize * 2 * mul;\
    adj_line += adj_linesize * 2 * mul;\
    out_line += out_linesize * 2 * mul;


static int adaptive_deint_plane(AVFrame* prev, AVFrame* cur, AVFrame* next, AVFrame* out,int plane,int width,int height, int jobnr, int nb_jobs, int is_second, int field_order, int thres1, int thres2,int enable_ela) {
    int is_first_in;
    int32_t m1, m2, m3, m4, diff1, diff2, diff3, max_diff;
    int32_t inter_deint_value, intra_deint_value;

    uint8_t* cur_line,*prev_line,*next_line,*adj_line,*out_line;
    
    int64_t prev_linesize, cur_linesize, next_linesize,adj_linesize, out_linesize;
    float factor;
    int y_out;
    const int start = (height * jobnr) / nb_jobs;
    const int end = (height * (jobnr + 1)) / nb_jobs;
    int32_t scores[ELA_TAP];
    int32_t spatial_max = INT16_MAX;
    int32_t spatial_score = spatial_max;
    int spatial_min_index;

    int interpolate_top_line = !!(field_order ^ is_second);

    prev_linesize = prev->linesize[plane];
    cur_linesize = cur->linesize[plane];
    next_linesize = next->linesize[plane];
    out_linesize = out->linesize[plane];
    adj_linesize = field_order != interpolate_top_line ? next_linesize : prev_linesize;

    y_out = start + (interpolate_top_line ^ (start & 1));
    cur_line = &((cur->data[plane])[y_out * cur_linesize]);
    out_line = &((out->data[plane])[y_out * out_linesize]);

    while (y_out < end) {
        memcpy(out_line, cur_line, cur_linesize > out_linesize ? out_linesize : cur_linesize);
        y_out += 2;
        cur_line += cur_linesize * 2;
        out_line += out_linesize * 2;
    }

    y_out = start + ((!interpolate_top_line) ^ (start & 1));
    prev_line = &((prev->data[plane])[y_out * prev_linesize]);
    cur_line = &((cur->data[plane])[y_out * cur_linesize]);
    next_line = &((next->data[plane])[y_out * next_linesize]);
    out_line = &((out->data[plane])[y_out * out_linesize]);
    adj_line = field_order != interpolate_top_line ? next_line : prev_line;

    if (y_out == 0) {
        GET_X_POINT_START(uint8_t*)

        m3 = FFABS(cur_p[+cur_linesize] - prev_p[+prev_linesize]);
        m4 = FFABS(cur_p[+cur_linesize] - next_p[+next_linesize]);
        diff1 = m3;
        diff2 = m4;
        diff3 = FFABS(*adj_p - *cur_p);
        max_diff = FFMAX3(diff1, diff2, diff3);

        if (max_diff < thres1)
            *out_p = (*adj_p + *cur_p) >> 1;
        else if (max_diff >= thres2)
        {
            *out_p = cur_p[+cur_linesize];
        }
        else
        {
            factor = (max_diff - thres1) / ((float)thres2 - thres1);
            inter_deint_value = (*adj_p + *cur_p) >> 1;
            intra_deint_value = cur_p[+cur_linesize];
            *out_p = (factor * intra_deint_value) + ((1 - factor) * inter_deint_value);
        }

        GET_X_POINT_END(sizeof(uint8_t));
    }


    while (y_out < (end >= height ? height - 1 : end)) {
        GET_X_POINT_START(uint8_t*)

        m1 = FFABS(cur_p[-cur_linesize] - prev_p[-prev_linesize]);
        m2 = FFABS(cur_p[-cur_linesize] - next_p[-next_linesize]);
        m3 = FFABS(cur_p[+cur_linesize] - prev_p[+prev_linesize]);
        m4 = FFABS(cur_p[+cur_linesize] - next_p[+next_linesize]);

        diff1 = (m1 + m3) >> 1;
        diff2 = (m2 + m4) >> 1;

        diff3 = FFABS(*adj_p - *cur_p);
        max_diff = FFMAX3(diff1, diff2, diff3);

        if (max_diff < thres1)
            *out_p = (*adj_p + *cur_p) >> 1;
        else if (max_diff >= thres2)
        {
            ELA_INTRA_DEINT(*out_p);
        }
        else
        {
            factor = (max_diff - thres1) / ((float)thres2 - thres1);
            inter_deint_value = (*adj_p + *cur_p) >> 1;

            ELA_INTRA_DEINT(intra_deint_value);

            *out_p = (factor * intra_deint_value) + ((1 - factor) * inter_deint_value);
        }

        GET_X_POINT_END(sizeof(uint8_t));
    }

    if (y_out == height - 1)
    {
        GET_X_POINT_START(uint8_t*)

        m1 = FFABS(cur_p[-cur_linesize] - prev_p[-prev_linesize]);
        m2 = FFABS(cur_p[-cur_linesize] - next_p[-next_linesize]);
        diff1 = m1;
        diff2 = m2;
        diff3 = FFABS(*adj_p - *cur_p);
        max_diff = FFMAX3(diff1, diff2, diff3);

        if (max_diff < thres1)
            *out_p = (*adj_p + *cur_p) >> 1;
        else if (max_diff >= thres2)
        {
            *out_p = cur_p[-cur_linesize];
        }
        else
        {
            factor = (max_diff - thres1) / ((float)thres2 - thres1);
            inter_deint_value = (*adj_p + *cur_p) >> 1;
            intra_deint_value = cur_p[-cur_linesize];
            *out_p = (factor * intra_deint_value) + ((1 - factor) * inter_deint_value);
        }

        GET_X_POINT_END(sizeof(uint8_t));
    }

    return 0;
}


static int adaptive_deint_plane_16bit(AVFrame* prev, AVFrame* cur, AVFrame* next, AVFrame* out, int plane, int width, int height, int jobnr, int nb_jobs, int is_second, int field_order, int thres1, int thres2, int enable_ela) {
    int is_first_in;
    int32_t m1, m2, m3, m4, diff1, diff2, diff3, max_diff;
    int32_t inter_deint_value, intra_deint_value;
    uint8_t* cur_line, * prev_line, * next_line, * adj_line, * out_line;
    int64_t prev_linesize, cur_linesize, next_linesize, adj_linesize, out_linesize;
    float factor;
    int y_out;
    const int start = (height * jobnr) / nb_jobs;
    const int end = (height * (jobnr + 1)) / nb_jobs;
    int32_t scores[ELA_TAP];
    int32_t spatial_max = INT16_MAX;
    int32_t spatial_score = spatial_max;
    int spatial_min_index;

    int interpolate_top_line = !!(field_order ^ is_second);

    prev_linesize = prev->linesize[plane];
    cur_linesize = cur->linesize[plane];
    next_linesize = next->linesize[plane];
    out_linesize = out->linesize[plane];
    adj_linesize = field_order != interpolate_top_line ? next_linesize : prev_linesize;

    y_out = start + (interpolate_top_line ^ (start & 1));
    cur_line = &((cur->data[plane])[y_out * cur_linesize]);
    out_line = &((out->data[plane])[y_out * out_linesize]);

    while (y_out < end) {
        memcpy(out_line, cur_line, cur_linesize > out_linesize ? out_linesize : cur_linesize);
        y_out += 2;
        cur_line += cur_linesize * 2;
        out_line += out_linesize * 2;
    }

    y_out = start + ((!interpolate_top_line) ^ (start & 1));
    prev_line = &((prev->data[plane])[y_out * prev_linesize]);
    cur_line = &((cur->data[plane])[y_out * cur_linesize]);
    next_line = &((next->data[plane])[y_out * next_linesize]);
    out_line = &((out->data[plane])[y_out * out_linesize]);
    adj_line = field_order != interpolate_top_line ? next_line : prev_line;

    prev_linesize /= 2;
    cur_linesize /= 2;
    next_linesize /= 2;
    adj_linesize /= 2;
    out_linesize /= 2;

    if (y_out == 0) {
        GET_X_POINT_START(uint16_t*)

        m3 = FFABS(cur_p[+cur_linesize] - prev_p[+prev_linesize]);
        m4 = FFABS(cur_p[+cur_linesize] - next_p[+next_linesize]);
        diff1 = m3;
        diff2 = m4;
        diff3 = FFABS(*adj_p - *cur_p);
        max_diff = FFMAX3(diff1, diff2, diff3);

        if (max_diff < thres1)
            *out_p = (*adj_p + *cur_p) >> 1;
        else if (max_diff >= thres2)
        {
            *out_p = cur_p[+cur_linesize];
        }
        else
        {
            factor = (max_diff - thres1) / ((float)thres2 - thres1);
            inter_deint_value = (*adj_p + *cur_p) >> 1;
            intra_deint_value = cur_p[+cur_linesize];
            *out_p = (factor * intra_deint_value) + ((1 - factor) * inter_deint_value);
        }

        GET_X_POINT_END(sizeof(uint16_t))
    }


    while (y_out < (end >= height ? height - 1 : end)) {
        GET_X_POINT_START(uint16_t*)


        m1 = FFABS(cur_p[-cur_linesize] - prev_p[-prev_linesize]);
        m2 = FFABS(cur_p[-cur_linesize] - next_p[-next_linesize]);
        m3 = FFABS(cur_p[+cur_linesize] - prev_p[+prev_linesize]);
        m4 = FFABS(cur_p[+cur_linesize] - next_p[+next_linesize]);

        diff1 = (m1 + m3) >> 1;
        diff2 = (m2 + m4) >> 1;

        diff3 = FFABS(*adj_p - *cur_p);
        max_diff = FFMAX3(diff1, diff2, diff3);

        if (max_diff < thres1)
            *out_p = (*adj_p + *cur_p) >> 1;
        else if (max_diff >= thres2)
        {
            ELA_INTRA_DEINT(*out_p);
        }
        else
        {
            factor = (max_diff - thres1) / ((float)thres2 - thres1);
            inter_deint_value = (*adj_p + *cur_p) >> 1;

            ELA_INTRA_DEINT(intra_deint_value);

            *out_p = (factor * intra_deint_value) + ((1 - factor) * inter_deint_value);
        }

        GET_X_POINT_END(sizeof(uint16_t))
    }

    if (y_out == height - 1)
    {
        GET_X_POINT_START(uint16_t*)

        m1 = FFABS(cur_p[-cur_linesize] - prev_p[-prev_linesize]);
        m2 = FFABS(cur_p[-cur_linesize] - next_p[-next_linesize]);
        diff1 = m1;
        diff2 = m2;
        diff3 = FFABS(*adj_p - *cur_p);
        max_diff = FFMAX3(diff1, diff2, diff3);

        if (max_diff < thres1)
            *out_p = (*adj_p + *cur_p) >> 1;
        else if (max_diff >= thres2)
        {
            *out_p = cur_p[-cur_linesize];
        }
        else
        {
            factor = (max_diff - thres1) / ((float)thres2 - thres1);
            inter_deint_value = (*adj_p + *cur_p) >> 1;
            intra_deint_value = cur_p[-cur_linesize];
            *out_p = (factor * intra_deint_value) + ((1 - factor) * inter_deint_value);
        }

        GET_X_POINT_END(sizeof(uint16_t))
    }

    return 0;
}

static int config_input(AVFilterLink *inlink)
{
    AVFilterContext *ctx = inlink->dst;
    ADEINTContext *s = ctx->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);
    int depth,depth_max;

    if(s->lthres2 < s->lthres1){
        av_log(ctx, AV_LOG_ERROR, "Failed to process command. lthres2 must greater than or equal to lthres1\n");
        return AVERROR(EINVAL);
    }
    if(s->cthres2 < s->cthres1){
        av_log(ctx, AV_LOG_ERROR, "Failed to process command. cthres2 must greater than or equal to cthres1\n");
        return AVERROR(EINVAL);
    }

    if (inlink->h % 2 == 1) {
        av_log(ctx, AV_LOG_ERROR, "Video of odd lines is not supported\n");
        return AVERROR(EINVAL);
    }

    s->nb_threads = ff_filter_get_nb_threads(ctx);

    depth = desc->comp[0].depth;
    if (depth <= 8) {
        s->deint_func = adaptive_deint_plane;
    } else {
        s->deint_func = adaptive_deint_plane_16bit;
    }
    depth_max = (1 << depth) - 1;
    s->internal_lthres1 = s->lthres1 * depth_max;
    s->internal_lthres2 = s->lthres2 * depth_max;
    s->internal_cthres1 = s->cthres1 * depth_max;
    s->internal_cthres2 = s->cthres2 * depth_max;
    
    s->format_desc = desc;
    
    av_log(ctx, AV_LOG_DEBUG, "Deinterlace thresholds(for %d bit depth): lthres1:%d lthres2:%d cthres1:%d cthres2:%d ela:%d\n",
           depth,
           s->internal_lthres1,
           s->internal_lthres2,
           s->internal_cthres1,
           s->internal_cthres2,
           s->ela
        );

    return 0;
}

static int config_output(AVFilterLink *outlink)
{
    AVFilterLink *inlink = outlink->src->inputs[0];
    AVFilterContext *ctx = outlink->src;
    ADEINTContext *s = ctx->priv;

    outlink->time_base.num = inlink->time_base.num;
    outlink->time_base.den = inlink->time_base.den * 2;
    if (s->mode == 1)
        outlink->frame_rate.num = inlink->frame_rate.num * 2;
    else
        outlink->frame_rate.num = inlink->frame_rate.num;
    outlink->frame_rate.den = inlink->frame_rate.den;

    return 0;
}

typedef struct ThreadData {
    AVFrame *out;
    int is_second;
} ThreadData;


static int deinterlace_slice(AVFilterContext *ctx, void *arg,
                             int jobnr, int nb_jobs)
{
    ADEINTContext *s = ctx->priv;
    ThreadData *td = arg;    
    AVFrame *out = td->out;
    int is_second = td->is_second;
    int plane;

    int interlaced = !!(s->cur->flags & AV_FRAME_FLAG_INTERLACED);
    int field_order = s->parity == -1 ? interlaced ? !(s->cur->flags & AV_FRAME_FLAG_TOP_FIELD_FIRST) : 0 :
                                 s->parity;
    
    for (plane = 0; plane < s->format_desc->nb_components; plane++) {
        int width, height, thres1, thres2;
        if (plane == 0 || plane == 3) {
            width = out->width;
            height = out->height;
            thres1 = s->internal_lthres1;
            thres2 = s->internal_lthres2;
        }
        else {
            width = AV_CEIL_RSHIFT(out->width, s->format_desc->log2_chroma_w);
            height = AV_CEIL_RSHIFT(out->height, s->format_desc->log2_chroma_h);
            thres1 = s->internal_cthres1;
            thres2 = s->internal_cthres2;
        }
        s->deint_func(s->prev, s->cur, s->next, out, plane, width, height, jobnr, nb_jobs, is_second, field_order, thres1, thres2, s->ela);
    }
    
    return 0;
}


static int filter(AVFilterContext *ctx, int is_second)
{
    ADEINTContext *s = ctx->priv;
    AVFilterLink *outlink = ctx->outputs[0];
    AVFrame *out;
    ThreadData td;

    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!out)
        return AVERROR(ENOMEM);
    av_frame_copy_props(out, s->cur);
    
#if FF_API_INTERLACED_FRAME
FF_DISABLE_DEPRECATION_WARNINGS
    out->interlaced_frame = 0;
FF_ENABLE_DEPRECATION_WARNINGS
#endif
    out->flags &= ~AV_FRAME_FLAG_INTERLACED;

    if (!is_second) {
        if (out->pts != AV_NOPTS_VALUE)
            out->pts *= 2;
    } else {
        int64_t cur_pts  = s->cur->pts;
        int64_t next_pts = s->next->pts;

        if (next_pts != AV_NOPTS_VALUE && cur_pts != AV_NOPTS_VALUE) {
            out->pts = cur_pts + next_pts;
        } else {
            out->pts = AV_NOPTS_VALUE;
        }
    }

    td.out = out;
    td.is_second = is_second;

    ctx->internal->execute(ctx, deinterlace_slice, &td, NULL, FFMIN(out->height, s->nb_threads));

    return ff_filter_frame(outlink, out);
}


static int filter_frame(AVFilterLink *inlink, AVFrame *frame)
{
    AVFilterContext *ctx = inlink->dst;
    ADEINTContext *s = ctx->priv;
    int ret;

    av_frame_free(&s->prev);
    s->prev = s->cur;
    s->cur  = s->next;
    s->next = frame;

    if (!s->cur) {
        s->cur = av_frame_clone(s->next);
        if (!s->cur)
            return AVERROR(ENOMEM);
    }

    if ((s->deint && !(s->cur->flags & AV_FRAME_FLAG_INTERLACED)) || ctx->is_disabled) {
        AVFrame *out = av_frame_clone(s->cur);
        if (!out)
            return AVERROR(ENOMEM);

        av_frame_free(&s->prev);
        if (out->pts != AV_NOPTS_VALUE)
            out->pts *= 2;
        return ff_filter_frame(ctx->outputs[0], out);
    }

    if (!s->prev)
        return 0;


    ret = filter(ctx, 0);
    if (ret < 0 || s->mode == 0)
        return ret;

    return filter(ctx, 1);
}


static int request_frame(AVFilterLink *outlink)
{
    AVFilterContext *ctx = outlink->src;
    ADEINTContext *s = ctx->priv;
    int ret;

    if (s->eof)
        return AVERROR_EOF;

    ret = ff_request_frame(ctx->inputs[0]);

    if (ret == AVERROR_EOF && s->cur) {
        AVFrame *next = av_frame_clone(s->next);
        if (!next)
            return AVERROR(ENOMEM);
        next->pts = s->next->pts * 2 - s->cur->pts;
        filter_frame(ctx->inputs[0], next);
        s->eof = 1;
    } else if (ret < 0) {
        return ret;
    }

    return 0;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    ADEINTContext *s = ctx->priv;

    av_frame_free(&s->prev);
    av_frame_free(&s->cur );
    av_frame_free(&s->next);
}

static const AVFilterPad adeint_inputs[] = {
    {
        .name          = "default",
        .type          = AVMEDIA_TYPE_VIDEO,
        .filter_frame  = filter_frame,
        .config_props  = config_input,
    },
};

static const AVFilterPad adeint_outputs[] = {
    {
        .name          = "default",
        .type          = AVMEDIA_TYPE_VIDEO,
        .config_props  = config_output,
        .request_frame = request_frame,
    },
};

const AVFilter ff_vf_adeint = {
    .name          = "adeint",
    .description   = NULL_IF_CONFIG_SMALL("Apply simple 5 fields motion adaptive deinterlace."),
    .priv_size     = sizeof(ADEINTContext),
    .priv_class    = &adeint_class,
    .uninit        = uninit,
    FILTER_PIXFMTS_ARRAY(pix_fmts),
    FILTER_INPUTS(adeint_inputs),
    FILTER_OUTPUTS(adeint_outputs),
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_INTERNAL|AVFILTER_FLAG_SLICE_THREADS,
    .process_command = ff_filter_process_command,
};

