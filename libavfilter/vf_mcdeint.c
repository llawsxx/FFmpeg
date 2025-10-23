/*
 * Copyright (c) 2006 Michael Niedermayer <michaelni@gmx.at>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with FFmpeg; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/**
 * @file
 * Motion Compensation Deinterlacer
 * Ported from MPlayer libmpcodecs/vf_mcdeint.c.
 *
 * Known Issues:
 *
 * The motion estimation is somewhat at the mercy of the input, if the
 * input frames are created purely based on spatial interpolation then
 * for example a thin black line or another random and not
 * interpolateable pattern will cause problems.
 * Note: completely ignoring the "unavailable" lines during motion
 * estimation did not look any better, so the most obvious solution
 * would be to improve tfields or penalize problematic motion vectors.
 *
 * If non iterative ME is used then snow currently ignores the OBMC
 * window and as a result sometimes creates artifacts.
 *
 * Only past frames are used, we should ideally use future frames too,
 * something like filtering the whole movie in forward and then
 * backward direction seems like an interesting idea but the current
 * filter framework is FAR from supporting such things.
 *
 * Combining the motion compensated image with the input image also is
 * not as trivial as it seems, simple blindly taking even lines from
 * one and odd ones from the other does not work at all as ME/MC
 * sometimes has nothing in the previous frames which matches the
 * current. The current algorithm has been found by trial and error
 * and almost certainly can be improved...
 */

#include "libavutil/opt.h"
#include "libavutil/mem.h"
#include "libavcodec/avcodec.h"
#include "libavutil/pixdesc.h"
#include "libavcodec/motion_est.h"
#include "avfilter.h"
#include "filters.h"
#include "video.h"

#define FF_ME_ITER 3

enum MCDeintParity {
    PARITY_TFF  =  0, ///< top field first
    PARITY_BFF  =  1, ///< bottom field first
};

typedef struct MCDeintContext {
    const AVClass *class;
    int parity;         ///< MCDeintParity
    int pixel_ajust_mode;
    int motion_est;
    int iterative_dia_size;
    int ref_frames;
    int four_mv;
    int qpel;
    int output_memc;
    int qp;
    int log2_chroma_w;
    int log2_chroma_h;
    int nb_planes;
    int32_t *last_diff;
    AVPacket *pkt;
    AVFrame *frame_dec;
    AVCodecContext *enc_ctx;
} MCDeintContext;

#define OFFSET(x) offsetof(MCDeintContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM
#define CONST(name, help, val, u) { name, help, 0, AV_OPT_TYPE_CONST, {.i64=val}, INT_MIN, INT_MAX, FLAGS, .unit = u }

static const AVOption mcdeint_options[] = {
    { "parity", "set the assumed picture field parity", OFFSET(parity), AV_OPT_TYPE_INT, {.i64=PARITY_TFF}, PARITY_TFF, PARITY_BFF, FLAGS, "parity" },
    CONST("tff", "assume top field first",    PARITY_TFF, "parity"),
    CONST("bff", "assume bottom field first", PARITY_BFF, "parity"),

    { "qp", "set qp", OFFSET(qp), AV_OPT_TYPE_INT, {.i64 = 1}, INT_MIN, INT_MAX, FLAGS },

    { "motion_est", "set motion estimation algorithm", OFFSET(motion_est), AV_OPT_TYPE_INT, {.i64 = FF_ME_EPZS}, FF_ME_ZERO, FF_ME_ITER, FLAGS, .unit = "motion_est" },
    CONST("zero",       NULL, FF_ME_ZERO, "motion_est"),
    CONST("epzs",       NULL, FF_ME_EPZS, "motion_est"),
    CONST("xone",       NULL, FF_ME_XONE, "motion_est"),
    CONST("iter",       NULL, FF_ME_ITER, "motion_est"),

    { "dia_size", "dia size for the iterative ME", OFFSET(iterative_dia_size), AV_OPT_TYPE_INT, { .i64 = 0 }, 0, INT_MAX, FLAGS },
    { "refs", "set ref frames count", OFFSET(ref_frames), AV_OPT_TYPE_INT, { .i64 = 1 }, 1, 8, FLAGS },
    { "four_mv", "4 MV per MB allowed / advanced prediction for H.263.", OFFSET(four_mv), AV_OPT_TYPE_INT, { .i64 = 0 }, 0, 1, FLAGS },
    { "qpel", "use qpel motion compensation", OFFSET(qpel), AV_OPT_TYPE_INT, { .i64 = 1 }, 0, 1, FLAGS },
    { "adj", "pixel ajust mode", OFFSET(pixel_ajust_mode), AV_OPT_TYPE_INT, { .i64 = 0 }, 0, 1, FLAGS ,"pixel_adjust"},
    CONST("minor",       NULL, 0, "pixel_adjust"),
    CONST("normal",      NULL, 1, "pixel_adjust"),
    { "output_memc", "output memc image", OFFSET(output_memc), AV_OPT_TYPE_INT, { .i64 = 0 }, 0, 1, FLAGS },

    { NULL }
};

AVFILTER_DEFINE_CLASS(mcdeint);

static int config_props(AVFilterLink *inlink)
{
    AVFilterContext *ctx = inlink->dst;
    MCDeintContext *mcdeint = ctx->priv;
    const AVCodec *enc;
    AVCodecContext *enc_ctx;
    AVDictionary *opts = NULL;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);
    int ret;

    if((mcdeint->last_diff = av_mallocz(sizeof(*mcdeint->last_diff) * inlink->w)) == NULL)
        return AVERROR(ENOMEM);

    if (!(enc = avcodec_find_encoder(AV_CODEC_ID_SNOW))) {
        av_log(ctx, AV_LOG_ERROR, "Snow encoder is not enabled in libavcodec\n");
        return AVERROR(EINVAL);
    }

    mcdeint->log2_chroma_h = desc->log2_chroma_h;
    mcdeint->log2_chroma_w = desc->log2_chroma_w;
    mcdeint->nb_planes = av_pix_fmt_count_planes(inlink->format);

    mcdeint->pkt = av_packet_alloc();
    if (!mcdeint->pkt)
    return AVERROR(ENOMEM);
    mcdeint->frame_dec = av_frame_alloc();
    if (!mcdeint->frame_dec)
        return AVERROR(ENOMEM);
    mcdeint->enc_ctx = avcodec_alloc_context3(enc);
    if (!mcdeint->enc_ctx)
        return AVERROR(ENOMEM);
    enc_ctx = mcdeint->enc_ctx;
    enc_ctx->width  = inlink->w;
    enc_ctx->height = inlink->h;
    enc_ctx->time_base = (AVRational){1,25};  // meaningless
    enc_ctx->gop_size = INT_MAX;
    enc_ctx->max_b_frames = 0;
    enc_ctx->pix_fmt = inlink->format;
    enc_ctx->flags = AV_CODEC_FLAG_QSCALE | AV_CODEC_FLAG_LOW_DELAY | AV_CODEC_FLAG_RECON_FRAME;
    enc_ctx->strict_std_compliance = FF_COMPLIANCE_EXPERIMENTAL;
    enc_ctx->global_quality = 1;
    enc_ctx->me_cmp = enc_ctx->me_sub_cmp = FF_CMP_SAD;
    enc_ctx->mb_cmp = FF_CMP_SSE;
    av_dict_set(&opts, "memc_only", "1", 0);
    av_dict_set(&opts, "no_bitstream", "1", 0);

    switch (mcdeint->motion_est) {
    case FF_ME_ZERO:
        av_dict_set(&opts, "motion_est", "zero", 0);
        break;
    case FF_ME_EPZS:
        av_dict_set(&opts, "motion_est", "epzs", 0);
        break;
    case FF_ME_XONE:
        av_dict_set(&opts, "motion_est", "xone", 0);
        break;
    case FF_ME_ITER:
        av_dict_set(&opts, "motion_est", "iter", 0);
        break;
    }

    enc_ctx->dia_size = mcdeint->iterative_dia_size;
    enc_ctx->refs = mcdeint->ref_frames;

    if(mcdeint->four_mv)
        enc_ctx->flags |= AV_CODEC_FLAG_4MV;

    if(mcdeint->qpel)
        enc_ctx->flags |= AV_CODEC_FLAG_QPEL;

    ret = avcodec_open2(enc_ctx, enc, &opts);
    av_dict_free(&opts);
    if (ret < 0)
        return ret;

    return 0;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    MCDeintContext *mcdeint = ctx->priv;

    av_free(mcdeint->last_diff);
    av_packet_free(&mcdeint->pkt);
    avcodec_free_context(&mcdeint->enc_ctx);
    av_frame_free(&mcdeint->frame_dec);
}



static int filter_frame(AVFilterLink* inlink, AVFrame* inpic)
{
    MCDeintContext *mcdeint = inlink->dst->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    AVFrame *outpic, *frame_dec = mcdeint->frame_dec;
    AVPacket *pkt = mcdeint->pkt;
    int32_t *last_diff = mcdeint->last_diff;
    int x, y, i, ret;

    outpic = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!outpic) {
        av_frame_free(&inpic);
        return AVERROR(ENOMEM);
    }
    av_frame_copy_props(outpic, inpic);
    inpic->quality = mcdeint->qp * FF_QP2LAMBDA;

    ret = avcodec_send_frame(mcdeint->enc_ctx, inpic);
    if (ret < 0) {
        av_log(mcdeint->enc_ctx, AV_LOG_ERROR, "Error sending a frame for encoding\n");
        goto end;
    }
    ret = avcodec_receive_packet(mcdeint->enc_ctx, pkt);
    if (ret < 0) {
        av_log(mcdeint->enc_ctx, AV_LOG_ERROR, "Error receiving a packet from encoding\n");
        goto end;
    }
    av_packet_unref(pkt);
    ret = avcodec_receive_frame(mcdeint->enc_ctx, frame_dec);
    if (ret < 0) {
        av_log(mcdeint->enc_ctx, AV_LOG_ERROR, "Error receiving a frame from encoding\n");
        goto end;
    }
    
    if(mcdeint->output_memc)
        av_frame_copy(outpic, frame_dec);

    for (i = 0; i < mcdeint->nb_planes; i++) {
        int w = inlink->w;
        int h = inlink->h;
        int fils = frame_dec->linesize[i];
        int srcs = inpic->linesize[i];
        int dsts = outpic->linesize[i];
        int start_line;
        
        if(i > 0){
            w = AV_CEIL_RSHIFT(w, mcdeint->log2_chroma_w);
            h = AV_CEIL_RSHIFT(h, mcdeint->log2_chroma_h);
        }

        if (mcdeint->parity == 0) {
            start_line = 1;
            if(!mcdeint->output_memc)
                memcpy(&outpic->data[i][(h - 1) * dsts], &frame_dec->data[i][(h - 1) * fils], sizeof(uint8_t) * w);
            for (x = 0; x < w; x++) {
                uint8_t* filp = &frame_dec->data[i][x];
                uint8_t* srcp = &inpic->data[i][x];
                last_diff[x] = filp[0] - srcp[0];
            }
        }
        else {
            start_line = 2;
            if(!mcdeint->output_memc)
                memcpy(&outpic->data[i][0], &frame_dec->data[i][0], sizeof(uint8_t) * w);
            for (x = 0; x < w; x++) {
                uint8_t* filp = &frame_dec->data[i][x + fils];
                uint8_t* srcp = &inpic->data[i][x + srcs];
                last_diff[x] = filp[0] - srcp[0];
            }
        }

        for (y = start_line; y < h - 1; y+=2) {
            for (x = 0; x < w; x++) {
                    uint8_t* filp = &frame_dec->data[i][x + y * fils];
                    uint8_t* srcp = &inpic->data[i][x + y * srcs];
                    uint8_t* dstp = &outpic->data[i][x + y * dsts];
                    int diff0 = last_diff[x];
                    int diff1 = last_diff[x] = filp[+fils] - srcp[+srcs];
                    int temp = filp[0];
                if(mcdeint->pixel_ajust_mode){
                    temp -= (diff0 + diff1) / 2;
                }else{
                    if (diff0 + diff1 > 0)
                        temp -= (diff0 + diff1 - FFABS(FFABS(diff0) - FFABS(diff1)) / 2) / 2;
                    else
                        temp -= (diff0 + diff1 + FFABS(FFABS(diff0) - FFABS(diff1)) / 2) / 2;
                }
                *filp = temp > 255U ? ~(temp >> 31) : temp;
                if(!mcdeint->output_memc)
                    *dstp = *filp;
            }
        }

        for (y = start_line - 1; y < h; y+=2) {
            memcpy(&frame_dec->data[i][y * fils], &inpic->data[i][y * srcs],sizeof(uint8_t) * w);
            if(!mcdeint->output_memc)
                memcpy(&outpic->data[i][y * dsts], &inpic->data[i][y * srcs], sizeof(uint8_t) * w);
        }
    }
    mcdeint->parity ^= 1;

end:
    av_packet_unref(pkt);
    av_frame_free(&inpic);
    if (ret < 0) {
        av_frame_free(&outpic);
        return ret;
    }
    return ff_filter_frame(outlink, outpic);
}

static const AVFilterPad mcdeint_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .filter_frame = filter_frame,
        .config_props = config_props,
    },
};

const FFFilter ff_vf_mcdeint = {
    .p.name        = "mcdeint",
    .p.description = NULL_IF_CONFIG_SMALL("Apply motion compensating deinterlacing."),
    .p.priv_class  = &mcdeint_class,
    .priv_size     = sizeof(MCDeintContext),
    .uninit        = uninit,
    FILTER_INPUTS(mcdeint_inputs),
    FILTER_OUTPUTS(ff_video_default_filterpad),
    FILTER_PIXFMTS(AV_PIX_FMT_YUV420P, AV_PIX_FMT_YUV444P, AV_PIX_FMT_YUV410P, AV_PIX_FMT_GRAY8),
};
