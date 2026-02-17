/*
 * Adaptive Bitrate FIFO pseudo-muxer
 * Copyright (c) 2016 Jan Sebechlebsky
 * Copyright (c) 2024 Adaptive Bitrate Extension
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with FFmpeg; if not, write to the Free Software * Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include <stdatomic.h>

#include "libavutil/avassert.h"
#include "libavutil/opt.h"
#include "libavutil/time.h"
#include "libavutil/thread.h"
#include "libavutil/threadmessage.h"
#include "libavutil/pixdesc.h"
#include "libavutil/parseutils.h"
#include "libavutil/eval.h"
#include "libavutil/avstring.h"
#include "libavutil/mem.h"
#include "libavcodec/avcodec.h"
#include "libswscale/swscale.h"
#include "avformat.h"
#include "internal.h"
#include "mux.h"

#define FIFO_DEFAULT_QUEUE_SIZE              300
#define FIFO_DEFAULT_MAX_RECOVERY_ATTEMPTS   0
#define FIFO_DEFAULT_RECOVERY_WAIT_TIME_USEC 5000000 // 5 seconds
#define MAX_BITRATE_LEVELS                   8

typedef struct BitrateLevel {
    int64_t bitrate;        // Target bitrate in bits per second
} BitrateLevel;

typedef struct FifoContext {
    const AVClass* class;
    AVFormatContext* avf;

    char* format;
    AVDictionary* format_options;

    int queue_size;
    AVThreadMessageQueue* queue;

    pthread_t writer_thread;

    /* Return value of last write_trailer_call */
    int write_trailer_ret;

    /* Time to wait before next recovery attempt */
    int64_t recovery_wait_time;

    /* Maximal number of unsuccessful successive recovery attempts */
    int max_recovery_attempts;

    /* Whether to attempt recovery from failure */
    int attempt_recovery;

    /* If >0 stream time will be used when waiting for recovery */
    int recovery_wait_streamtime;

    /* If >0 recovery will be attempted regardless of error code */
    int recover_any_error;

    /* Whether to drop packets in case the queue is full */
    int drop_pkts_on_overflow;

    /* Whether to wait for keyframe when recovering */
    int restart_with_keyframe;

    pthread_mutex_t overflow_flag_lock;
    int overflow_flag_lock_initialized;
    volatile uint8_t overflow_flag;

    /* Adaptive bitrate control */
    int enable_abr;                     // Enable adaptive bitrate encoding
    char* abr_codec;                    // Encoder name (e.g., "libx264")
    AVDictionary* abr_codec_options;    // Options for encoder
    char* abr_bitrates;                 // Comma-separated bitrates: "500k,1M,2M"
    int abr_levels;                     // Number of bitrate levels parsed

    /* ABR switching parameters */
    int abr_low_threshold;              // Queue percentage to switch up (default 25)
    int abr_high_threshold;             // Queue percentage to switch down (default 75)
    int64_t abr_down_switch_interval;        // Min time between level switches (ms)
    int64_t abr_up_switch_interval;        // Min time between level switches (ms)

    BitrateLevel bitrate_levels[MAX_BITRATE_LEVELS]; // Bitrate configurations

    /* Video stream index for ABR */
    int video_stream_idx;

    /* Encoder contexts for each bitrate level (pre-initialized) */
    AVCodecContext* enc_ctx;
    int current_bitrate_idx;            // Currently active bitrate index

    /* ABR state tracking */
    int64_t last_switch_time;           // Timestamp of last bitrate switch
    int last_qsize;                     // Last queue size
    double qsize_since_switch_trend;    // Queue size trend

    /* Decoder context for input video */
    AVCodecContext* dec_ctx;

    /* Current encoding parameters */
    int enc_width;
    int enc_height;
    enum AVPixelFormat enc_pix_fmt;
    AVRational enc_framerate;
    AVRational enc_time_base;
    AVRational enc_sample_aspect_ratio;

    enum AVColorPrimaries enc_color_primaries;
    enum AVColorTransferCharacteristic enc_color_trc;
    enum AVColorSpace enc_colorspace;
    enum AVColorRange enc_color_range;
    enum AVChromaLocation enc_chroma_sample_location;
    enum AVFieldOrder enc_field_order;

    /* Frame buffer for decoding/encoding */
    AVFrame* decoded_frame;
    AVPacket* enc_pkt;
    AVPacket* pkt_for_write;
    const AVOutputFormat* oformat;
} FifoContext;

typedef struct FifoThreadContext {
    AVFormatContext* avf;

    /* Timestamp of last failure */
    int64_t last_recovery_ts;

    /* Number of current recovery process */
    int recovery_nr;

    /* If > 0 all frames will be dropped until keyframe is received */
    uint8_t drop_until_keyframe;

    /* Value > 0 means previous write_header was successful */
    uint8_t header_written;

    int64_t last_received_dts;
} FifoThreadContext;

typedef enum FifoMessageType {
    FIFO_NOOP,
    FIFO_WRITE_HEADER,
    FIFO_WRITE_PACKET,
    FIFO_FLUSH_OUTPUT
} FifoMessageType;

typedef struct FifoMessage {
    FifoMessageType type;
    AVPacket pkt;
} FifoMessage;

/* Parse bitrate string with k/M/G suffix (e.g., "500k", "1M", "2.5M") */
static int64_t parse_bitrate_string(const char* str)
{
    char* end;
    double val = av_strtod(str, &end);

    if (end == str) {
        return -1;
    }

    if (*end) {
        switch (av_tolower(*end)) {
        case 'k':
            val *= 1000;
            break;
        case 'm':
            val *= 1000000;
            break;
        case 'g':
            val *= 1000000000;
            break;
        default:
            return -1;
        }
    }

    return (int64_t)val;
}

/* Parse comma-separated bitrate string into levels */
static int parse_abr_bitrates(FifoContext* fifo)
{
    if (!fifo->abr_bitrates || !*fifo->abr_bitrates) {
        av_log(fifo, AV_LOG_ERROR, "ABR enabled but no bitrates specified\n");
        return AVERROR(EINVAL);
    }

    char* str = av_strdup(fifo->abr_bitrates);
    if (!str)
        return AVERROR(ENOMEM);

    char* saveptr;
    char* token = av_strtok(str, ",", &saveptr);
    int count = 0;

    while (token && count < MAX_BITRATE_LEVELS) {
        // Trim whitespace
        while (av_isspace(*token)) token++;
        char* end = token + strlen(token) - 1;
        while (end > token && av_isspace(*end)) *end-- = '\0';

        int64_t bitrate = parse_bitrate_string(token);
        if (bitrate < 0) {
            av_log(fifo, AV_LOG_ERROR, "Invalid bitrate format: %s\n", token);
            av_free(str);
            return AVERROR(EINVAL);
        }

        fifo->bitrate_levels[count].bitrate = bitrate;
        av_log(fifo, AV_LOG_VERBOSE, "ABR level %d: bitrate=%" PRId64 " bps\n",
            count, bitrate);

        count++;
        token = av_strtok(NULL, ",", &saveptr);
    }

    av_free(str);

    if (count < 2) {
        av_log(fifo, AV_LOG_ERROR, "At least 2 bitrate levels required for ABR\n");
        return AVERROR(EINVAL);
    }

    if (count > MAX_BITRATE_LEVELS) {
        av_log(fifo, AV_LOG_WARNING, "Too many bitrate levels, using first %d\n",
            MAX_BITRATE_LEVELS);
        count = MAX_BITRATE_LEVELS;
    }

    fifo->abr_levels = count;

    /* Sort bitrates in ascending order (lowest first) for easier indexing */
    /* Level 0 = lowest bitrate (used when queue is full), Level N = highest (queue empty) */

    return 0;
}

/* Initialize decoder for the first video stream */
static int init_decoder(FifoContext* fifo, AVStream* st)
{
    const AVCodec* dec = avcodec_find_decoder(st->codecpar->codec_id);
    if (!dec) {
        av_log(fifo, AV_LOG_ERROR, "Failed to find decoder for stream %d\n", st->index);
        return AVERROR_DECODER_NOT_FOUND;
    }

    fifo->dec_ctx = avcodec_alloc_context3(dec);
    if (!fifo->dec_ctx) {
        return AVERROR(ENOMEM);
    }

    int ret = avcodec_parameters_to_context(fifo->dec_ctx, st->codecpar);
    if (ret < 0) {
        avcodec_free_context(&fifo->dec_ctx);
        return ret;
    }

    ret = avcodec_open2(fifo->dec_ctx, dec, NULL);
    if (ret < 0) {
        av_log(fifo, AV_LOG_ERROR, "Failed to open decoder: %s\n", av_err2str(ret));
        avcodec_free_context(&fifo->dec_ctx);
        return ret;
    }

    return 0;
}

static int init_encoder(FifoContext* fifo, int i, const AVOutputFormat* oformat) {
    const AVCodec* enc = avcodec_find_encoder_by_name(fifo->abr_codec);
    if (!enc) {
        av_log(fifo, AV_LOG_ERROR, "Encoder '%s' not found\n", fifo->abr_codec);
        return AVERROR_ENCODER_NOT_FOUND;
    }
    if (fifo->enc_ctx) {
        avcodec_free_context(&fifo->enc_ctx);
    }

    fifo->enc_ctx = avcodec_alloc_context3(enc);
    if (!fifo->enc_ctx) {
        return AVERROR(ENOMEM);
    }

    /* Set encoding parameters */
    fifo->enc_ctx->width = fifo->enc_width;
    fifo->enc_ctx->height = fifo->enc_height;
    fifo->enc_ctx->pix_fmt = fifo->enc_pix_fmt;
    fifo->enc_ctx->sample_aspect_ratio = fifo->enc_sample_aspect_ratio;
    fifo->enc_ctx->time_base = fifo->enc_time_base;
    fifo->enc_ctx->framerate = fifo->enc_framerate;
    fifo->enc_ctx->field_order = fifo->enc_field_order;
    fifo->enc_ctx->color_range = fifo->enc_color_range;
    fifo->enc_ctx->color_primaries = fifo->enc_color_primaries;
    fifo->enc_ctx->color_trc = fifo->enc_color_trc;
    fifo->enc_ctx->colorspace = fifo->enc_colorspace;
    fifo->enc_ctx->chroma_sample_location = fifo->enc_chroma_sample_location;

    if (oformat && oformat->flags & AVFMT_GLOBALHEADER)
        fifo->enc_ctx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

    /* Validate pixel format support */
    if (fifo->enc_pix_fmt == AV_PIX_FMT_NONE) {
        av_log(fifo, AV_LOG_ERROR, "Unknown pixel format for video stream\n");
        return AVERROR_INVALIDDATA;
    }

    fifo->enc_ctx->bit_rate = fifo->bitrate_levels[i].bitrate;
    fifo->enc_ctx->rc_min_rate = fifo->bitrate_levels[i].bitrate * 0.8;
    fifo->enc_ctx->rc_max_rate = fifo->bitrate_levels[i].bitrate * 1.2;
    fifo->enc_ctx->rc_buffer_size = (int)(fifo->bitrate_levels[i].bitrate * 2);

    /* Apply user options */
    AVDictionary* opts = NULL;
    if (fifo->abr_codec_options) {
        av_dict_copy(&opts, fifo->abr_codec_options, 0);
    }

    int ret = avcodec_open2(fifo->enc_ctx, enc, &opts);
    av_dict_free(&opts);

    if (ret < 0) {
        av_log(fifo, AV_LOG_ERROR, "Failed to open encoder at bitrate %" PRId64 ": %s\n",
            fifo->bitrate_levels[i].bitrate, av_err2str(ret));
        return ret;
    }

    av_log(fifo, AV_LOG_INFO, "Initialized encoder level %d: bitrate=%" PRId64 " bps\n",
        i, fifo->bitrate_levels[i].bitrate);
    return 0;
}

static int init_decoder_and_encoder(FifoContext* fifo, AVFormatContext* in_avf, const AVOutputFormat* oformat)
{
    int ret;
    /* Get source stream parameters */
    AVStream* src_st = NULL;

    for (unsigned int i = 0; i < in_avf->nb_streams; i++) {
        if (in_avf->streams[i]->codecpar->codec_type == AVMEDIA_TYPE_VIDEO) {
            src_st = in_avf->streams[i];
            fifo->video_stream_idx = i;
            ret = init_decoder(fifo, in_avf->streams[i]);
            if (ret < 0)
                return ret;
            break;
        }
    }

    if (!src_st) {
        av_log(fifo, AV_LOG_ERROR, "No video stream found for ABR\n");
        return AVERROR_INVALIDDATA;
    }

    /* Store encoding parameters */
    fifo->enc_width = src_st->codecpar->width;
    fifo->enc_height = src_st->codecpar->height;
    fifo->enc_pix_fmt = src_st->codecpar->format;
    fifo->enc_sample_aspect_ratio = src_st->codecpar->sample_aspect_ratio;
    fifo->enc_framerate = src_st->avg_frame_rate;
    fifo->enc_time_base = src_st->time_base;
    fifo->enc_color_primaries = src_st->codecpar->color_primaries;
    fifo->enc_color_trc = src_st->codecpar->color_trc;
    fifo->enc_colorspace = src_st->codecpar->color_space;
    fifo->enc_color_range = src_st->codecpar->color_range;
    fifo->enc_chroma_sample_location = src_st->codecpar->chroma_location;
    fifo->enc_field_order = src_st->codecpar->field_order;

    fifo->current_bitrate_idx = fifo->abr_levels - 1;

    ret = init_encoder(fifo, fifo->current_bitrate_idx, oformat);
    if (ret < 0) {
        return ret;
    }
    /* Allocate frames */
    fifo->decoded_frame = av_frame_alloc();
    fifo->enc_pkt = av_packet_alloc();

    if (!fifo->decoded_frame || !fifo->enc_pkt) {
        return AVERROR(ENOMEM);
    }

    fifo->last_switch_time = av_gettime_relative() / 1000; // Convert to ms

    return 0;
}


static int flush_write_packet(AVFormatContext* avf, AVPacket* pkt) {
    FifoContext* fifo = avf->priv_data;
    FifoMessage msg = { FIFO_WRITE_PACKET };
    int ret;

    ret = av_packet_ref(&msg.pkt, pkt);
    if (ret < 0) return ret;

    ret = av_thread_message_queue_send(fifo->queue, &msg,
        fifo->drop_pkts_on_overflow ?
        AV_THREAD_MESSAGE_NONBLOCK : 0);
    if (ret == AVERROR(EAGAIN)) {
        uint8_t overflow_set = 0;

        /* Queue is full, set fifo->overflow_flag to 1
         * to let consumer thread know the queue should
         * be flushed. */
        pthread_mutex_lock(&fifo->overflow_flag_lock);
        if (!fifo->overflow_flag)
            fifo->overflow_flag = overflow_set = 1;
        pthread_mutex_unlock(&fifo->overflow_flag_lock);

        if (overflow_set)
            av_log(avf, AV_LOG_WARNING, "FIFO queue full when flush encoder packets\n");
        ret = 0;
        goto fail;
    }
    else if (ret < 0) {
        goto fail;
    }

    return ret;
fail:
    av_packet_unref(&msg.pkt);
    return ret;
}


/* Switch to a different bitrate encoder */
static int switch_bitrate(AVFormatContext* avf, int new_idx)
{
    FifoContext* fifo = avf->priv_data;
    if (new_idx == fifo->current_bitrate_idx)
        return 0;

    if (new_idx < 0 || new_idx >= fifo->abr_levels)
        return AVERROR(EINVAL);

    av_log(fifo, AV_LOG_INFO, "Switching bitrate from level %d to level %d (%" PRId64 " -> %" PRId64 " bps)\n",
        fifo->current_bitrate_idx, new_idx,
        fifo->bitrate_levels[fifo->current_bitrate_idx].bitrate,
        fifo->bitrate_levels[new_idx].bitrate);

    /* Flush current encoder */
    AVCodecContext* old_ctx = fifo->enc_ctx;
    int ret = avcodec_send_frame(old_ctx, NULL);
    if (ret < 0 && ret != AVERROR_EOF) {
        av_log(fifo, AV_LOG_WARNING, "Error flushing encoder: %s\n", av_err2str(ret));
    }

    /* Receive remaining packets */
    while (ret >= 0) {
        ret = avcodec_receive_packet(old_ctx, fifo->enc_pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            break;
        if (ret < 0) {
            av_log(fifo, AV_LOG_WARNING, "Error receiving packet during flush: %s\n", av_err2str(ret));
            break;
        }
        flush_write_packet(avf, fifo->enc_pkt);
        av_packet_unref(fifo->enc_pkt);
    }
    ret = init_encoder(fifo, new_idx, fifo->oformat);
    if (ret < 0) {
        return ret;
    }

    fifo->current_bitrate_idx = new_idx;
    fifo->last_switch_time = av_gettime_relative() / 1000;

    return 0;
}

/* Determine target bitrate level based on current queue size and time */
static int determine_bitrate_level(FifoContext* fifo)
{
    int64_t now = av_gettime_relative() / 1000; // Current time in ms
    int64_t time_since_switch = now - fifo->last_switch_time;
    /* Check if we can switch (enough time passed since last switch) */

    /* Get current queue size and calculate percentage */
    int current_qsize = av_thread_message_queue_nb_elems(fifo->queue);
    //av_log(fifo, AV_LOG_INFO, "Current queue size %d \n", current_qsize);

    /* Calculate threshold values in queue size units */
    int low_threshold = (fifo->abr_low_threshold * fifo->queue_size) / 100;
    int high_threshold = (fifo->abr_high_threshold * fifo->queue_size) / 100;
    int new_idx = fifo->current_bitrate_idx;
    if (current_qsize <= low_threshold && time_since_switch >= fifo->abr_up_switch_interval) {
        if (fifo->current_bitrate_idx < fifo->abr_levels - 1) {
            /* Not at highest bitrate yet */
            /* Switch up one level */
            new_idx = fifo->current_bitrate_idx + 1;
        }
    }
    else if (current_qsize >= high_threshold && fifo->qsize_since_switch_trend > 0.1 && time_since_switch >= fifo->abr_down_switch_interval) {
        if (fifo->current_bitrate_idx > 0) {
            /* Not at lowest bitrate yet */
            /* Switch down one level */
            new_idx = fifo->current_bitrate_idx - 1;
        }
    }

    if (new_idx != fifo->current_bitrate_idx) {
        av_log(fifo, AV_LOG_INFO, "Current queue size: %d, queue size trend: %lf\n", current_qsize, fifo->qsize_since_switch_trend);
        fifo->qsize_since_switch_trend = 0;
    }

    fifo->qsize_since_switch_trend = (fifo->qsize_since_switch_trend * 9 + (current_qsize - fifo->last_qsize)) / 10;
    fifo->last_qsize = current_qsize;

    return new_idx;
}


static int packet_enqueue(AVFormatContext* avf, AVPacket* pkt) {
    int ret;
    FifoContext* fifo = avf->priv_data;
    FifoMessage msg = { .type = pkt ? FIFO_WRITE_PACKET : FIFO_FLUSH_OUTPUT };
    if (pkt) {
        ret = av_packet_ref(&msg.pkt, pkt);
        if (ret < 0) return ret;
    }
    ret = av_thread_message_queue_send(fifo->queue, &msg,
        fifo->drop_pkts_on_overflow ?
        AV_THREAD_MESSAGE_NONBLOCK : 0);
    if (ret == AVERROR(EAGAIN)) {
        uint8_t overflow_set = 0;

        pthread_mutex_lock(&fifo->overflow_flag_lock);
        if (!fifo->overflow_flag)
            fifo->overflow_flag = overflow_set = 1;
        pthread_mutex_unlock(&fifo->overflow_flag_lock);

        if (overflow_set)
            av_log(avf, AV_LOG_WARNING, "FIFO queue full\n");
        ret = 0;
        goto fail;
    }
    else if (ret < 0) {
        goto fail;
    }
    return ret;
fail:
    if (pkt)
        av_packet_unref(&msg.pkt);
    return ret;
}


static int encode_video_packet_and_enqueue(AVFormatContext* avf, AVPacket* pkt)
{
    int ret;
    FifoContext* fifo = avf->priv_data;
    AVCodecContext* enc_ctx = fifo->enc_ctx;
    AVPacket* enc_pkt = fifo->enc_pkt;
    /* Send packet to decoder */
    if (pkt) {
        ret = avcodec_send_packet(fifo->dec_ctx, pkt);
        if (ret < 0) {
            if (ret == AVERROR(EAGAIN)) {
                /* Drain decoded frames first */
                av_log(fifo, AV_LOG_DEBUG, "Decoder needs draining, will continue\n");
            }
            else {
                av_log(fifo, AV_LOG_ERROR, "Error sending packet to decoder: %s\n", av_err2str(ret));
                return ret;
            }
        }
    }
    else {
        /* Flush decoder */
        ret = avcodec_send_packet(fifo->dec_ctx, NULL);
        if (ret < 0) {
            av_log(fifo, AV_LOG_ERROR, "Error flushing decoder: %s\n", av_err2str(ret));
            return ret;
        }
    }

    /* Try to receive decoded frames and encode them */
    while (1) {
        /* Receive decoded frame */
        ret = avcodec_receive_frame(fifo->dec_ctx, fifo->decoded_frame);
        if (ret == AVERROR(EAGAIN)) {
            /* Need more packets for decoder */
            break;
        }
        if (ret == AVERROR_EOF) {
            /* Decoder has been flushed, now flush encoder */
            av_log(fifo, AV_LOG_DEBUG, "Decoder EOF reached, flushing encoder\n");

            /* Send flush to encoder */
            ret = avcodec_send_frame(enc_ctx, NULL);
            if (ret < 0) {
                av_log(fifo, AV_LOG_ERROR, "Error flushing encoder: %s\n", av_err2str(ret));
                return ret;
            }
            /* Receive all remaining packets from encoder */
            while (1) {
                ret = avcodec_receive_packet(enc_ctx, enc_pkt);
                if (ret == AVERROR(EAGAIN)) {
                    continue;  // Should not happen after flush, but just in case
                }
                if (ret == AVERROR_EOF) {
                    av_log(fifo, AV_LOG_DEBUG, "Encoder flushed successfully\n");
                    break;
                }
                if (ret < 0) {
                    av_log(fifo, AV_LOG_ERROR, "Error receiving packet from encoder after flush: %s\n", av_err2str(ret));
                    return ret;
                }
                enc_pkt->stream_index = fifo->video_stream_idx;
                packet_enqueue(avf, enc_pkt);
                av_packet_unref(enc_pkt);
            }
            break;  // Exit the while loop after handling EOF
        }
        if (ret < 0) {
            av_log(fifo, AV_LOG_ERROR, "Error receiving frame from decoder: %s\n", av_err2str(ret));
            return ret;
        }

        /* Prepare frame for encoder */
        fifo->decoded_frame->pts = fifo->decoded_frame->best_effort_timestamp;
        fifo->decoded_frame->pict_type = AV_PICTURE_TYPE_NONE; /* Let encoder decide */

        /* Send frame to encoder */
        ret = avcodec_send_frame(enc_ctx, fifo->decoded_frame);
        av_frame_unref(fifo->decoded_frame);

        if (ret < 0) {
            av_log(fifo, AV_LOG_ERROR, "Error sending frame to encoder: %s\n", av_err2str(ret));
            return ret;
        }

        /* Try to receive any pending packets from encoder */
        while (1) {
            ret = avcodec_receive_packet(enc_ctx, enc_pkt);
            if (ret == AVERROR(EAGAIN)) {
                /* Encoder needs more frames */
                break;
            }
            if (ret == AVERROR_EOF) {
                /* Encoder flushed - should not happen here as we haven't sent flush */
                av_log(fifo, AV_LOG_DEBUG, "Unexpected encoder EOF\n");
                break;
            }
            if (ret < 0) {
                av_log(fifo, AV_LOG_ERROR, "Error receiving packet from encoder: %s\n", av_err2str(ret));
                return ret;
            }

            /* Got an encoded packet */
            enc_pkt->stream_index = fifo->video_stream_idx;
            packet_enqueue(avf, enc_pkt);
            av_packet_unref(enc_pkt);
        }
    }

    return 0;
}

static int fifo_mux_init(AVFormatContext* avf, const AVOutputFormat* oformat,
    const char* filename)
{
    FifoContext* fifo = avf->priv_data;
    AVFormatContext* avf2;
    int ret = 0, i;
    if (fifo->avf) return 0;

    ret = avformat_alloc_output_context2(&avf2, oformat, NULL, filename);
    if (ret < 0)
        return ret;

    fifo->avf = avf2;

    avf2->interrupt_callback = avf->interrupt_callback;
    avf2->max_delay = avf->max_delay;
    ret = av_dict_copy(&avf2->metadata, avf->metadata, 0);
    if (ret < 0)
        return ret;
    avf2->opaque = avf->opaque;
    avf2->io_close2 = avf->io_close2;
    avf2->io_open = avf->io_open;
    avf2->flags = avf->flags;

    /* For ABR mode, don't clone video stream - we'll create it from encoder */
    if (fifo->enable_abr) {
        int video_stream_found = 0;

        for (i = 0; i < avf->nb_streams; i++) {
            AVStream* src_st = avf->streams[i];

            if (src_st->codecpar->codec_type == AVMEDIA_TYPE_VIDEO && !video_stream_found) {
                /* Create stream from encoder context (use highest bitrate encoder) */
                AVCodecContext* enc_ctx = fifo->enc_ctx;
                AVStream* st = avformat_new_stream(avf2, NULL);
                if (!st)
                    return AVERROR(ENOMEM);

                ret = avcodec_parameters_from_context(st->codecpar, enc_ctx);
                if (ret < 0)
                    return ret;

                st->time_base = enc_ctx->time_base;
                st->avg_frame_rate = enc_ctx->framerate;
                st->sample_aspect_ratio = src_st->sample_aspect_ratio;

                video_stream_found = 1;
                av_log(avf, AV_LOG_INFO, "Created video stream from encoder context (bitrate=%" PRId64 ")\n",
                    enc_ctx->bit_rate);
            }
            else {
                /* Clone non-video streams as-is */
                AVStream* st = ff_stream_clone(avf2, src_st);
                if (!st)
                    return AVERROR(ENOMEM);
            }
        }
    }
    else {
        /* Standard mode: clone all streams */
        for (i = 0; i < avf->nb_streams; ++i) {
            AVStream* st = ff_stream_clone(avf2, avf->streams[i]);
            if (!st)
                return AVERROR(ENOMEM);
        }
    }

    return 0;
}


static int fifo_thread_write_header(FifoThreadContext* ctx)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    AVFormatContext* avf2;
    AVDictionary* format_options = NULL;
    int ret, i;

    ret = av_dict_copy(&format_options, fifo->format_options, 0);
    if (ret < 0)
        goto end;

    ret = fifo_mux_init(avf, fifo->oformat, avf->url);
    if (ret < 0) {
        av_log(avf, AV_LOG_ERROR, "Error opening %s: %s\n", avf->url,
            av_err2str(ret));
        goto end;
    }
    avf2 = fifo->avf;

    ret = ff_format_output_open(avf2, avf->url, &format_options);
    if (ret < 0) {
        av_log(avf, AV_LOG_ERROR, "Error opening %s: %s\n", avf->url,
            av_err2str(ret));
        goto end;
    }

    ret = avformat_write_header(avf2, &format_options);
    if (!ret)
        ctx->header_written = 1;

    if (format_options) {
        const AVDictionaryEntry* entry = NULL;
        while ((entry = av_dict_iterate(format_options, entry)))
            av_log(avf2, AV_LOG_ERROR, "Unknown option '%s'\n", entry->key);
        ret = AVERROR(EINVAL);
    }

end:
    av_dict_free(&format_options);
    return ret;
}

static int fifo_thread_flush_output(FifoThreadContext* ctx)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    AVFormatContext* avf2 = fifo->avf;

    return av_interleaved_write_frame(avf2, NULL);
}

static int fifo_thread_write_packet(FifoThreadContext* ctx, AVPacket* pkt)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    AVFormatContext* avf2 = fifo->avf;
    AVRational src_tb, dst_tb;
    int ret, s_idx;
    int64_t orig_pts, orig_dts, orig_duration;

    if (ctx->drop_until_keyframe) {
        if (pkt->flags & AV_PKT_FLAG_KEY) {
            ctx->drop_until_keyframe = 0;
            av_log(avf, AV_LOG_VERBOSE, "Keyframe received, recovering...\n");
        }
        else {
            av_log(avf, AV_LOG_VERBOSE, "Dropping non-keyframe packet\n");
            av_packet_unref(pkt);
            return 0;
        }
    }
    ret = av_packet_ref(fifo->pkt_for_write, pkt);
    if (ret < 0) return ret;

    s_idx = pkt->stream_index;
    src_tb = avf->streams[s_idx]->time_base;
    dst_tb = avf2->streams[s_idx]->time_base;
    av_packet_rescale_ts(fifo->pkt_for_write, src_tb, dst_tb);
    ret = av_interleaved_write_frame(avf2, fifo->pkt_for_write);
    if (ret >= 0) {
        av_packet_unref(pkt);
    }
    return ret;
}

static int fifo_thread_write_trailer(FifoThreadContext* ctx)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    AVFormatContext* avf2 = fifo->avf;
    int ret;

    if (!ctx->header_written)
        return 0;

    ret = av_write_trailer(avf2);
    ff_format_io_close(avf2, &avf2->pb);
    avformat_free_context(fifo->avf);
    fifo->avf = NULL;
    return ret;
}

static int fifo_thread_dispatch_message(FifoThreadContext* ctx, FifoMessage* msg)
{
    int ret = AVERROR(EINVAL);

    if (msg->type == FIFO_NOOP)
        return 0;

    if (!ctx->header_written) {
        ret = fifo_thread_write_header(ctx);
        if (ret < 0)
            return ret;
    }

    switch (msg->type) {
    case FIFO_WRITE_HEADER:
        av_assert0(ret >= 0);
        return ret;
    case FIFO_WRITE_PACKET:
        return fifo_thread_write_packet(ctx, &msg->pkt);
    case FIFO_FLUSH_OUTPUT:
        return fifo_thread_flush_output(ctx);
    }

    av_assert0(0);
    return AVERROR(EINVAL);
}

static int is_recoverable(const FifoContext* fifo, int err_no) {
    if (!fifo->attempt_recovery)
        return 0;

    if (fifo->recover_any_error)
        return err_no != AVERROR_EXIT;

    switch (err_no) {
    case AVERROR(EINVAL):
    case AVERROR(ENOSYS):
    case AVERROR_EOF:
    case AVERROR_EXIT:
    case AVERROR_PATCHWELCOME:
        return 0;
    default:
        return 1;
    }
}

static void free_message(void* msg)
{
    FifoMessage* fifo_msg = msg;

    if (fifo_msg->type == FIFO_WRITE_PACKET)
        av_packet_unref(&fifo_msg->pkt);
}

static int fifo_thread_process_recovery_failure(FifoThreadContext* ctx, AVPacket* pkt,
    int err_no)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    int ret;

    av_log(avf, AV_LOG_INFO, "Recovery failed: %s\n",
        av_err2str(err_no));

    if (fifo->recovery_wait_streamtime) {
        if (pkt->pts == AV_NOPTS_VALUE)
            av_log(avf, AV_LOG_WARNING, "Packet does not contain presentation"
                " timestamp, recovery will be attempted immediately");
        ctx->last_recovery_ts = pkt->pts;
    }
    else {
        ctx->last_recovery_ts = av_gettime_relative();
    }

    if (fifo->max_recovery_attempts &&
        ctx->recovery_nr >= fifo->max_recovery_attempts) {
        av_log(avf, AV_LOG_ERROR,
            "Maximal number of %d recovery attempts reached.\n",
            fifo->max_recovery_attempts);
        ret = err_no;
    }
    else {
        ret = AVERROR(EAGAIN);
    }

    return ret;
}

static int fifo_thread_attempt_recovery(FifoThreadContext* ctx, FifoMessage* msg, int err_no)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    AVPacket* pkt = &msg->pkt;
    int64_t time_since_recovery;
    int ret;

    if (!is_recoverable(fifo, err_no)) {
        ret = err_no;
        goto fail;
    }

    if (ctx->header_written) {
        fifo->write_trailer_ret = fifo_thread_write_trailer(ctx);
        ctx->header_written = 0;
    }

    if (!ctx->recovery_nr) {
        ctx->last_recovery_ts = fifo->recovery_wait_streamtime ?
            AV_NOPTS_VALUE : 0;
    }
    else {
        if (fifo->recovery_wait_streamtime) {
            if (ctx->last_recovery_ts == AV_NOPTS_VALUE) {
                AVRational tb = avf->streams[pkt->stream_index]->time_base;
                time_since_recovery = av_rescale_q(pkt->pts - ctx->last_recovery_ts,
                    tb, AV_TIME_BASE_Q);
            }
            else {
                time_since_recovery = fifo->recovery_wait_time;
            }
        }
        else {
            time_since_recovery = av_gettime_relative() - ctx->last_recovery_ts;
        }

        if (time_since_recovery < fifo->recovery_wait_time)
            return AVERROR(EAGAIN);
    }

    ctx->recovery_nr++;

    if (fifo->max_recovery_attempts) {
        av_log(avf, AV_LOG_VERBOSE, "Recovery attempt #%d/%d\n",
            ctx->recovery_nr, fifo->max_recovery_attempts);
    }
    else {
        av_log(avf, AV_LOG_VERBOSE, "Recovery attempt #%d\n",
            ctx->recovery_nr);
    }

    if (fifo->restart_with_keyframe && fifo->drop_pkts_on_overflow)
        ctx->drop_until_keyframe = 1;

    ret = fifo_thread_dispatch_message(ctx, msg);
    if (ret < 0) {
        if (is_recoverable(fifo, ret)) {
            return fifo_thread_process_recovery_failure(ctx, pkt, ret);
        }
        else {
            goto fail;
        }
    }
    else {
        av_log(avf, AV_LOG_INFO, "Recovery successful\n");
        ctx->recovery_nr = 0;
    }

    return 0;

fail:
    free_message(msg);
    return ret;
}

static int fifo_thread_recover(FifoThreadContext* ctx, FifoMessage* msg, int err_no)
{
    AVFormatContext* avf = ctx->avf;
    FifoContext* fifo = avf->priv_data;
    int ret;

    do {
        if (!fifo->recovery_wait_streamtime && ctx->recovery_nr > 0) {
            int64_t time_since_recovery = av_gettime_relative() - ctx->last_recovery_ts;
            int64_t time_to_wait = FFMAX(0, fifo->recovery_wait_time - time_since_recovery);
            if (time_to_wait)
                av_usleep(FFMIN(10000, time_to_wait));
        }

        ret = fifo_thread_attempt_recovery(ctx, msg, err_no);
    } while (ret == AVERROR(EAGAIN) && !fifo->drop_pkts_on_overflow);

    if (ret == AVERROR(EAGAIN) && fifo->drop_pkts_on_overflow) {
        if (msg->type == FIFO_WRITE_PACKET)
            av_packet_unref(&msg->pkt);
        ret = 0;
    }

    return ret;
}

static void* fifo_consumer_thread(void* data)
{
    AVFormatContext* avf = data;
    FifoContext* fifo = avf->priv_data;
    AVThreadMessageQueue* queue = fifo->queue;
    FifoMessage msg = { FIFO_WRITE_HEADER, {0} };
    int ret;

    FifoThreadContext fifo_thread_ctx;
    memset(&fifo_thread_ctx, 0, sizeof(FifoThreadContext));
    fifo_thread_ctx.avf = avf;
    fifo_thread_ctx.last_received_dts = AV_NOPTS_VALUE;

    ff_thread_setname("fifo-consumer");

    while (1) {
        uint8_t just_flushed = 0;

        if (!fifo_thread_ctx.recovery_nr)
            ret = fifo_thread_dispatch_message(&fifo_thread_ctx, &msg);

        if (ret < 0 || fifo_thread_ctx.recovery_nr > 0) {
            int rec_ret = fifo_thread_recover(&fifo_thread_ctx, &msg, ret);
            if (rec_ret < 0) {
                av_thread_message_queue_set_err_send(queue, rec_ret);
                break;
            }
        }

        pthread_mutex_lock(&fifo->overflow_flag_lock);
        if (fifo->overflow_flag) {
            av_thread_message_flush(queue);
            if (fifo->restart_with_keyframe)
                fifo_thread_ctx.drop_until_keyframe = 1;
            fifo->overflow_flag = 0;
            just_flushed = 1;
        }
        pthread_mutex_unlock(&fifo->overflow_flag_lock);

        if (just_flushed)
            av_log(avf, AV_LOG_INFO, "FIFO queue flushed\n");

        ret = av_thread_message_queue_recv(queue, &msg, 0);
        if (ret < 0) {
            av_thread_message_queue_set_err_send(queue, ret);
            break;
        }
    }

    fifo->write_trailer_ret = fifo_thread_write_trailer(&fifo_thread_ctx);

    return NULL;
}


static int fifo_init(AVFormatContext* avf)
{
    FifoContext* fifo = avf->priv_data;
    const AVOutputFormat* oformat;
    int ret = 0;

    if (fifo->recovery_wait_streamtime && !fifo->drop_pkts_on_overflow) {
        av_log(avf, AV_LOG_ERROR, "recovery_wait_streamtime can be turned on"
            " only when drop_pkts_on_overflow is also turned on\n");
        return AVERROR(EINVAL);
    }

    oformat = av_guess_format(fifo->format, avf->url, NULL);
    if (!oformat) {
        ret = AVERROR_MUXER_NOT_FOUND;
        return ret;
    }
    fifo->oformat = oformat;

    /* Initialize ABR encoders if enabled */
    if (fifo->enable_abr) {
        if (!fifo->abr_codec) {
            av_log(avf, AV_LOG_ERROR, "ABR enabled but no codec specified\n");
            return AVERROR(EINVAL);
        }

        /* Parse bitrate string */
        ret = parse_abr_bitrates(fifo);
        if (ret < 0)
            return ret;

        ret = init_decoder_and_encoder(fifo, avf, oformat);
        if (ret < 0) {
            av_log(avf, AV_LOG_ERROR, "Failed to initialize ABR encoders: %s\n", av_err2str(ret));
            return ret;
        }
    }

    ret = fifo_mux_init(avf, oformat, avf->url);
    if (ret < 0)
        return ret;

    ret = av_thread_message_queue_alloc(&fifo->queue, (unsigned)fifo->queue_size,
        sizeof(FifoMessage));
    if (ret < 0)
        return ret;

    av_thread_message_queue_set_free_func(fifo->queue, free_message);

    fifo->pkt_for_write = av_packet_alloc();
    if (!fifo->pkt_for_write)
        return AVERROR(ENOMEM);

    ret = pthread_mutex_init(&fifo->overflow_flag_lock, NULL);
    if (ret < 0)
        return AVERROR(ret);
    fifo->overflow_flag_lock_initialized = 1;

    return 0;
}

static int fifo_write_header(AVFormatContext* avf)
{
    FifoContext* fifo = avf->priv_data;
    int ret;

    ret = pthread_create(&fifo->writer_thread, NULL, fifo_consumer_thread, avf);
    if (ret) {
        av_log(avf, AV_LOG_ERROR, "Failed to start thread: %s\n",
            av_err2str(AVERROR(ret)));
        ret = AVERROR(ret);
    }

    return ret;
}

static int fifo_write_packet(AVFormatContext* avf, AVPacket* pkt)
{
    FifoContext* fifo = avf->priv_data;
    int ret;

    if (fifo->enable_abr) {
        /* Check if we need to switch bitrate based on queue state */
        int target_level = determine_bitrate_level(fifo);
        if (target_level != fifo->current_bitrate_idx) {
            ret = switch_bitrate(avf, target_level);
            if (ret < 0) {
                return ret;
            }
        }
    }
    /* ABR encoding logic for video packets */
    if (fifo->enable_abr && pkt && pkt->stream_index == fifo->video_stream_idx) {
        ret = encode_video_packet_and_enqueue(avf, pkt);
    }
    else {
        /* Non-ABR mode or non-video stream: pass through */
        ret = packet_enqueue(avf, pkt);
    }
    return ret;
}

static int fifo_write_trailer(AVFormatContext* avf)
{
    FifoContext* fifo = avf->priv_data;
    int ret;

    if (fifo->enable_abr) {
        encode_video_packet_and_enqueue(avf, NULL);
    }
    packet_enqueue(avf, NULL);

    av_thread_message_queue_set_err_recv(fifo->queue, AVERROR_EOF);

    ret = pthread_join(fifo->writer_thread, NULL);
    if (ret < 0) {
        av_log(avf, AV_LOG_ERROR, "pthread join error: %s\n",
            av_err2str(AVERROR(ret)));
        return AVERROR(ret);
    }

    ret = fifo->write_trailer_ret;
    return ret;
}

static void fifo_deinit(AVFormatContext* avf)
{
    FifoContext* fifo = avf->priv_data;

    /* Cleanup ABR resources */
    if (fifo->enable_abr) {
        /* Free encoders */
        if (fifo->enc_ctx) {
            avcodec_free_context(&fifo->enc_ctx);
        }

        if (fifo->dec_ctx)
            avcodec_free_context(&fifo->dec_ctx);

        av_frame_free(&fifo->decoded_frame);
        av_packet_free(&fifo->enc_pkt);
    }

    if (fifo->pkt_for_write)
        av_packet_free(&fifo->pkt_for_write);

    avformat_free_context(fifo->avf);
    av_thread_message_queue_free(&fifo->queue);
    if (fifo->overflow_flag_lock_initialized)
        pthread_mutex_destroy(&fifo->overflow_flag_lock);
}

#define OFFSET(x) offsetof(FifoContext, x)
static const AVOption options[] = {
        {"attempt_recovery", "Attempt recovery in case of failure", OFFSET(attempt_recovery),
        AV_OPT_TYPE_BOOL, {.i64 = 0}, 0, 1, AV_OPT_FLAG_ENCODING_PARAM},

        {"drop_pkts_on_overflow", "Drop packets on fifo queue overflow not to block encoder", OFFSET(drop_pkts_on_overflow),
         AV_OPT_TYPE_BOOL, {.i64 = 0}, 0, 1, AV_OPT_FLAG_ENCODING_PARAM},

        {"fifo_format", "Target muxer", OFFSET(format),
         AV_OPT_TYPE_STRING, {.str = NULL}, 0, 0, AV_OPT_FLAG_ENCODING_PARAM},

        {"format_opts", "Options to be passed to underlying muxer", OFFSET(format_options),
         AV_OPT_TYPE_DICT, {.str = NULL}, 0, 0, AV_OPT_FLAG_ENCODING_PARAM},

        {"max_recovery_attempts", "Maximal number of recovery attempts", OFFSET(max_recovery_attempts),
         AV_OPT_TYPE_INT, {.i64 = FIFO_DEFAULT_MAX_RECOVERY_ATTEMPTS}, 0, INT_MAX, AV_OPT_FLAG_ENCODING_PARAM},

        {"queue_size", "Size of fifo queue", OFFSET(queue_size),
         AV_OPT_TYPE_INT, {.i64 = FIFO_DEFAULT_QUEUE_SIZE}, 1, INT_MAX, AV_OPT_FLAG_ENCODING_PARAM},

        {"recovery_wait_streamtime", "Use stream time instead of real time while waiting for recovery",
         OFFSET(recovery_wait_streamtime), AV_OPT_TYPE_BOOL, {.i64 = 0}, 0, 1, AV_OPT_FLAG_ENCODING_PARAM},

        {"recovery_wait_time", "Waiting time between recovery attempts", OFFSET(recovery_wait_time),
         AV_OPT_TYPE_DURATION, {.i64 = FIFO_DEFAULT_RECOVERY_WAIT_TIME_USEC}, 0, INT64_MAX, AV_OPT_FLAG_ENCODING_PARAM},

        {"recover_any_error", "Attempt recovery regardless of type of the error", OFFSET(recover_any_error),
         AV_OPT_TYPE_BOOL, {.i64 = 0}, 0, 1, AV_OPT_FLAG_ENCODING_PARAM},

        {"restart_with_keyframe", "Wait for keyframe when restarting output", OFFSET(restart_with_keyframe),
         AV_OPT_TYPE_BOOL, {.i64 = 0}, 0, 1, AV_OPT_FLAG_ENCODING_PARAM},



         /* Adaptive Bitrate Options */
         {"enable_abr", "Enable adaptive bitrate encoding", OFFSET(enable_abr),
          AV_OPT_TYPE_BOOL, {.i64 = 0}, 0, 1, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_codec", "Encoder for adaptive bitrate (e.g., libx264)", OFFSET(abr_codec),
          AV_OPT_TYPE_STRING, {.str = NULL}, 0, 0, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_codec_opts", "Options for ABR encoder", OFFSET(abr_codec_options),
          AV_OPT_TYPE_DICT, {.str = NULL}, 0, 0, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_bitrates", "Comma-separated bitrates low-to-high (e.g., '500k,1M,2M,4M')", OFFSET(abr_bitrates),
          AV_OPT_TYPE_STRING, {.str = NULL}, 0, 0, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_low_threshold", "Queue percentage to switch to higher bitrate (default 5)", OFFSET(abr_low_threshold),
          AV_OPT_TYPE_INT, {.i64 = 5}, 0, 100, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_high_threshold", "Queue percentage to switch to lower bitrate (default 30)", OFFSET(abr_high_threshold),
          AV_OPT_TYPE_INT, {.i64 = 30}, 0, 100, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_up_switch_interval", "Milliseconds between bitrate up switches (default 30000)", OFFSET(abr_up_switch_interval),
          AV_OPT_TYPE_INT64, {.i64 = 30000}, 0, 60000, AV_OPT_FLAG_ENCODING_PARAM},

         {"abr_down_switch_interval", "Milliseconds between bitrate down switches (default 5000)", OFFSET(abr_down_switch_interval),
          AV_OPT_TYPE_INT64, {.i64 = 5000}, 0, 60000, AV_OPT_FLAG_ENCODING_PARAM},

         {NULL},
};

static const AVClass fifo_adaptive_muxer_class = {
    .class_name = "adaptive fifo muxer",
    .item_name = av_default_item_name,
    .option = options,
    .version = LIBAVUTIL_VERSION_INT,
};

const FFOutputFormat ff_fifo_adaptive_muxer = {
    .p.name = "fifo_adaptive",
    .p.long_name = NULL_IF_CONFIG_SMALL("Adaptive bitrate FIFO queue muxer"),
    .p.priv_class = &fifo_adaptive_muxer_class,
    .p.flags = AVFMT_NOFILE | AVFMT_TS_NEGATIVE,
    .priv_data_size = sizeof(FifoContext),
    .init = fifo_init,
    .write_header = fifo_write_header,
    .write_packet = fifo_write_packet,
    .write_trailer = fifo_write_trailer,
    .deinit = fifo_deinit,
    .flags_internal = FF_OFMT_FLAG_ALLOW_FLUSH,
};