/*
 * HEVC Supplementary Enhancement Information messages
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef AVCODEC_HEVC_SEI_H
#define AVCODEC_HEVC_SEI_H

#include <stdint.h>

#include "libavutil/buffer.h"

#include "libavcodec/get_bits.h"
#include "libavcodec/h2645_sei.h"
#include "libavcodec/sei.h"

#include "hevc.h"


typedef enum {
    // SEI Picture Timing Picture Structure.
    // From the ITU-T H.265 Standards Document v3 (04/2015):
    // Table D.2: Interpretation of pic_struct:
    // When present, pic_struct is constrained to use one of the following:
    //   - all pictures in CSV are one of: 0, 7 or 8.
    //   - all pictures in CSV are one of: 1, 2, 9, 10, 11, or 12..
    //   - all pictures in CSV are one of: 3, 4, 5 or 6...

    // progressive frame.
    HEVC_SEI_PIC_STRUCT_FRAME_PROGRESSIVE = 0,

    // top field.
    HEVC_SEI_PIC_STRUCT_FIELD_TOP         = 1,
    // bottom field.
    HEVC_SEI_PIC_STRUCT_FIELD_BOTTOM      = 2,

    // top field, bottom field, in that order. Top Field First.
    HEVC_SEI_PIC_STRUCT_FRAME_TFBF        = 3,
    // bottom Field, top field, in that order. Bottom Field First.
    HEVC_SEI_PIC_STRUCT_FRAME_BFTF        = 4,

    // top field, bottom field, top field repeated, Top Field First.
    HEVC_SEI_PIC_STRUCT_FRAME_TFBFTF      = 5,
    // bottom field, top field, bottom field repeated, Bottom Field First.
    HEVC_SEI_PIC_STRUCT_FRAME_BFTFBF      = 6,

    // frame doubling.
    HEVC_SEI_PIC_STRUCT_FRAME_DOUBLING    = 7,
    // frame trippling.
    HEVC_SEI_PIC_STRUCT_FRAME_TRIPLING    = 8,

    // top field paired with previous bottom field. Bottom Field First.
    HEVC_SEI_PIC_STRUCT_FIELD_TFPBF       = 9,
    // bottom field paired with previous top field. Top Field First.
    HEVC_SEI_PIC_STRUCT_FIELD_BFPTF       = 10,

    // top field paired with next bottom field. Top Field First.
    HEVC_SEI_PIC_STRUCT_FIELD_TFNBF       = 11,
    // bottom field paired with next top field. Bottom Field First.
    HEVC_SEI_PIC_STRUCT_FIELD_BFNTF       = 12,
} HEVC_SEI_PicStructType;

// Returns 1 - when type is interlaced, 0 - otherwise.
static inline int ff_hevc_sei_pic_struct_is_interlaced(HEVC_SEI_PicStructType type)
{
    switch (type) {
    case HEVC_SEI_PIC_STRUCT_FIELD_TOP:
    case HEVC_SEI_PIC_STRUCT_FIELD_BOTTOM:
    case HEVC_SEI_PIC_STRUCT_FRAME_TFBF:
    case HEVC_SEI_PIC_STRUCT_FRAME_BFTF:
    case HEVC_SEI_PIC_STRUCT_FRAME_TFBFTF:
    case HEVC_SEI_PIC_STRUCT_FRAME_BFTFBF:
    case HEVC_SEI_PIC_STRUCT_FIELD_TFPBF:
    case HEVC_SEI_PIC_STRUCT_FIELD_BFPTF:
    case HEVC_SEI_PIC_STRUCT_FIELD_TFNBF:
    case HEVC_SEI_PIC_STRUCT_FIELD_BFNTF:
        return 1;
    default:
        return 0;
    }
}

// Returns 1 - when type is top field first, 0 - otherwise.
static inline int ff_hevc_sei_pic_struct_is_tff(HEVC_SEI_PicStructType type)
{
    switch (type) {
    case HEVC_SEI_PIC_STRUCT_FRAME_TFBF:
    case HEVC_SEI_PIC_STRUCT_FRAME_TFBFTF:
    case HEVC_SEI_PIC_STRUCT_FIELD_BFPTF:
    case HEVC_SEI_PIC_STRUCT_FIELD_TFNBF:
        return 1;
    default:
        return 0;
    }
}

// Returns 1 - when type is bottom field first, 0 - otherwise.
static inline int ff_hevc_sei_pic_struct_is_bff(HEVC_SEI_PicStructType type)
{
    switch (type) {
    case HEVC_SEI_PIC_STRUCT_FRAME_BFTF:
    case HEVC_SEI_PIC_STRUCT_FRAME_BFTFBF:
    case HEVC_SEI_PIC_STRUCT_FIELD_TFPBF:
    case HEVC_SEI_PIC_STRUCT_FIELD_BFNTF:
        return 1;
    default:
        return 0;
    }
}

// Returns 1 - when type is top field, 0 - otherwise.
static inline int ff_hevc_sei_pic_struct_is_tf(HEVC_SEI_PicStructType type)
{
    switch (type) {
    case HEVC_SEI_PIC_STRUCT_FIELD_TOP:
    case HEVC_SEI_PIC_STRUCT_FIELD_TFPBF:
    case HEVC_SEI_PIC_STRUCT_FIELD_TFNBF:
        return 1;
    default:
        return 0;
    }
}

// Returns 1 - when type is bottom field, 0 - otherwise.
static inline int ff_hevc_sei_pic_struct_is_bf(HEVC_SEI_PicStructType type)
{
    switch (type) {
    case HEVC_SEI_PIC_STRUCT_FIELD_BOTTOM:
    case HEVC_SEI_PIC_STRUCT_FIELD_BFPTF:
    case HEVC_SEI_PIC_STRUCT_FIELD_BFNTF:
        return 1;
    default:
        return 0;
    }
}

// Returns 1 - when type is a field picture, 0 - otherwise.
static inline int ff_hevc_sei_pict_struct_is_field_picture(HEVC_SEI_PicStructType type)
{
    return (ff_hevc_sei_pic_struct_is_tf(type) || ff_hevc_sei_pic_struct_is_bf(type)) ? 1 : 0;
}

// Returns 1 - when type is a frame picture, 0 - otherwise.
static inline int ff_hevc_sei_pict_struct_is_frame_picture(HEVC_SEI_PicStructType type)
{
    return ff_hevc_sei_pict_struct_is_field_picture(type) ? 0 : 1;
}

typedef struct HEVCSEIPictureHash {
    uint8_t       md5[3][16];
    uint8_t is_md5;
} HEVCSEIPictureHash;

typedef struct HEVCSEIFramePacking {
    int present;
    int arrangement_type;
    int content_interpretation_type;
    int quincunx_subsampling;
    int current_frame_is_frame0_flag;
} HEVCSEIFramePacking;

typedef struct HEVCSEIPictureTiming {
    int picture_struct;
} HEVCSEIPictureTiming;

typedef struct HEVCSEIAlternativeTransfer {
    int present;
    int preferred_transfer_characteristics;
} HEVCSEIAlternativeTransfer;

typedef struct HEVCSEITimeCode {
    int      present;
    uint8_t  num_clock_ts;
    uint8_t  clock_timestamp_flag[3];
    uint8_t  units_field_based_flag[3];
    uint8_t  counting_type[3];
    uint8_t  full_timestamp_flag[3];
    uint8_t  discontinuity_flag[3];
    uint8_t  cnt_dropped_flag[3];
    uint16_t n_frames[3];
    uint8_t  seconds_value[3];
    uint8_t  minutes_value[3];
    uint8_t  hours_value[3];
    uint8_t  seconds_flag[3];
    uint8_t  minutes_flag[3];
    uint8_t  hours_flag[3];
    uint8_t  time_offset_length[3];
    int32_t  time_offset_value[3];
} HEVCSEITimeCode;

typedef struct HEVCSEITDRDI {
    uint8_t prec_ref_display_width;
    uint8_t ref_viewing_distance_flag;
    uint8_t prec_ref_viewing_dist;
    uint8_t num_ref_displays;
    uint16_t left_view_id[32];
    uint16_t right_view_id[32];
    uint8_t exponent_ref_display_width[32];
    uint64_t mantissa_ref_display_width[32];
    uint8_t exponent_ref_viewing_distance[32];
    uint64_t mantissa_ref_viewing_distance[32];
    uint8_t additional_shift_present_flag[32];
    int16_t num_sample_shift[32];
    uint8_t three_dimensional_reference_displays_extension_flag;
    int present;
} HEVCSEITDRDI;

typedef struct HEVCSEIRecoveryPoint {
    int16_t recovery_poc_cnt;
    uint8_t exact_match_flag;
    uint8_t broken_link_flag;
    uint8_t has_recovery_poc;
} HEVCSEIRecoveryPoint;

typedef struct HEVCSEI {
    H2645SEI common;
    HEVCSEIPictureHash picture_hash;
    HEVCSEIPictureTiming picture_timing;
    int active_seq_parameter_set_id;
    HEVCSEITimeCode timecode;
    HEVCSEITDRDI tdrdi;
    HEVCSEIRecoveryPoint recovery_point;
} HEVCSEI;

struct HEVCParamSets;

int ff_hevc_decode_nal_sei(GetBitContext *gb, void *logctx, HEVCSEI *s,
                           const struct HEVCParamSets *ps, enum HEVCNALUnitType type);

/**
 * Reset SEI values that are stored on the Context.
 * e.g. Caption data that was extracted during NAL
 * parsing.
 *
 * @param sei HEVCSEI.
 */
static inline void ff_hevc_reset_sei(HEVCSEI *sei)
{
    sei->timecode.present = 0;
    sei->tdrdi.present = 0;
    ff_h2645_sei_reset(&sei->common);
}

#endif /* AVCODEC_HEVC_SEI_H */
