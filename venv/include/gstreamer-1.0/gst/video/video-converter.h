/* Video conversion api function
 * Copyright (C) 2014 Wim Taymans <wim.taymans@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __GST_VIDEO_CONVERTER_H__
#define __GST_VIDEO_CONVERTER_H__

#include <gst/video/video.h>

G_BEGIN_DECLS

/**
 * GST_VIDEO_CONVERTER_OPT_RESAMPLER_METHOD:
 *
 * #GstVideoResamplerMethod, The resampler method to use for
 * resampling. Other options for the resampler can be used, see
 * the #GstVideoResampler. Default is #GST_VIDEO_RESAMPLER_METHOD_CUBIC
 */
#define GST_VIDEO_CONVERTER_OPT_RESAMPLER_METHOD   "GstVideoConverter.resampler-method"
/**
 * GST_VIDEO_CONVERTER_OPT_CHROMA_RESAMPLER_METHOD:
 *
 * #GstVideoChromaMethod, The resampler method to use for
 * chroma resampling. Other options for the resampler can be used, see
 * the #GstVideoResampler. Default is #GST_VIDEO_RESAMPLER_METHOD_LINEAR
 */
#define GST_VIDEO_CONVERTER_OPT_CHROMA_RESAMPLER_METHOD   "GstVideoConverter.chroma-resampler-method"
/**
 * GST_VIDEO_CONVERTER_OPT_RESAMPLER_TAPS:
 *
 * #G_TYPE_UINT, The number of taps for the resampler.
 * Default is 0: let the resampler choose a good value.
 */
#define GST_VIDEO_CONVERTER_OPT_RESAMPLER_TAPS   "GstVideoConverter.resampler-taps"

/**
 * GST_VIDEO_CONVERTER_OPT_DITHER_METHOD:
 *
 * #GstVideoDitherMethod, The dither method to use when
 * changing bit depth.
 * Default is #GST_VIDEO_DITHER_BAYER.
 */
#define GST_VIDEO_CONVERTER_OPT_DITHER_METHOD   "GstVideoConverter.dither-method"

/**
 * GST_VIDEO_CONVERTER_OPT_DITHER_QUANTIZATION:
 *
 * #G_TYPE_UINT, The quantization amount to dither to. Components will be
 * quantized to multiples of this value.
 * Default is 1
 */
#define GST_VIDEO_CONVERTER_OPT_DITHER_QUANTIZATION   "GstVideoConverter.dither-quantization"

/**
 * GST_VIDEO_CONVERTER_OPT_SRC_X:
 *
 * #G_TYPE_INT, source x position to start conversion, default 0
 */
#define GST_VIDEO_CONVERTER_OPT_SRC_X   "GstVideoConverter.src-x"
/**
 * GST_VIDEO_CONVERTER_OPT_SRC_Y:
 *
 * #G_TYPE_INT, source y position to start conversion, default 0
 */
#define GST_VIDEO_CONVERTER_OPT_SRC_Y   "GstVideoConverter.src-y"
/**
 * GST_VIDEO_CONVERTER_OPT_SRC_WIDTH:
 *
 * #G_TYPE_INT, source width to convert, default source width
 */
#define GST_VIDEO_CONVERTER_OPT_SRC_WIDTH   "GstVideoConverter.src-width"
/**
 * GST_VIDEO_CONVERTER_OPT_SRC_HEIGHT:
 *
 * #G_TYPE_INT, source height to convert, default source height
 */
#define GST_VIDEO_CONVERTER_OPT_SRC_HEIGHT   "GstVideoConverter.src-height"

/**
 * GST_VIDEO_CONVERTER_OPT_DEST_X:
 *
 * #G_TYPE_INT, x position in the destination frame, default 0
 */
#define GST_VIDEO_CONVERTER_OPT_DEST_X   "GstVideoConverter.dest-x"
/**
 * GST_VIDEO_CONVERTER_OPT_DEST_Y:
 *
 * #G_TYPE_INT, y position in the destination frame, default 0
 */
#define GST_VIDEO_CONVERTER_OPT_DEST_Y   "GstVideoConverter.dest-y"
/**
 * GST_VIDEO_CONVERTER_OPT_DEST_WIDTH:
 *
 * #G_TYPE_INT, width in the destination frame, default destination width
 */
#define GST_VIDEO_CONVERTER_OPT_DEST_WIDTH   "GstVideoConverter.dest-width"
/**
 * GST_VIDEO_CONVERTER_OPT_DEST_HEIGHT:
 *
 * #G_TYPE_INT, height in the destination frame, default destination height
 */
#define GST_VIDEO_CONVERTER_OPT_DEST_HEIGHT   "GstVideoConverter.dest-height"

/**
 * GST_VIDEO_CONVERTER_OPT_FILL_BORDER:
 *
 * #G_TYPE_BOOLEAN, if the destination rectangle does not fill the complete
 * destination image, render a border with
 * #GST_VIDEO_CONVERTER_OPT_BORDER_ARGB. Otherwise the unusded pixels in the
 * destination are untouched. Default %TRUE.
 */
#define GST_VIDEO_CONVERTER_OPT_FILL_BORDER   "GstVideoConverter.fill-border"
/**
 * GST_VIDEO_CONVERTER_OPT_ALPHA_VALUE:
 *
 * #G_TYPE_DOUBLE, the alpha color value to use.
 * Default to 1.0
 */
#define GST_VIDEO_CONVERTER_OPT_ALPHA_VALUE   "GstVideoConverter.alpha-value"
/**
 * GstVideoAlphaMode:
 * @GST_VIDEO_ALPHA_MODE_COPY: When input and output have alpha, it will be copied.
 *         When the input has no alpha, alpha will be set to
 *         #GST_VIDEO_CONVERTER_OPT_ALPHA_VALUE
 * @GST_VIDEO_ALPHA_MODE_SET: set all alpha to
 *	   #GST_VIDEO_CONVERTER_OPT_ALPHA_VALUE
 * @GST_VIDEO_ALPHA_MODE_MULT:  multiply all alpha with
 *         #GST_VIDEO_CONVERTER_OPT_ALPHA_VALUE.
 *         When the input format has no alpha but the output format has, the
 *         alpha value will be set to #GST_VIDEO_CONVERTER_OPT_ALPHA_VALUE
 *
 * Different alpha modes.
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_ALPHA_MODE_COPY,
  GST_VIDEO_ALPHA_MODE_SET,
  GST_VIDEO_ALPHA_MODE_MULT
} GstVideoAlphaMode;
/**
 * GST_VIDEO_CONVERTER_OPT_ALPHA_MODE:
 *
 * #GstVideoAlphaMode, the alpha mode to use.
 * Default is #GST_VIDEO_ALPHA_MODE_COPY.
 */
#define GST_VIDEO_CONVERTER_OPT_ALPHA_MODE   "GstVideoConverter.alpha-mode"
/**
 * GST_VIDEO_CONVERTER_OPT_BORDER_ARGB:
 *
 * #G_TYPE_UINT, the border color to use if #GST_VIDEO_CONVERTER_OPT_FILL_BORDER
 * is set to %TRUE. The color is in ARGB format.
 * Default 0xff000000
 */
#define GST_VIDEO_CONVERTER_OPT_BORDER_ARGB   "GstVideoConverter.border-argb"

/**
 * GstVideoChromaMode:
 * @GST_VIDEO_CHROMA_MODE_FULL: do full chroma up and down sampling
 * @GST_VIDEO_CHROMA_MODE_UPSAMPLE_ONLY: only perform chroma upsampling
 * @GST_VIDEO_CHROMA_MODE_DOWNSAMPLE_ONLY: only perform chroma downsampling
 * @GST_VIDEO_CHROMA_MODE_NONE: disable chroma resampling
 *
 * Different chroma downsampling and upsampling modes
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_CHROMA_MODE_FULL,
  GST_VIDEO_CHROMA_MODE_UPSAMPLE_ONLY,
  GST_VIDEO_CHROMA_MODE_DOWNSAMPLE_ONLY,
  GST_VIDEO_CHROMA_MODE_NONE
} GstVideoChromaMode;

/**
 * GST_VIDEO_CONVERTER_OPT_CHROMA_MODE:
 *
 * #GstVideoChromaMode, set the chroma resample mode subsampled
 * formats. Default is #GST_VIDEO_CHROMA_MODE_FULL.
 */
#define GST_VIDEO_CONVERTER_OPT_CHROMA_MODE   "GstVideoConverter.chroma-mode"

/**
 *GstVideoMatrixMode:
 * @GST_VIDEO_MATRIX_MODE_FULL: do conversion between color matrices
 * @GST_VIDEO_MATRIX_MODE_INPUT_ONLY:  use the input color matrix to convert
 *	  to and from R'G'B
 * @GST_VIDEO_MATRIX_MODE_OUTPUT_ONLY: use the output color matrix to convert
 *	  to and from R'G'B
 * @GST_VIDEO_MATRIX_MODE_NONE: disable color matrix conversion.
 *
 * Different color matrix conversion modes
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_MATRIX_MODE_FULL,
  GST_VIDEO_MATRIX_MODE_INPUT_ONLY,
  GST_VIDEO_MATRIX_MODE_OUTPUT_ONLY,
  GST_VIDEO_MATRIX_MODE_NONE
} GstVideoMatrixMode;
/**
 * GST_VIDEO_CONVERTER_OPT_MATRIX_MODE:
 *
 * #GstVideoMatrixMode, set the color matrix conversion mode for
 *	converting between Y'PbPr and non-linear RGB (R'G'B').
 * Default is #GST_VIDEO_MATRIX_MODE_FULL.
 */
#define GST_VIDEO_CONVERTER_OPT_MATRIX_MODE   "GstVideoConverter.matrix-mode"
/**
 * GstVideoGammaMode:
 * @GST_VIDEO_GAMMA_MODE_NONE: disable gamma handling
 * @GST_VIDEO_GAMMA_MODE_REMAP: convert between input and output gamma
 * Different gamma conversion modes
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_GAMMA_MODE_NONE,
  GST_VIDEO_GAMMA_MODE_REMAP
} GstVideoGammaMode;
/**
 * GST_VIDEO_CONVERTER_OPT_GAMMA_MODE:
 *
 * #GstVideoGammaMode, set the gamma mode.
 * Default is #GST_VIDEO_GAMMA_MODE_NONE.
 */
#define GST_VIDEO_CONVERTER_OPT_GAMMA_MODE   "GstVideoConverter.gamma-mode"
/**
 * GstVideoPrimariesMode:
 * @GST_VIDEO_PRIMARIES_MODE_NONE: disable conversion between primaries
 * @GST_VIDEO_PRIMARIES_MODE_MERGE_ONLY: do conversion between primaries only
 *	  when it can be merged with color matrix conversion.
 * @GST_VIDEO_PRIMARIES_MODE_FAST: fast conversion between primaries
 *
 * Different primaries conversion modes
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_PRIMARIES_MODE_NONE,
  GST_VIDEO_PRIMARIES_MODE_MERGE_ONLY,
  GST_VIDEO_PRIMARIES_MODE_FAST
} GstVideoPrimariesMode;
/**
 * GST_VIDEO_CONVERTER_OPT_PRIMARIES_MODE:
 *
 * #GstVideoPrimariesMode, set the primaries conversion mode.
 * Default is #GST_VIDEO_PRIMARIES_MODE_NONE.
 */
#define GST_VIDEO_CONVERTER_OPT_PRIMARIES_MODE   "GstVideoConverter.primaries-mode"

/**
 * GST_VIDEO_CONVERTER_OPT_THREADS:
 *
 * #G_TYPE_UINT, maximum number of threads to use. Default 1, 0 for the number
 * of cores.
 */
#define GST_VIDEO_CONVERTER_OPT_THREADS   "GstVideoConverter.threads"

/**
 * GST_VIDEO_CONVERTER_OPT_ASYNC_TASKS:
 *
 * #G_TYPE_BOOLEAN, whether gst_video_converter_frame() will return immediately
 * without waiting for the conversion to complete.  A subsequent
 * gst_video_converter_frame_finish() must be performed to ensure completion of the
 * conversion before subsequent use.  Default %FALSE
 *
 * Since: 1.20
 */
#define GST_VIDEO_CONVERTER_OPT_ASYNC_TASKS   "GstVideoConverter.async-tasks"

typedef struct _GstVideoConverter GstVideoConverter;

GST_VIDEO_API
GstVideoConverter *  gst_video_converter_new            (const GstVideoInfo *in_info,
                                                         const GstVideoInfo *out_info,
                                                         GstStructure *config);

GST_VIDEO_API
GstVideoConverter * gst_video_converter_new_with_pool   (const GstVideoInfo * in_info,
                                                         const GstVideoInfo * out_info,
                                                         GstStructure * config,
                                                         GstTaskPool  * pool);

GST_VIDEO_API
void                 gst_video_converter_free           (GstVideoConverter * convert);

GST_VIDEO_API
gboolean             gst_video_converter_set_config     (GstVideoConverter * convert, GstStructure *config);

GST_VIDEO_API
const GstStructure * gst_video_converter_get_config     (GstVideoConverter * convert);

GST_VIDEO_API
void                 gst_video_converter_frame          (GstVideoConverter * convert,
                                                         const GstVideoFrame *src, GstVideoFrame *dest);
GST_VIDEO_API
void                 gst_video_converter_frame_finish   (GstVideoConverter * convert);

GST_VIDEO_API
const GstVideoInfo * gst_video_converter_get_in_info    (GstVideoConverter * convert);

GST_VIDEO_API
const GstVideoInfo * gst_video_converter_get_out_info   (GstVideoConverter * convert);

G_END_DECLS

#endif /* __GST_VIDEO_CONVERTER_H__ */
