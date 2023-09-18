/* GStreamer
 * Copyright (C) <2014> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_VIDEO_RESAMPLER_H__
#define __GST_VIDEO_RESAMPLER_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

typedef struct _GstVideoResampler GstVideoResampler;

/**
 * GstVideoResamplerMethod:
 * @GST_VIDEO_RESAMPLER_METHOD_NEAREST: Duplicates the samples when
 *    upsampling and drops when downsampling
 * @GST_VIDEO_RESAMPLER_METHOD_LINEAR: Uses linear interpolation to reconstruct
 *    missing samples and averaging to downsample
 * @GST_VIDEO_RESAMPLER_METHOD_CUBIC: Uses cubic interpolation
 * @GST_VIDEO_RESAMPLER_METHOD_SINC: Uses sinc interpolation
 * @GST_VIDEO_RESAMPLER_METHOD_LANCZOS: Uses lanczos interpolation
 *
 * Different subsampling and upsampling methods
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_RESAMPLER_METHOD_NEAREST,
  GST_VIDEO_RESAMPLER_METHOD_LINEAR,
  GST_VIDEO_RESAMPLER_METHOD_CUBIC,
  GST_VIDEO_RESAMPLER_METHOD_SINC,
  GST_VIDEO_RESAMPLER_METHOD_LANCZOS
} GstVideoResamplerMethod;

/**
 * GST_VIDEO_RESAMPLER_OPT_CUBIC_B:
 *
 * G_TYPE_DOUBLE, B parameter of the cubic filter. The B
 * parameter controls the bluriness. Values between 0.0 and
 * 2.0 are accepted. 1/3 is the default.
 *
 * Below are some values of popular filters:
 *                    B       C
 * Hermite           0.0     0.0
 * Spline            1.0     0.0
 * Catmull-Rom       0.0     1/2
 * Mitchell          1/3     1/3
 * Robidoux          0.3782  0.3109
 * Robidoux
 *  Sharp            0.2620  0.3690
 * Robidoux
 *  Soft             0.6796  0.1602
 */
#define GST_VIDEO_RESAMPLER_OPT_CUBIC_B      "GstVideoResampler.cubic-b"
/**
 * GST_VIDEO_RESAMPLER_OPT_CUBIC_C:
 *
 * G_TYPE_DOUBLE, C parameter of the cubic filter. The C
 * parameter controls the Keys alpha value. Values between 0.0 and
 * 2.0 are accepted. 1/3 is the default.
 *
 * See #GST_VIDEO_RESAMPLER_OPT_CUBIC_B for some more common values
 */
#define GST_VIDEO_RESAMPLER_OPT_CUBIC_C      "GstVideoResampler.cubic-c"

/**
 * GST_VIDEO_RESAMPLER_OPT_ENVELOPE:
 *
 * G_TYPE_DOUBLE, specifies the size of filter envelope for
 * @GST_VIDEO_RESAMPLER_METHOD_LANCZOS. values are clamped between
 * 1.0 and 5.0. 2.0 is the default.
 */
#define GST_VIDEO_RESAMPLER_OPT_ENVELOPE      "GstVideoResampler.envelope"

/**
 * GST_VIDEO_RESAMPLER_OPT_SHARPNESS:
 *
 * G_TYPE_DOUBLE, specifies sharpness of the filter for
 * @GST_VIDEO_RESAMPLER_METHOD_LANCZOS. values are clamped between
 * 0.5 and 1.5. 1.0 is the default.
 */
#define GST_VIDEO_RESAMPLER_OPT_SHARPNESS     "GstVideoResampler.sharpness"

/**
 * GST_VIDEO_RESAMPLER_OPT_SHARPEN:
 *
 * G_TYPE_DOUBLE, specifies sharpening of the filter for
 * @GST_VIDEO_RESAMPLER_METHOD_LANCZOS. values are clamped between
 * 0.0 and 1.0. 0.0 is the default.
 */
#define GST_VIDEO_RESAMPLER_OPT_SHARPEN      "GstVideoResampler.sharpen"
/**
 * GST_VIDEO_RESAMPLER_OPT_MAX_TAPS:
 *
 * G_TYPE_INT, limits the maximum number of taps to use.
 * 16 is the default.
 */
#define GST_VIDEO_RESAMPLER_OPT_MAX_TAPS     "GstVideoResampler.max-taps"

/**
 * GstVideoResamplerFlags:
 * @GST_VIDEO_RESAMPLER_FLAG_NONE: no flags
 * @GST_VIDEO_RESAMPLER_FLAG_HALF_TAPS: when no taps are given, half the
 *              number of calculated taps. This can be used when making scalers
 *              for the different fields of an interlaced picture. Since: 1.10
 *
 * Different resampler flags.
 *
 * Since: 1.6
 */
typedef enum {
  GST_VIDEO_RESAMPLER_FLAG_NONE                 = (0),
  GST_VIDEO_RESAMPLER_FLAG_HALF_TAPS            = (1 << 0),
} GstVideoResamplerFlags;

/**
 * GstVideoResampler:
 * @in_size: the input size
 * @out_size: the output size
 * @max_taps: the maximum number of taps
 * @n_phases: the number of phases
 * @offset: array with the source offset for each output element
 * @phase: array with the phase to use for each output element
 * @n_taps: array with new number of taps for each phase
 * @taps: the taps for all phases
 *
 * A structure holding resampler information.
 *
 * Since: 1.6
 */
struct _GstVideoResampler
{
  gint in_size;
  gint out_size;
  guint max_taps;
  guint n_phases;
  guint32 *offset;
  guint32 *phase;
  guint32 *n_taps;
  gdouble *taps;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};


GST_VIDEO_API
gboolean    gst_video_resampler_init      (GstVideoResampler *resampler,
                                           GstVideoResamplerMethod method,
                                           GstVideoResamplerFlags flags,
                                           guint n_phases, guint n_taps,
                                           gdouble shift,
                                           guint in_size, guint out_size,
                                           GstStructure *options);

GST_VIDEO_API
void        gst_video_resampler_clear     (GstVideoResampler *resampler);

G_END_DECLS

#endif /* __GST_VIDEO_RESAMPLER_H__ */
