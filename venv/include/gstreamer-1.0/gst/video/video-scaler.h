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

#ifndef __GST_VIDEO_SCALER_H__
#define __GST_VIDEO_SCALER_H__

#include <gst/gst.h>

#include <gst/video/video-format.h>
#include <gst/video/video-color.h>
#include <gst/video/video-resampler.h>

G_BEGIN_DECLS

/**
 * GST_VIDEO_SCALER_OPT_DITHER_METHOD:
 *
 * #GstVideoDitherMethod, The dither method to use for propagating
 * quatization errors.
 */
#define GST_VIDEO_SCALER_OPT_DITHER_METHOD   "GstVideoScaler.dither-method"

/**
 * GstVideoScalerFlags:
 * @GST_VIDEO_SCALER_FLAG_NONE: no flags
 * @GST_VIDEO_SCALER_FLAG_INTERLACED: Set up a scaler for interlaced content
 *
 * Different scale flags.
 */
typedef enum {
  GST_VIDEO_SCALER_FLAG_NONE                 = (0),
  GST_VIDEO_SCALER_FLAG_INTERLACED           = (1 << 0),
} GstVideoScalerFlags;

typedef struct _GstVideoScaler GstVideoScaler;

GST_VIDEO_API
GstVideoScaler *      gst_video_scaler_new            (GstVideoResamplerMethod method,
                                                       GstVideoScalerFlags flags,
                                                       guint n_taps,
                                                       guint in_size, guint out_size,
                                                       GstStructure * options);

GST_VIDEO_API
void                  gst_video_scaler_free           (GstVideoScaler *scale);

GST_VIDEO_API
guint                 gst_video_scaler_get_max_taps   (GstVideoScaler *scale);

GST_VIDEO_API
const gdouble *       gst_video_scaler_get_coeff      (GstVideoScaler *scale,
                                                       guint out_offset,
                                                       guint *in_offset,
                                                       guint *n_taps);

GST_VIDEO_API
void                  gst_video_scaler_horizontal     (GstVideoScaler *scale,
                                                       GstVideoFormat format,
                                                       gpointer src, gpointer dest,
                                                       guint dest_offset, guint width);

GST_VIDEO_API
void                  gst_video_scaler_vertical       (GstVideoScaler *scale,
                                                       GstVideoFormat format,
                                                       gpointer src_lines[], gpointer dest,
                                                       guint dest_offset, guint width);

GST_VIDEO_API
GstVideoScaler *      gst_video_scaler_combine_packed_YUV (GstVideoScaler * y_scale,
                                                           GstVideoScaler *uv_scale,
                                                           GstVideoFormat in_format,
                                                           GstVideoFormat out_format);

GST_VIDEO_API
void                  gst_video_scaler_2d             (GstVideoScaler *hscale,
                                                       GstVideoScaler *vscale,
                                                       GstVideoFormat format,
                                                       gpointer src, gint src_stride,
                                                       gpointer dest, gint dest_stride,
                                                       guint x, guint y,
                                                       guint width, guint height);

G_END_DECLS

#endif /* __GST_VIDEO_SCALER_H__ */
