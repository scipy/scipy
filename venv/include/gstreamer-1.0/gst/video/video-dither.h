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

#ifndef __GST_VIDEO_DITHER_H__
#define __GST_VIDEO_DITHER_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

/**
 * GstVideoDitherMethod:
 * @GST_VIDEO_DITHER_NONE: no dithering
 * @GST_VIDEO_DITHER_VERTERR: propagate rounding errors downwards
 * @GST_VIDEO_DITHER_FLOYD_STEINBERG: Dither with floyd-steinberg error diffusion
 * @GST_VIDEO_DITHER_SIERRA_LITE: Dither with Sierra Lite error diffusion
 * @GST_VIDEO_DITHER_BAYER: ordered dither using a bayer pattern
 *
 * Different dithering methods to use.
 */
typedef enum {
  GST_VIDEO_DITHER_NONE,
  GST_VIDEO_DITHER_VERTERR,
  GST_VIDEO_DITHER_FLOYD_STEINBERG,
  GST_VIDEO_DITHER_SIERRA_LITE,
  GST_VIDEO_DITHER_BAYER,
} GstVideoDitherMethod;

/**
 * GstVideoDitherFlags:
 * @GST_VIDEO_DITHER_FLAG_NONE: no flags
 * @GST_VIDEO_DITHER_FLAG_INTERLACED: the input is interlaced
 * @GST_VIDEO_DITHER_FLAG_QUANTIZE: quantize values in addition to adding dither.
 *
 * Extra flags that influence the result from gst_video_chroma_resample_new().
 */
typedef enum {
  GST_VIDEO_DITHER_FLAG_NONE       = 0,
  GST_VIDEO_DITHER_FLAG_INTERLACED = (1 << 0),
  GST_VIDEO_DITHER_FLAG_QUANTIZE   = (1 << 1),
} GstVideoDitherFlags;

typedef struct _GstVideoDither GstVideoDither;

/* circular dependency, need to include this after defining the enums */
#include <gst/video/video-format.h>

GST_VIDEO_API
GstVideoDither    * gst_video_dither_new      (GstVideoDitherMethod method,
                                               GstVideoDitherFlags flags,
                                               GstVideoFormat format,
                                               guint quantizer[GST_VIDEO_MAX_COMPONENTS],
                                               guint width);

GST_VIDEO_API
void                gst_video_dither_free     (GstVideoDither *dither);

GST_VIDEO_API
void                gst_video_dither_line     (GstVideoDither *dither,
                                               gpointer line, guint x, guint y, guint width);

G_END_DECLS

#endif /* __GST_VIDEO_DITHER_H__ */
