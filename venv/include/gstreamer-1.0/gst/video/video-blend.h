/* Gstreamer video blending utility functions
 *
 * Copyright (C) <2011> Intel Corporation
 * Copyright (C) <2011> Collabora Ltd.
 * Copyright (C) <2011> Thibault Saunier <thibault.saunier@collabora.com>
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


#ifndef  __GST_VIDEO_BLEND__
#define  __GST_VIDEO_BLEND__

#include <gst/gst.h>
#include <gst/video/video.h>

GST_VIDEO_API
void       gst_video_blend_scale_linear_RGBA  (GstVideoInfo * src, GstBuffer * src_buffer,
                                               gint dest_height, gint dest_width,
                                               GstVideoInfo * dest, GstBuffer ** dest_buffer);

GST_VIDEO_API
gboolean   gst_video_blend                    (GstVideoFrame * dest,
                                               GstVideoFrame * src,
                                               gint x, gint y,
                                               gfloat global_alpha);

#endif
