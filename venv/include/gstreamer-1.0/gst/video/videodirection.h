/* GStreamer
 * Copyright (C) 2016 Igalia <calvaris@igalia.com>
 *
 * videodirection.h: video rotation and flipping interface
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

#ifndef __GST_VIDEO_DIRECTION_H__
#define __GST_VIDEO_DIRECTION_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS
#define GST_TYPE_VIDEO_DIRECTION \
  (gst_video_direction_get_type ())
#define GST_VIDEO_DIRECTION(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_VIDEO_DIRECTION, GstVideoDirection))
#define GST_IS_VIDEO_DIRECTION(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_VIDEO_DIRECTION))
#define GST_VIDEO_DIRECTION_GET_INTERFACE(inst) \
  (G_TYPE_INSTANCE_GET_INTERFACE ((inst), GST_TYPE_VIDEO_DIRECTION, GstVideoDirectionInterface))
/**
 * GstVideoDirection:
 *
 * Opaque #GstVideoDirection data structure.
 *
 * Since: 1.10
 */
typedef struct _GstVideoDirection GstVideoDirection;
typedef struct _GstVideoDirectionInterface GstVideoDirectionInterface;

/**
 * GstVideoDirectionInterface:
 * @iface: parent interface type.
 *
 * #GstVideoDirectionInterface interface.
 *
 * Since: 1.10
 */
struct _GstVideoDirectionInterface
{
  GTypeInterface iface;
};

GST_VIDEO_API
GType gst_video_direction_get_type (void);

G_END_DECLS
#endif /* __GST_VIDEO_DIRECTION_H__ */
