/* GStreamer
 * Copyright (C) <2011> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_VIDEO_POOL_H__
#define __GST_VIDEO_POOL_H__

#include <gst/gst.h>

#include <gst/video/video.h>

G_BEGIN_DECLS

/**
 * GST_BUFFER_POOL_OPTION_VIDEO_META:
 *
 * An option that can be activated on bufferpool to request video metadata
 * on buffers from the pool.
 */
#define GST_BUFFER_POOL_OPTION_VIDEO_META "GstBufferPoolOptionVideoMeta"

/**
 * GST_BUFFER_POOL_OPTION_VIDEO_ALIGNMENT:
 *
 * A bufferpool option to enable extra padding. When a bufferpool supports this
 * option, gst_buffer_pool_config_set_video_alignment() can be called.
 *
 * When this option is enabled on the bufferpool,
 * #GST_BUFFER_POOL_OPTION_VIDEO_META should also be enabled.
 */
#define GST_BUFFER_POOL_OPTION_VIDEO_ALIGNMENT "GstBufferPoolOptionVideoAlignment"

/* setting a bufferpool config */

GST_VIDEO_API
void             gst_buffer_pool_config_set_video_alignment  (GstStructure *config, const GstVideoAlignment *align);

GST_VIDEO_API
gboolean         gst_buffer_pool_config_get_video_alignment  (GstStructure *config, GstVideoAlignment *align);

/* video bufferpool */
typedef struct _GstVideoBufferPool GstVideoBufferPool;
typedef struct _GstVideoBufferPoolClass GstVideoBufferPoolClass;
typedef struct _GstVideoBufferPoolPrivate GstVideoBufferPoolPrivate;

#define GST_TYPE_VIDEO_BUFFER_POOL      (gst_video_buffer_pool_get_type())
#define GST_IS_VIDEO_BUFFER_POOL(obj)   (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_VIDEO_BUFFER_POOL))
#define GST_VIDEO_BUFFER_POOL(obj)      (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_VIDEO_BUFFER_POOL, GstVideoBufferPool))
#define GST_VIDEO_BUFFER_POOL_CAST(obj) ((GstVideoBufferPool*)(obj))

struct _GstVideoBufferPool
{
  GstBufferPool bufferpool;

  GstVideoBufferPoolPrivate *priv;
};

struct _GstVideoBufferPoolClass
{
  GstBufferPoolClass parent_class;
};

GST_VIDEO_API
GType             gst_video_buffer_pool_get_type      (void);

GST_VIDEO_API
GstBufferPool *   gst_video_buffer_pool_new           (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstVideoBufferPool, gst_object_unref)

G_END_DECLS

#endif /* __GST_VIDEO_POOL_H__ */
