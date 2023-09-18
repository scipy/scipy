/*
 * GStreamer
 * Copyright (C) 2012 Matthew Waters <ystreet00@gmail.com>
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

#ifndef _GST_GL_BUFFER_POOL_H_
#define _GST_GL_BUFFER_POOL_H_

#include <gst/video/gstvideometa.h>
#include <gst/video/gstvideopool.h>

#include <gst/gl/gstglbasememory.h>

G_BEGIN_DECLS

/* buffer pool functions */
GST_GL_API
GType gst_gl_buffer_pool_get_type (void);
#define GST_TYPE_GL_BUFFER_POOL      (gst_gl_buffer_pool_get_type())
#define GST_IS_GL_BUFFER_POOL(obj)   (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_GL_BUFFER_POOL))
#define GST_GL_BUFFER_POOL(obj)      (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_GL_BUFFER_POOL, GstGLBufferPool))
#define GST_GL_BUFFER_POOL_CAST(obj) ((GstGLBufferPool*)(obj))

/**
 * GstGLBufferPool:
 *
 * Opaque GstGLBufferPool struct
 */
struct _GstGLBufferPool
{
  GstBufferPool bufferpool;

  GstGLContext *context;

  /*< private >*/
  GstGLBufferPoolPrivate *priv;

  gpointer _padding[GST_PADDING];
};

/**
 * GstGLBufferPoolClass:
 *
 * The #GstGLBufferPoolClass structure contains only private data
 */
struct _GstGLBufferPoolClass
{
  GstBufferPoolClass parent_class;

  /*< private >*/
  gpointer _padding[GST_PADDING];
};

GST_GL_API
GstBufferPool *gst_gl_buffer_pool_new (GstGLContext * context);

GST_GL_API
GstGLAllocationParams * gst_gl_buffer_pool_get_gl_allocation_params        (GstGLBufferPool * pool);

GST_GL_API
GstGLAllocationParams * gst_buffer_pool_config_get_gl_allocation_params    (GstStructure * config);
GST_GL_API
void                    gst_buffer_pool_config_set_gl_allocation_params    (GstStructure * config,
                                                                            const GstGLAllocationParams * params);

G_END_DECLS

#endif /* _GST_GL_BUFFER_POOL_H_ */
