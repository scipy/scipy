/*
 * GStreamer
 * Copyright (C) 2012 Matthew Waters <ystreet00@gmail.com>
 * Copyright (C) 2015 Matthew Waters <matthew@centricular.com>
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

#ifndef _GST_GL_MEMORY_PBO_H_
#define _GST_GL_MEMORY_PBO_H_

#include <gst/gl/gstglmemory.h>

G_BEGIN_DECLS

#define GST_TYPE_GL_MEMORY_PBO_ALLOCATOR (gst_gl_memory_pbo_allocator_get_type())
GST_GL_API
GType gst_gl_memory_pbo_allocator_get_type(void);

#define GST_IS_GL_MEMORY_PBO_ALLOCATOR(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_GL_MEMORY_PBO_ALLOCATOR))
#define GST_IS_GL_MEMORY_PBO_ALLOCATOR_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_GL_MEMORY_PBO_ALLOCATOR))
#define GST_GL_MEMORY_PBO_ALLOCATOR_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_GL_MEMORY_PBO_ALLOCATOR, GstGLMemoryPBOAllocatorClass))
#define GST_GL_MEMORY_PBO_ALLOCATOR(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_GL_MEMORY_PBO_ALLOCATOR, GstGLMemoryPBOAllocator))
#define GST_GL_MEMORY_PBO_ALLOCATOR_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_GL_MEMORY_PBO_ALLOCATOR, GstGLMemoryPBOAllocatorClass))
#define GST_GL_MEMORY_PBO_ALLOCATOR_CAST(obj)            ((GstGLMemoryPBOAllocator *)(obj))

/**
 * GstGLMemoryPBO:
 *
 * Private instance
 */
struct _GstGLMemoryPBO
{
  /*< private >*/
  GstGLMemory      mem;

  GstGLBuffer          *pbo;

  gpointer                  _padding[GST_PADDING];
};

/**
 * GST_GL_MEMORY_PBO_ALLOCATOR_NAME:
 *
 * The name of the GL Memory PBO allocator
 */
#define GST_GL_MEMORY_PBO_ALLOCATOR_NAME   "GLMemoryPBO"

/**
 * GST_TYPE_GL_MEMORY_PBO
 *
 * Since: 1.20
 * Deprecated: 1.22: This type has no use.
 */
#define GST_TYPE_GL_MEMORY_PBO (gst_gl_memory_pbo_get_type())
GST_GL_DEPRECATED
GType gst_gl_memory_pbo_get_type(void);

GST_GL_API
void          gst_gl_memory_pbo_init_once               (void);
GST_GL_API
gboolean      gst_is_gl_memory_pbo                      (GstMemory * mem);

GST_GL_API
void          gst_gl_memory_pbo_download_transfer       (GstGLMemoryPBO * gl_mem);
GST_GL_API
void          gst_gl_memory_pbo_upload_transfer         (GstGLMemoryPBO * gl_mem);

GST_GL_API
gboolean      gst_gl_memory_pbo_copy_into_texture       (GstGLMemoryPBO *gl_mem,
                                                         guint tex_id,
                                                         GstGLTextureTarget target,
                                                         GstGLFormat tex_format,
                                                         gint width,
                                                         gint height,
                                                         gint stride,
                                                         gboolean respecify);

/**
 * GstGLMemoryPBOAllocator:
 *
 * Opaque #GstGLMemoryPBOAllocator struct
 */
struct _GstGLMemoryPBOAllocator
{
  GstGLMemoryAllocator parent;

  /*< private >*/
  gpointer             _padding[GST_PADDING];
};

/**
 * GstGLMemoryPBOAllocatorClass:
 *
 * Only contains private data
 */
struct _GstGLMemoryPBOAllocatorClass
{
  GstGLMemoryAllocatorClass parent_class;

  /*< private >*/
  gpointer                  _padding[GST_PADDING];
};

G_END_DECLS

#endif /* _GST_GL_MEMORY_PBO_H_ */
