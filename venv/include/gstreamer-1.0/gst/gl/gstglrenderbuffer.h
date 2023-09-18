/*
 * GStreamer
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

#ifndef _GST_GL_RENDERBUFFER_H_
#define _GST_GL_RENDERBUFFER_H_

#include <gst/gl/gstglbasememory.h>

G_BEGIN_DECLS

#define GST_TYPE_GL_RENDERBUFFER_ALLOCATOR (gst_gl_renderbuffer_allocator_get_type())
GST_GL_API GType gst_gl_renderbuffer_allocator_get_type(void);

#define GST_IS_GL_RENDERBUFFER_ALLOCATOR(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_GL_RENDERBUFFER_ALLOCATOR))
#define GST_IS_GL_RENDERBUFFER_ALLOCATOR_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_GL_RENDERBUFFER_ALLOCATOR))
#define GST_GL_RENDERBUFFER_ALLOCATOR_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_GL_RENDERBUFFER_ALLOCATOR, GstGLRenderbufferAllocatorClass))
#define GST_GL_RENDERBUFFER_ALLOCATOR(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_GL_RENDERBUFFER_ALLOCATOR, GstGLRenderbufferAllocator))
#define GST_GL_RENDERBUFFER_ALLOCATOR_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_GL_RENDERBUFFER_ALLOCATOR, GstGLRenderbufferAllocatorClass))
#define GST_GL_RENDERBUFFER_ALLOCATOR_CAST(obj)            ((GstGLRenderbufferAllocator *)(obj))

#define GST_GL_RENDERBUFFER_CAST(obj) ((GstGLRenderbuffer *) obj)

/**
 * GST_GL_RENDERBUFFER_ALLOCATOR_NAME:
 *
 * The name of the GL renderbuffer allocator
 */
#define GST_GL_RENDERBUFFER_ALLOCATOR_NAME   "GLRenderbuffer"

/**
 * GstGLRenderbuffer:
 * @renderbuffer_id: the GL texture id for this memory
 * @renderbuffer_format: the texture type
 * @width: the width
 * @height: the height
 *
 * Represents information about a GL renderbuffer
 */
struct _GstGLRenderbuffer
{
  /*< private >*/
  GstGLBaseMemory           mem;

  /*< public >*/
  guint                     renderbuffer_id;
  GstGLFormat               renderbuffer_format;
  guint                     width;
  guint                     height;

  /*< protected >*/
  gboolean                  renderbuffer_wrapped;

  /*< private >*/
  gpointer                  _padding[GST_PADDING];
};

/**
 * GstGLRenderbufferAllocator:
 *
 * Opaque #GstGLRenderbufferAllocator struct
 */
struct _GstGLRenderbufferAllocator
{
  GstGLBaseMemoryAllocator parent;

  /*< private >*/
  gpointer                  _padding[GST_PADDING];
};

/**
 * GstGLRenderbufferAllocatorClass:
 *
 * The #GstGLRenderbufferAllocatorClass only contains private data
 */
struct _GstGLRenderbufferAllocatorClass
{
  GstGLBaseMemoryAllocatorClass             parent_class;

  /*< private >*/
  gpointer                  _padding[GST_PADDING];
};

#include <gst/gl/gstglbasememory.h>

typedef struct _GstGLRenderbufferAllocationParams GstGLRenderbufferAllocationParams;

GST_GL_API GType gst_gl_renderbuffer_allocation_params_get_type (void);
#define GST_TYPE_RENDERBUFFER_ALLOCATION_PARAMS (gst_gl_renderbuffer_allocation_params_get_type)

/**
 * GstGLRenderbufferAllocationParams:
 * @renderbuffer_format: the #GstGLFormat
 * @width: the width
 * @height: the height
 *
 * Allocation parameters
 */
struct _GstGLRenderbufferAllocationParams
{
  /*< private >*/
  GstGLAllocationParams parent;

  /*< public >*/
  GstGLFormat renderbuffer_format;
  guint width;
  guint height;

  /*< private >*/
  gpointer _padding[GST_PADDING];
};

GST_GL_API
GstGLRenderbufferAllocationParams *     gst_gl_renderbuffer_allocation_params_new           (GstGLContext * context,
                                                                                             const GstAllocationParams * alloc_params,
                                                                                             GstGLFormat renderbuffer_format,
                                                                                             guint width,
                                                                                             guint height);

GST_GL_API
GstGLRenderbufferAllocationParams *     gst_gl_renderbuffer_allocation_params_new_wrapped   (GstGLContext * context,
                                                                                             const GstAllocationParams * alloc_params,
                                                                                             GstGLFormat renderbuffer_format,
                                                                                             guint width,
                                                                                             guint height,
                                                                                             gpointer gl_handle,
                                                                                             gpointer user_data,
                                                                                             GDestroyNotify notify);

/**
 * GST_TYPE_GL_RENDERBUFFER:
 *
 * Since: 1.20
 * Deprecated: 1.22: This type has no use.
 */
#define GST_TYPE_GL_RENDERBUFFER (gst_gl_renderbuffer_get_type())
GST_GL_DEPRECATED
GType gst_gl_renderbuffer_get_type(void);

GST_GL_API
void            gst_gl_renderbuffer_init_once   (void);

GST_GL_API
gboolean        gst_is_gl_renderbuffer          (GstMemory * mem);

/* accessors */
GST_GL_API
gint                    gst_gl_renderbuffer_get_width     (GstGLRenderbuffer * gl_mem);

GST_GL_API
gint                    gst_gl_renderbuffer_get_height    (GstGLRenderbuffer * gl_mem);

GST_GL_API
GstGLFormat             gst_gl_renderbuffer_get_format    (GstGLRenderbuffer * gl_mem);

GST_GL_API
guint                   gst_gl_renderbuffer_get_id        (GstGLRenderbuffer * gl_mem);

G_END_DECLS

#endif /* _GST_GL_RENDERBUFFER_H_ */
