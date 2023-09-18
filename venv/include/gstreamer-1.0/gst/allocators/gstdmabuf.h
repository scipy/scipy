/* GStreamer dmabuf allocator
 * Copyright (C) 2013 Linaro SA
 * Author: Benjamin Gaignard <benjamin.gaignard@linaro.org> for Linaro.
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __GST_DMABUF_H__
#define __GST_DMABUF_H__

#include <gst/gst.h>
#include <gst/allocators/gstfdmemory.h>

G_BEGIN_DECLS

/**
 * GST_CAPS_FEATURE_MEMORY_DMABUF:
 *
 * Constant that defines the caps feature name for DMA buffer sharing.
 *
 * It has to be used for non-mappable dma-buf only, i.e. when the underlying
 * memory is not mappable to user space. Or when the mapped memory contains
 * non meaningful data. It can be the case for protected content or when the
 * user wants explicitly avoid any software post processing.
 *
 * In these cases all elements between the exported and the importer has to work
 * in passthrough mode. This is done by adding this caps feature.
 *
 * When the memory is mappable for read and write requests then it is assumes
 * to be a fast path and so this caps feature should not be used. Though
 * according to the dma-buf protocol, while it is mapped it prevents the
 * exporter to migrate the buffer.
 *
 * This caps feature should not serve at all the purpose of selecting the
 * @GST_ALLOCATOR_DMABUF allocator during caps negotiation.
 * When the exporter is the upstream element from the importer point of view,
 * the exporter should try to map the dma buffer at runtime (preferably during
 * decide_allocation phase). When it succeeds for #GST_MAP_READWRITE this caps
 * feature should not be used. This allows scalers, color converts and any image
 * processing filters to work directly on the dma buffer.
 * In this case the importer element should check all incoming memory using
 * gst_is_dmabuf_memory().
 *
 * Since: 1.12
 */
#define GST_CAPS_FEATURE_MEMORY_DMABUF "memory:DMABuf"

#define GST_ALLOCATOR_DMABUF "dmabuf"

#define GST_TYPE_DMABUF_ALLOCATOR              (gst_dmabuf_allocator_get_type())
#define GST_IS_DMABUF_ALLOCATOR(obj)           (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_DMABUF_ALLOCATOR))
#define GST_IS_DMABUF_ALLOCATOR_CLASS(klass)   (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_DMABUF_ALLOCATOR))
#define GST_DMABUF_ALLOCATOR_GET_CLASS(obj)    (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_DMABUF_ALLOCATOR, GstDmaBufAllocatorClass))
#define GST_DMABUF_ALLOCATOR(obj)              (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_DMABUF_ALLOCATOR, GstDmaBufAllocator))
#define GST_DMABUF_ALLOCATOR_CLASS(klass)      (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_DMABUF_ALLOCATOR, GstDmaBufAllocatorClass))
#define GST_DMABUF_ALLOCATOR_CAST(obj)         ((GstDmaBufAllocator *)(obj))

typedef struct _GstDmaBufAllocator GstDmaBufAllocator;
typedef struct _GstDmaBufAllocatorClass GstDmaBufAllocatorClass;

/**
 * GstDmaBufAllocator:
 *
 * Base class for allocators with dmabuf-backed memory
 *
 * Since: 1.12
 */
struct _GstDmaBufAllocator
{
  GstFdAllocator parent;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

struct _GstDmaBufAllocatorClass
{
  GstFdAllocatorClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};


GST_ALLOCATORS_API
GType          gst_dmabuf_allocator_get_type (void);

GST_ALLOCATORS_API
GstAllocator * gst_dmabuf_allocator_new (void);

GST_ALLOCATORS_API
GstMemory    * gst_dmabuf_allocator_alloc (GstAllocator * allocator, gint fd, gsize size);

GST_ALLOCATORS_API
GstMemory    * gst_dmabuf_allocator_alloc_with_flags (GstAllocator * allocator, gint fd, gsize size, GstFdMemoryFlags flags);

GST_ALLOCATORS_API
gint           gst_dmabuf_memory_get_fd (GstMemory * mem);

GST_ALLOCATORS_API
gboolean       gst_is_dmabuf_memory (GstMemory * mem);


G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstDmaBufAllocator, gst_object_unref)

G_END_DECLS
#endif /* __GST_DMABUF_H__ */
