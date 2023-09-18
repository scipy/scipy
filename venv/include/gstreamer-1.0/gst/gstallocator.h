/* GStreamer
 * Copyright (C) 2009 Wim Taymans <wim.taymans@gmail.be>
 *
 * gstallocator.h: Header for memory allocation
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


#ifndef __GST_ALLOCATOR_H__
#define __GST_ALLOCATOR_H__

#include <gst/gstmemory.h>
#include <gst/gstobject.h>

G_BEGIN_DECLS

typedef struct _GstAllocatorPrivate GstAllocatorPrivate;
typedef struct _GstAllocatorClass GstAllocatorClass;

#define GST_TYPE_ALLOCATOR                 (gst_allocator_get_type())
#define GST_IS_ALLOCATOR(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_ALLOCATOR))
#define GST_IS_ALLOCATOR_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_ALLOCATOR))
#define GST_ALLOCATOR_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_ALLOCATOR, GstAllocatorClass))
#define GST_ALLOCATOR(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_ALLOCATOR, GstAllocator))
#define GST_ALLOCATOR_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_ALLOCATOR, GstAllocatorClass))
#define GST_ALLOCATOR_CAST(obj)            ((GstAllocator *)(obj))

#define GST_TYPE_ALLOCATION_PARAMS (gst_allocation_params_get_type())

GST_API
GType gst_allocation_params_get_type(void);

typedef struct _GstAllocationParams GstAllocationParams;

/**
 * gst_memory_alignment:
 *
 * The default memory alignment in bytes - 1
 * an alignment of 7 would be the same as what malloc() guarantees.
 */

GST_API gsize gst_memory_alignment;

/**
 * GST_ALLOCATOR_SYSMEM:
 *
 * The allocator name for the default system memory allocator
 */
#define GST_ALLOCATOR_SYSMEM   "SystemMemory"

/**
 * GstAllocationParams:
 * @flags: flags to control allocation
 * @align: the desired alignment of the memory
 * @prefix: the desired prefix
 * @padding: the desired padding
 *
 * Parameters to control the allocation of memory
 */
struct _GstAllocationParams {
  GstMemoryFlags flags;
  gsize          align;
  gsize          prefix;
  gsize          padding;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAllocatorFlags:
 * @GST_ALLOCATOR_FLAG_CUSTOM_ALLOC: The allocator has a custom alloc function.
 * @GST_ALLOCATOR_FLAG_LAST: first flag that can be used for custom purposes
 *
 * Flags for allocators.
 */
typedef enum {
  GST_ALLOCATOR_FLAG_CUSTOM_ALLOC  = (GST_OBJECT_FLAG_LAST << 0),

  GST_ALLOCATOR_FLAG_LAST          = (GST_OBJECT_FLAG_LAST << 16)
} GstAllocatorFlags;

/**
 * GstAllocator:
 * @mem_map: the implementation of the GstMemoryMapFunction
 * @mem_unmap: the implementation of the GstMemoryUnmapFunction
 * @mem_copy: the implementation of the GstMemoryCopyFunction
 * @mem_share: the implementation of the GstMemoryShareFunction
 * @mem_is_span: the implementation of the GstMemoryIsSpanFunction
 * @mem_map_full: the implementation of the GstMemoryMapFullFunction.
 *      Will be used instead of @mem_map if present. (Since: 1.6)
 * @mem_unmap_full: the implementation of the GstMemoryUnmapFullFunction.
 *      Will be used instead of @mem_unmap if present. (Since: 1.6)
 *
 * The #GstAllocator is used to create new memory.
 */
struct _GstAllocator
{
  GstObject  object;

  const gchar               *mem_type;

  /*< public >*/
  GstMemoryMapFunction       mem_map;
  GstMemoryUnmapFunction     mem_unmap;

  GstMemoryCopyFunction      mem_copy;
  GstMemoryShareFunction     mem_share;
  GstMemoryIsSpanFunction    mem_is_span;

  GstMemoryMapFullFunction   mem_map_full;
  GstMemoryUnmapFullFunction mem_unmap_full;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING - 2];

  GstAllocatorPrivate *priv;
};

/**
 * GstAllocatorClass:
 * @object_class:  Object parent class
 * @alloc: implementation that acquires memory
 * @free: implementation that releases memory
 *
 * The #GstAllocator is used to create new memory.
 */
struct _GstAllocatorClass {
  GstObjectClass object_class;

  /*< public >*/
  GstMemory *  (*alloc)      (GstAllocator *allocator, gsize size,
                              GstAllocationParams *params);
  void         (*free)       (GstAllocator *allocator, GstMemory *memory);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType          gst_allocator_get_type (void);

/* allocators */

GST_API
void           gst_allocator_register        (const gchar *name, GstAllocator *allocator);

GST_API
GstAllocator * gst_allocator_find            (const gchar *name);

GST_API
void           gst_allocator_set_default     (GstAllocator * allocator);

/* allocation parameters */

GST_API
GstAllocationParams * gst_allocation_params_new (void) G_GNUC_MALLOC;

GST_API
void           gst_allocation_params_init    (GstAllocationParams *params);

GST_API
GstAllocationParams *
               gst_allocation_params_copy    (const GstAllocationParams *params) G_GNUC_MALLOC;

GST_API
void           gst_allocation_params_free    (GstAllocationParams *params);

/* allocating memory blocks */

GST_API
GstMemory *    gst_allocator_alloc           (GstAllocator * allocator, gsize size,
                                              GstAllocationParams *params);

GST_API
void           gst_allocator_free            (GstAllocator * allocator, GstMemory *memory);

GST_API
GstMemory *    gst_memory_new_wrapped  (GstMemoryFlags flags, gpointer data, gsize maxsize,
                                        gsize offset, gsize size, gpointer user_data,
                                        GDestroyNotify notify);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAllocationParams, gst_allocation_params_free)

G_END_DECLS

#endif /* __GST_ALLOCATOR_H__ */
