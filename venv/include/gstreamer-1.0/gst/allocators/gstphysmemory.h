/* GStreamer
 * Copyright (C) 2017 Sebastian Dröge <sebastian@centricular.com>
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

#ifndef __GST_PHYS_MEMORY_H__
#define __GST_PHYS_MEMORY_H__

#include <gst/gst.h>
#include <gst/allocators/allocators-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_PHYS_MEMORY_ALLOCATOR (gst_phys_memory_allocator_get_type())
GST_ALLOCATORS_API
G_DECLARE_INTERFACE (GstPhysMemoryAllocator, gst_phys_memory_allocator,
    GST, PHYS_MEMORY_ALLOCATOR, GstAllocator)

#define GST_PHYS_MEMORY_ALLOCATOR_GET_INTERFACE(obj) (GST_PHYS_MEMORY_ALLOCATOR_GET_IFACE(obj))
#define GST_PHYS_MEMORY_ALLOCATOR_CAST(obj) ((GstPhysMemoryAllocator *)(obj))

/**
 * GstPhysMemoryAllocatorInterface:
 * @get_phys_addr: Implementations shall return the physicall memory address
 *    that is backing the provided memory, or 0 if none.
 *
 * Marker interface for allocators with physical address backed memory
 *
 * Since: 1.14
 */
struct _GstPhysMemoryAllocatorInterface
{
  /*< private >*/
  GTypeInterface parent_iface;

  /*< public >*/
  guintptr (*get_phys_addr) (GstPhysMemoryAllocator * allocator, GstMemory * mem);
};

GST_ALLOCATORS_API
gboolean gst_is_phys_memory (GstMemory *mem);

GST_ALLOCATORS_API
guintptr gst_phys_memory_get_phys_addr (GstMemory * mem);

G_END_DECLS

#endif /* __GST_PHYS_MEMORY_H__ */
