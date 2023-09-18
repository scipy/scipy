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

#ifndef _GST_GL_BASE_MEMORY_H_
#define _GST_GL_BASE_MEMORY_H_

#include <gst/gst.h>
#include <gst/gstallocator.h>
#include <gst/gstmemory.h>

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

/**
 * GST_TYPE_GL_BASE_MEMORY:
 *
 * Deprecated: 1.22: This type has no use.
 */
#define GST_TYPE_GL_BASE_MEMORY (gst_gl_base_memory_get_type())
GST_GL_DEPRECATED
GType gst_gl_base_memory_get_type(void);

#define GST_TYPE_GL_BASE_MEMORY_ALLOCATOR (gst_gl_base_memory_allocator_get_type())
GST_GL_API
GType gst_gl_base_memory_allocator_get_type(void);

#define GST_IS_GL_BASE_MEMORY_ALLOCATOR(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_GL_BASE_MEMORY_ALLOCATOR))
#define GST_IS_GL_BASE_MEMORY_ALLOCATOR_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_GL_BASE_MEMORY_ALLOCATOR))
#define GST_GL_BASE_MEMORY_ALLOCATOR_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_GL_BASE_MEMORY_ALLOCATOR, GstGLBaseMemoryAllocatorClass))
#define GST_GL_BASE_MEMORY_ALLOCATOR(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_GL_BASE_MEMORY_ALLOCATOR, GstGLBaseMemoryAllocator))
#define GST_GL_BASE_MEMORY_ALLOCATOR_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_GL_BASE_MEMORY_ALLOCATOR, GstGLBaseMemoryAllocatorClass))
#define GST_GL_BASE_MEMORY_ALLOCATOR_CAST(obj)            ((GstGLBaseMemoryAllocator *)(obj))

#define GST_GL_BASE_MEMORY_CAST(mem) ((GstGLBaseMemory *)mem)

GST_GL_API
GQuark gst_gl_base_memory_error_quark (void);
/**
 * GST_GL_BASE_MEMORY_ERROR:
 *
 * Error domain for GStreamer's GL memory module. Errors in this domain will be
 * from the #GstGLBaseMemoryError enumeration
 */
#define GST_GL_BASE_MEMORY_ERROR (gst_gl_base_memory_error_quark ())

/**
 * GstGLBaseMemoryError:
 * @GST_GL_BASE_MEMORY_ERROR_FAILED: generic failure
 * @GST_GL_BASE_MEMORY_ERROR_OLD_LIBS: the implementation is too old and doesn't
 *                                     implement enough features
 * @GST_GL_BASE_MEMORY_ERROR_RESOURCE_UNAVAILABLE: a resource could not be found
 */
typedef enum
{
  GST_GL_BASE_MEMORY_ERROR_FAILED,
  GST_GL_BASE_MEMORY_ERROR_OLD_LIBS,
  GST_GL_BASE_MEMORY_ERROR_RESOURCE_UNAVAILABLE,
} GstGLBaseMemoryError;

/**
 * GstGLBaseMemoryTransfer:
 * @GST_GL_BASE_MEMORY_TRANSFER_NEED_DOWNLOAD: the texture needs downloading
 *                                             to the data pointer
 * @GST_GL_BASE_MEMORY_TRANSFER_NEED_UPLOAD:   the data pointer needs uploading
 *                                             to the texture
 */
typedef enum
{
  GST_GL_BASE_MEMORY_TRANSFER_NEED_DOWNLOAD   = (GST_MEMORY_FLAG_LAST << 0),
  GST_GL_BASE_MEMORY_TRANSFER_NEED_UPLOAD     = (GST_MEMORY_FLAG_LAST << 1)
} GstGLBaseMemoryTransfer;

/**
 * GST_MAP_GL:
 *
 * Flag indicating that we should map the GL object instead of to system memory.
 *
 * Combining #GST_MAP_GL with #GST_MAP_WRITE has the same semantics as though
 * you are writing to OpenGL. Conversely, combining #GST_MAP_GL with
 * #GST_MAP_READ has the same semantics as though you are reading from OpenGL.
 */
#define GST_MAP_GL (GST_MAP_FLAG_LAST << 1)

/**
 * GstGLBaseMemory:
 * @mem: the parent object
 * @context: the #GstGLContext to use for GL operations
 *
 * Represents information about a GL memory object
 */
struct _GstGLBaseMemory
{
  GstMemory             mem;

  GstGLContext         *context;

  /*< protected >*/
  GMutex                lock;

  GstMapFlags           map_flags;       /* cumulative map flags */
  gint                  map_count;
  gint                  gl_map_count;

  gpointer              data;

  GstGLQuery           *query;

  /*< private >*/
  gsize                 alloc_size;     /* because maxsize is used for mapping */
  gpointer              alloc_data;

  GDestroyNotify        notify;
  gpointer              user_data;

  gpointer              _padding[GST_PADDING];
};

typedef struct _GstGLAllocationParams GstGLAllocationParams;
/**
 * GstGLAllocationParamsCopyFunc:
 * @src: the source #GstGLAllocationParams to copy from
 * @dest: the source #GstGLAllocationParams to copy
 *
 * Copies the parameters from @src into @dest.  The subclass must compose copy
 * functions from the superclass.
 */
typedef void    (*GstGLAllocationParamsCopyFunc)    (GstGLAllocationParams * src, GstGLAllocationParams * dest);
/**
 * GstGLAllocationParamsFreeFunc:
 * @params: a #GstGLAllocationParams
 *
 * Free any dynamically allocated data.  The subclass must call the superclass'
 * free.
 */
typedef void    (*GstGLAllocationParamsFreeFunc)    (gpointer params);

#define GST_TYPE_GL_ALLOCATION_PARAMS (gst_gl_allocation_params_get_type())
GST_GL_API
GType gst_gl_allocation_params_get_type (void);

/**
 * GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_ALLOC:
 *
 * GL Allocation flag indicating that the implementation should allocate the
 * necessary resources.
 */
#define GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_ALLOC (1 << 0)

/**
 * GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_WRAP_SYSMEM:
 *
 * GL Allocation flag for using the provided system memory data as storage.
 */
#define GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_WRAP_SYSMEM (1 << 1)

/**
 * GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_WRAP_GPU_HANDLE:
 *
 * GL Allocation flag for using the provided GPU handle as storage.
 */
#define GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_WRAP_GPU_HANDLE (1 << 2)

/**
 * GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_USER:
 *
 * Values >= than #GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_USER can be used for
 * user-defined purposes.
 */
#define GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_USER (1 << 16)

/**
 * GstGLAllocationParams:
 * @struct_size: the size of the struct (including and subclass data)
 * @copy: a #GstGLAllocationParamsCopyFunc
 * @free: a #GstGLAllocationParamsFreeFunc
 * @alloc_flags: allocation flags
 * @alloc_size: the allocation size
 * @alloc_params: the #GstAllocationParams
 * @context: a #GstGLContext
 * @notify: a #GDestroyNotify
 * @user_data: argument to call @notify with
 * @wrapped_data: the wrapped data pointer
 * @gl_handle: the wrapped OpenGL handle
 */
/* Because GstAllocationParams is not subclassable, start our own subclass
 * chain.  FIXME: 2.0 make GstAllocationParams subclassable */
struct _GstGLAllocationParams
{
  gsize                             struct_size;
  GstGLAllocationParamsCopyFunc     copy;
  GstGLAllocationParamsFreeFunc     free;

  guint                             alloc_flags;
  gsize                             alloc_size;
  GstAllocationParams              *alloc_params;
  GstGLContext                     *context;
  GDestroyNotify                    notify;
  gpointer                          user_data;

  /* GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_WRAP_SYSMEM only */
  gpointer                          wrapped_data;
  /* GST_GL_ALLOCATION_PARAMS_ALLOC_FLAG_WRAP_GPU_HANDLE only */
  gpointer                          gl_handle;

  /*< private >*/
  gpointer                          _padding[GST_PADDING];
};

GST_GL_API
gboolean                gst_gl_allocation_params_init       (GstGLAllocationParams * params,
                                                             gsize struct_size,
                                                             guint alloc_flags,
                                                             GstGLAllocationParamsCopyFunc copy,
                                                             GstGLAllocationParamsFreeFunc free,
                                                             GstGLContext * context,
                                                             gsize alloc_size,
                                                             const GstAllocationParams * alloc_params,
                                                             gpointer wrapped_data,
                                                             gpointer gl_handle,
                                                             gpointer user_data,
                                                             GDestroyNotify notify);

/* free with gst_gl_allocation_params_free */
GST_GL_API
GstGLAllocationParams * gst_gl_allocation_params_copy       (GstGLAllocationParams * src);

GST_GL_API
void                    gst_gl_allocation_params_free       (GstGLAllocationParams * params);

/* subclass usage */
GST_GL_API
void                    gst_gl_allocation_params_free_data  (GstGLAllocationParams * params);

/* subclass usage */
GST_GL_API
void                    gst_gl_allocation_params_copy_data  (GstGLAllocationParams * src,
                                                             GstGLAllocationParams * dest);

/**
 * GstGLBaseMemoryAllocatorAllocFunction:
 * @allocator: a #GstGLBaseMemoryAllocator
 * @params: the #GstGLAllocationParams to allocate the memory with
 *
 * Note: not called with a GL context current
 *
 * Returns: (transfer full) (nullable): a newly allocated #GstGLBaseMemory from @allocator and @params
 *
 * Since: 1.8
 */
typedef GstGLBaseMemory *   (*GstGLBaseMemoryAllocatorAllocFunction)        (GstGLBaseMemoryAllocator * allocator,
                                                                             GstGLAllocationParams * params);

/**
 * GstGLBaseMemoryAllocatorCreateFunction:
 * @mem: a #GstGLBaseMemory
 * @error: a #GError to use on failure
 *
 * As this virtual method is called with an OpenGL context current, use this
 * function to allocate and OpenGL resources needed for your application
 *
 * Returns: whether the creation succeeded
 *
 * Since: 1.8
 */
typedef gboolean            (*GstGLBaseMemoryAllocatorCreateFunction)       (GstGLBaseMemory * mem,
                                                                             GError ** error);

/**
 * GstGLBaseMemoryAllocatorMapFunction:
 * @mem: a #GstGLBaseMemory
 * @info: a #GstMapInfo to map with
 * @maxsize: the size to map
 *
 * Also see gst_memory_map();
 *
 * Returns: the mapped pointer
 *
 * Since: 1.8
 */
typedef gpointer            (*GstGLBaseMemoryAllocatorMapFunction)          (GstGLBaseMemory * mem,
                                                                             GstMapInfo * info,
                                                                             gsize maxsize);
/**
 * GstGLBaseMemoryAllocatorUnmapFunction:
 * @mem: a #GstGLBaseMemory
 * @info: a #GstMapInfo to map with
 *
 * Also see gst_memory_unmap();
 *
 * Since: 1.8
 */
typedef void                (*GstGLBaseMemoryAllocatorUnmapFunction)        (GstGLBaseMemory * mem,
                                                                             GstMapInfo * info);

/**
 * GstGLBaseMemoryAllocatorCopyFunction:
 * @mem: a #GstGLBaseMemory
 * @offset: the offset to copy from
 * @size: the number of bytes to copy
 *
 * Also see gst_memory_copy();
 *
 * Returns: (transfer full) (nullable): the newly copied #GstGLMemory or %NULL
 *
 * Since: 1.8
 */
typedef GstGLBaseMemory *   (*GstGLBaseMemoryAllocatorCopyFunction)         (GstGLBaseMemory * mem,
                                                                             gssize offset,
                                                                             gssize size);

/**
 * GstGLBaseMemoryAllocatorDestroyFunction:
 * @mem: a #GstGLBaseMemory
 *
 * Destroy any resources allocated throughout the lifetime of @mem
 *
 * Since: 1.8
 */
typedef void                (*GstGLBaseMemoryAllocatorDestroyFunction)      (GstGLBaseMemory * mem);

/**
 * GstGLBaseMemoryAllocator:
 *
 * Opaque #GstGLBaseMemoryAllocator struct
 *
 * Since: 1.8
 */
struct _GstGLBaseMemoryAllocator
{
  /*< private >*/
  GstAllocator parent;
  GstMemoryCopyFunction fallback_mem_copy;

  gpointer _padding[GST_PADDING];
};

/**
 * GstGLBaseMemoryAllocatorClass:
 * @parent_class: the parent class
 * @alloc: a #GstGLBaseMemoryAllocatorAllocFunction
 * @create: a #GstGLBaseMemoryAllocatorCreateFunction
 * @map: a #GstGLBaseMemoryAllocatorMapFunction
 * @unmap: a #GstGLBaseMemoryAllocatorUnmapFunction
 * @copy: a #GstGLBaseMemoryAllocatorCopyFunction
 * @destroy: a #GstGLBaseMemoryAllocatorDestroyFunction
 *
 * Since: 1.8
 */
struct _GstGLBaseMemoryAllocatorClass
{
  GstAllocatorClass parent_class;

  GstGLBaseMemoryAllocatorAllocFunction         alloc;

  GstGLBaseMemoryAllocatorCreateFunction        create;
  GstGLBaseMemoryAllocatorMapFunction           map;
  GstGLBaseMemoryAllocatorUnmapFunction         unmap;
  GstGLBaseMemoryAllocatorCopyFunction          copy;
  GstGLBaseMemoryAllocatorDestroyFunction       destroy;

  /*< private >*/
  gpointer                                      _padding[GST_PADDING];
};

#include <gst/gl/gstglconfig.h>
#include <gst/gl/gstglformat.h>

/**
 * GST_GL_BASE_MEMORY_ALLOCATOR_NAME:
 *
 * The name of the GL buffer allocator
 *
 * Since: 1.8
 */
#define GST_GL_BASE_MEMORY_ALLOCATOR_NAME   "GLBaseMemory"

GST_GL_API
void          gst_gl_base_memory_init_once (void);

GST_GL_API
gboolean      gst_is_gl_base_memory        (GstMemory * mem);

GST_GL_API
void          gst_gl_base_memory_init      (GstGLBaseMemory * mem,
                                            GstAllocator * allocator,
                                            GstMemory * parent,
                                            GstGLContext * context,
                                            const GstAllocationParams * params,
                                            gsize size,
                                            gpointer user_data,
                                            GDestroyNotify notify);

GST_GL_API
gboolean      gst_gl_base_memory_alloc_data (GstGLBaseMemory * gl_mem);

GST_GL_API
gboolean      gst_gl_base_memory_memcpy     (GstGLBaseMemory * src,
                                             GstGLBaseMemory * dest,
                                             gssize offset,
                                             gssize size);

GST_GL_API
GstGLBaseMemory *   gst_gl_base_memory_alloc    (GstGLBaseMemoryAllocator * allocator,
                                                 GstGLAllocationParams * params);

G_END_DECLS

#endif /* _GST_GL_BUFFER_H_ */
