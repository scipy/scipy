/*
 * GStreamer
 * Copyright (C) 2013 Julien Isorce <julien.isorce@gmail.com>
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

#ifndef __GST_GL_FWD_H__
#define __GST_GL_FWD_H__

#include <gst/gst.h>

#include <gst/gl/gstglapi.h>

G_BEGIN_DECLS

typedef struct _GstGLDisplay GstGLDisplay;
typedef struct _GstGLDisplayClass GstGLDisplayClass;
typedef struct _GstGLDisplayPrivate GstGLDisplayPrivate;

typedef struct _GstGLContext GstGLContext;
typedef struct _GstGLContextClass GstGLContextClass;
typedef struct _GstGLContextPrivate GstGLContextPrivate;

typedef struct _GstGLWindow        GstGLWindow;
typedef struct _GstGLWindowPrivate GstGLWindowPrivate;
typedef struct _GstGLWindowClass   GstGLWindowClass;

typedef struct _GstGLBaseMemory GstGLBaseMemory;
typedef struct _GstGLBaseMemoryAllocator GstGLBaseMemoryAllocator;
typedef struct _GstGLBaseMemoryAllocatorClass GstGLBaseMemoryAllocatorClass;

typedef struct _GstGLBuffer GstGLBuffer;
typedef struct _GstGLBufferAllocator GstGLBufferAllocator;
typedef struct _GstGLBufferAllocatorClass GstGLBufferAllocatorClass;

typedef struct _GstGLMemory GstGLMemory;
typedef struct _GstGLMemoryAllocator GstGLMemoryAllocator;
typedef struct _GstGLMemoryAllocatorClass GstGLMemoryAllocatorClass;

typedef struct _GstGLMemoryPBO GstGLMemoryPBO;
typedef struct _GstGLMemoryPBOAllocator GstGLMemoryPBOAllocator;
typedef struct _GstGLMemoryPBOAllocatorClass GstGLMemoryPBOAllocatorClass;

typedef struct _GstGLRenderbuffer GstGLRenderbuffer;
typedef struct _GstGLRenderbufferAllocator GstGLRenderbufferAllocator;
typedef struct _GstGLRenderbufferAllocatorClass GstGLRenderbufferAllocatorClass;

typedef struct _GstGLFramebuffer GstGLFramebuffer;
typedef struct _GstGLFramebufferClass GstGLFramebufferClass;

typedef struct _GstGLSLStage        GstGLSLStage;
typedef struct _GstGLSLStagePrivate GstGLSLStagePrivate;
typedef struct _GstGLSLStageClass   GstGLSLStageClass;

typedef struct _GstGLShader        GstGLShader;
typedef struct _GstGLShaderPrivate GstGLShaderPrivate;
typedef struct _GstGLShaderClass   GstGLShaderClass;

typedef struct _GstGLUpload GstGLUpload;
typedef struct _GstGLUploadClass GstGLUploadClass;
typedef struct _GstGLUploadPrivate GstGLUploadPrivate;

typedef struct _GstGLBufferPool GstGLBufferPool;
typedef struct _GstGLBufferPoolClass GstGLBufferPoolClass;
typedef struct _GstGLBufferPoolPrivate GstGLBufferPoolPrivate;

typedef struct _GstGLColorConvert GstGLColorConvert;
typedef struct _GstGLColorConvertClass GstGLColorConvertClass;
typedef struct _GstGLColorConvertPrivate GstGLColorConvertPrivate;

typedef struct _GstGLBaseFilter GstGLBaseFilter;
typedef struct _GstGLBaseFilterClass GstGLBaseFilterClass;
typedef struct _GstGLBaseFilterPrivate GstGLBaseFilterPrivate;

typedef struct _GstGLBaseSrc GstGLBaseSrc;
typedef struct _GstGLBaseSrcClass GstGLBaseSrcClass;
typedef struct _GstGLBaseSrcPrivate GstGLBaseSrcPrivate;

typedef struct _GstGLFilter GstGLFilter;
typedef struct _GstGLFilterClass GstGLFilterClass;

typedef struct _GstGLViewConvert GstGLViewConvert;
typedef struct _GstGLViewConvertClass GstGLViewConvertClass;
typedef struct _GstGLViewConvertPrivate GstGLViewConvertPrivate;

typedef struct _GstGLOverlayCompositor GstGLOverlayCompositor;
typedef struct _GstGLOverlayCompositorClass GstGLOverlayCompositorClass;

typedef struct _GstGLQuery GstGLQuery;

typedef struct _GstGLFuncs GstGLFuncs;

typedef struct _GstGLAsyncDebug GstGLAsyncDebug;

#include <gst/gl/gstgl_enums.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLBaseFilter, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLBaseMemoryAllocator, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLBaseSrc, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLBufferAllocator, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLBufferPool, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLColorConvert, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLContext, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLDisplay, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLFilter, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLMemoryAllocator, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLMemoryPBOAllocator, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLOverlayCompositor, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLSLStage, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLShader, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLUpload, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLViewConvert, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstGLWindow, gst_object_unref)

G_END_DECLS

#endif /* __GST_GL_FWD_H__ */
