/*
 * GStreamer
 * Copyright (C) 2012 Matthew Waters <ystree00@gmail.com>
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

#ifndef __GST_GL_COLOR_CONVERT_H__
#define __GST_GL_COLOR_CONVERT_H__

#include <gst/video/video.h>
#include <gst/gstmemory.h>

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_color_convert_get_type (void);
#define GST_TYPE_GL_COLOR_CONVERT (gst_gl_color_convert_get_type())
#define GST_GL_COLOR_CONVERT(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_GL_COLOR_CONVERT,GstGLColorConvert))
#define GST_GL_COLOR_CONVERT_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_GL_DISPLAY,GstGLColorConvertClass))
#define GST_IS_GL_COLOR_CONVERT(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_GL_COLOR_CONVERT))
#define GST_IS_GL_COLOR_CONVERT_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_GL_COLOR_CONVERT))
#define GST_GL_COLOR_CONVERT_CAST(obj) ((GstGLColorConvert*)(obj))

/**
 * GstGLColorConvert
 *
 * Opaque #GstGLColorConvert object
 */
struct _GstGLColorConvert
{
  /*< private >*/
  GstObject        parent;

  GstGLContext    *context;

  /* input data */
  GstVideoInfo     in_info;
  GstVideoInfo     out_info;

  gboolean         initted;
  gboolean         passthrough;

  GstBuffer *    inbuf;
  GstBuffer *    outbuf;

  /* used for the conversion */
  GstGLFramebuffer *fbo;
  GstGLShader     *shader;

  /*< private >*/
  GstGLColorConvertPrivate *priv;

  gpointer _reserved[GST_PADDING];
};

/**
 * GstGLColorConvertClass:
 *
 * The #GstGLColorConvertClass struct only contains private data
 */
struct _GstGLColorConvertClass
{
  /*< private >*/
  GstObjectClass object_class;

  gpointer _padding[GST_PADDING];
};

/**
 * GST_GL_COLOR_CONVERT_EXT_FORMATS: (skip)
 *
 */
#if G_BYTE_ORDER == G_LITTLE_ENDIAN
#define GST_GL_COLOR_CONVERT_EXT_FORMATS \
    ", BGR10A2_LE, RGB10A2_LE, P010_10LE, P012_LE, P016_LE, Y212_LE, Y412_LE"
#else
#define GST_GL_COLOR_CONVERT_EXT_FORMATS \
    ", P010_10BE, P012_BE, P016_BE, Y212_BE, Y412_BE"
#endif

/**
 * GST_GL_COLOR_CONVERT_FORMATS:
 *
 * The currently supported formats that can be converted
 */
#define GST_GL_COLOR_CONVERT_FORMATS "{ RGBA, RGB, RGBx, BGR, BGRx, BGRA, xRGB, " \
                               "xBGR, ARGB, ABGR, GBRA, GBR, RGBP, BGRP, Y444, I420, YV12, Y42B, " \
                               "Y41B, NV12, NV21, NV16, NV61, YUY2, UYVY, Y210, AYUV, " \
                               "VUYA, Y410, GRAY8, GRAY16_LE, GRAY16_BE, " \
                               "RGB16, BGR16, ARGB64, A420, AV12, NV12_16L32S, NV12_4L4" \
                               GST_GL_COLOR_CONVERT_EXT_FORMATS "}"

/**
 * GST_GL_COLOR_CONVERT_VIDEO_CAPS:
 *
 * The currently supported #GstCaps that can be converted
 */
#define GST_GL_COLOR_CONVERT_VIDEO_CAPS \
    "video/x-raw(" GST_CAPS_FEATURE_MEMORY_GL_MEMORY "), "              \
    "format = (string) " GST_GL_COLOR_CONVERT_FORMATS ", "              \
    "width = " GST_VIDEO_SIZE_RANGE ", "                                \
    "height = " GST_VIDEO_SIZE_RANGE ", "                               \
    "framerate = " GST_VIDEO_FPS_RANGE ", "                             \
    "texture-target = (string) { 2D, rectangle, external-oes } "        \
    " ; "                                                               \
    "video/x-raw(" GST_CAPS_FEATURE_MEMORY_GL_MEMORY ","                \
    GST_CAPS_FEATURE_META_GST_VIDEO_OVERLAY_COMPOSITION "), "           \
    "format = (string) " GST_GL_COLOR_CONVERT_FORMATS ", "              \
    "width = " GST_VIDEO_SIZE_RANGE ", "                                \
    "height = " GST_VIDEO_SIZE_RANGE ", "                               \
    "framerate = " GST_VIDEO_FPS_RANGE ", "                             \
    "texture-target = (string) { 2D, rectangle, external-oes }"

GST_GL_API
GstGLColorConvert * gst_gl_color_convert_new (GstGLContext * context);

GST_GL_API
GstCaps *   gst_gl_color_convert_transform_caps (GstGLContext * context,
                                                 GstPadDirection direction,
                                                 GstCaps * caps,
                                                 GstCaps * filter);
GST_GL_API
GstCaps *   gst_gl_color_convert_fixate_caps    (GstGLContext * context,
                                                 GstPadDirection direction,
                                                 GstCaps * caps,
                                                 GstCaps * other);
GST_GL_API
gboolean    gst_gl_color_convert_set_caps    (GstGLColorConvert * convert,
                                              GstCaps           * in_caps,
                                              GstCaps           * out_caps);
GST_GL_API
gboolean    gst_gl_color_convert_decide_allocation (GstGLColorConvert   * convert,
                                                    GstQuery            * query);

GST_GL_API
GstBuffer * gst_gl_color_convert_perform    (GstGLColorConvert * convert, GstBuffer * inbuf);

G_END_DECLS

#endif /* __GST_GL_COLOR_CONVERT_H__ */
