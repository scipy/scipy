/*
 * GStreamer
 * Copyright (C) 2014 Jan Schmidt <jan@centricular.com>
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

#ifndef _GST_GL_VIEW_CONVERT_H_
#define _GST_GL_VIEW_CONVERT_H_

#include <gst/gstmemory.h>
#include <gst/video/video.h>

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS
#define GST_TYPE_GL_VIEW_CONVERT            (gst_gl_view_convert_get_type())
#define GST_GL_VIEW_CONVERT(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_GL_VIEW_CONVERT,GstGLViewConvert))
#define GST_IS_GL_VIEW_CONVERT(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_GL_VIEW_CONVERT))
#define GST_GL_VIEW_CONVERT_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass) ,GST_TYPE_GL_VIEW_CONVERT,GstGLViewConvertClass))
#define GST_IS_GL_VIEW_CONVERT_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass) ,GST_TYPE_GL_VIEW_CONVERT))
#define GST_GL_VIEW_CONVERT_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj) ,GST_TYPE_GL_VIEW_CONVERT,GstGLViewConvertClass))

/**
 * GstGLStereoDownmix:
 * @GST_GL_STEREO_DOWNMIX_ANAGLYPH_GREEN_MAGENTA_DUBOIS: Dubois optimised Green-Magenta anaglyph
 * @GST_GL_STEREO_DOWNMIX_ANAGLYPH_RED_CYAN_DUBOIS: Dubois optimised Red-Cyan anaglyph
 * @GST_GL_STEREO_DOWNMIX_ANAGLYPH_AMBER_BLUE_DUBOIS: Dubois optimised Amber-Blue anaglyph
 *
 * Output anaglyph type to generate when downmixing to mono
 */
typedef enum {
  GST_GL_STEREO_DOWNMIX_ANAGLYPH_GREEN_MAGENTA_DUBOIS,
  GST_GL_STEREO_DOWNMIX_ANAGLYPH_RED_CYAN_DUBOIS,
  GST_GL_STEREO_DOWNMIX_ANAGLYPH_AMBER_BLUE_DUBOIS,
} GstGLStereoDownmix;

#ifndef GST_DISABLE_DEPRECATED
#define GST_TYPE_GL_STEREO_DOWNMIX_MODE_TYPE GST_TYPE_GL_STEREO_DOWNMIX

GST_GL_API
GType gst_gl_stereo_downmix_mode_get_type (void);
#endif

/**
 * GstGLViewConvert:
 *
 * #GstGLViewConvert is an opaque struct and should only be accessed through the
 * provided api.
 */
struct _GstGLViewConvert
{
  GstObject object;

  GstGLContext *context;

  GstGLShader *shader;

  GstVideoMultiviewMode input_mode_override;
  GstVideoMultiviewFlags input_flags_override;
  GstVideoMultiviewMode output_mode_override;
  GstVideoMultiviewFlags output_flags_override;

  GstGLStereoDownmix downmix_mode;

  GstVideoInfo in_info;
  GstVideoInfo out_info;

  GstGLTextureTarget from_texture_target;
  GstGLTextureTarget to_texture_target;
  gboolean caps_passthrough;

  gboolean initted;
  gboolean reconfigure;

  GstGLFramebuffer *fbo;

  /*< private >*/
  GstGLViewConvertPrivate *priv;

  gpointer _padding[GST_PADDING];
};

/**
 * GstGLViewConvertClass:
 *
 * Opaque #GstGLViewConvertClass struct
 */
struct _GstGLViewConvertClass
{
  /*< private >*/
  GstObjectClass object_class;

  gpointer                  _padding[GST_PADDING];
};

GST_GL_API
GType gst_gl_view_convert_get_type (void);
GST_GL_API
GstGLViewConvert * gst_gl_view_convert_new (void);

GST_GL_API
gboolean  gst_gl_view_convert_set_caps (GstGLViewConvert * viewconvert, GstCaps * in_caps, GstCaps * out_caps);
GST_GL_API
GstCaps * gst_gl_view_convert_transform_caps (GstGLViewConvert * viewconvert,
    GstPadDirection direction, GstCaps * caps, GstCaps * filter);
GST_GL_API
GstCaps * gst_gl_view_convert_fixate_caps (GstGLViewConvert *viewconvert,
    GstPadDirection direction, GstCaps * caps, GstCaps * othercaps);
GST_GL_API
GstFlowReturn gst_gl_view_convert_submit_input_buffer (GstGLViewConvert *viewconvert,
    gboolean is_discont, GstBuffer * input);
GST_GL_API
GstFlowReturn gst_gl_view_convert_get_output (GstGLViewConvert *viewconvert,
    GstBuffer ** outbuf_ptr);

GST_GL_API
GstBuffer * gst_gl_view_convert_perform (GstGLViewConvert * viewconvert, GstBuffer *inbuf);
GST_GL_API
void gst_gl_view_convert_reset (GstGLViewConvert * viewconvert);
GST_GL_API
void gst_gl_view_convert_set_context (GstGLViewConvert *viewconvert, GstGLContext * context);

G_END_DECLS
#endif /* _GST_GL_VIEW_CONVERT_H_ */
