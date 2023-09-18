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

#ifndef __GST_GLSL_STAGE_H__
#define __GST_GLSL_STAGE_H__

#include <gst/gl/gstglsl.h>

G_BEGIN_DECLS

#define GST_TYPE_GLSL_STAGE         (gst_glsl_stage_get_type())
#define GST_GLSL_STAGE(o)           (G_TYPE_CHECK_INSTANCE_CAST((o), GST_TYPE_GLSL_STAGE, GstGLSLStage))
#define GST_GLSL_STAGE_CLASS(k)     (G_TYPE_CHECK_CLASS((k), GST_TYPE_GLSL_STAGE, GstGLSLStageClass))
#define GST_IS_GLSL_STAGE(o)        (G_TYPE_CHECK_INSTANCE_TYPE((o), GST_TYPE_GLSL_STAGE))
#define GST_IS_GLSL_STAGE_CLASS(k)  (G_TYPE_CHECK_CLASS_TYPE((k), GST_TYPE_GLSL_STAGE))
#define GST_GLSL_STAGE_GET_CLASS(o) (G_TYPE_INSTANCE_GET_CLASS((o), GST_TYPE_GLSL_STAGE, GstGLSLStageClass))

/**
 * GstGLSLStage:
 *
 * Opaque #GstGLSLStage struct
 */
struct _GstGLSLStage
{
  /*< private >*/
  GstObject parent;

  GstGLContext *context;

  GstGLSLStagePrivate *priv;

  gpointer _padding[GST_PADDING];
};

/**
 * GstGLSLStageClass:
 *
 * Opaque #GstGLSLStageClass struct
 */
struct _GstGLSLStageClass
{
  /*< private >*/
  GstObjectClass parent;

  gpointer _padding[GST_PADDING];
};

GST_GL_API
GType          gst_glsl_stage_get_type          (void);
GST_GL_API
GstGLSLStage * gst_glsl_stage_new               (GstGLContext * context, guint type);
GST_GL_API
GstGLSLStage * gst_glsl_stage_new_with_string   (GstGLContext * context,
                                                 guint type,
                                                 GstGLSLVersion version,
                                                 GstGLSLProfile profile,
                                                 const gchar * str);
GST_GL_API
GstGLSLStage * gst_glsl_stage_new_with_strings  (GstGLContext * context,
                                                 guint type,
                                                 GstGLSLVersion version,
                                                 GstGLSLProfile profile,
                                                 gint n_strings,
                                                 const gchar ** str);

GST_GL_API
GstGLSLStage * gst_glsl_stage_new_default_fragment (GstGLContext * context);
GST_GL_API
GstGLSLStage * gst_glsl_stage_new_default_vertex   (GstGLContext * context);

GST_GL_API
guint          gst_glsl_stage_get_handle        (GstGLSLStage * stage);
GST_GL_API
GstGLSLProfile gst_glsl_stage_get_profile       (GstGLSLStage * stage);
GST_GL_API
GstGLSLVersion gst_glsl_stage_get_version       (GstGLSLStage * stage);
GST_GL_API
guint          gst_glsl_stage_get_shader_type   (GstGLSLStage * stage);
GST_GL_API
gboolean       gst_glsl_stage_set_strings       (GstGLSLStage * stage,
                                                 GstGLSLVersion version,
                                                 GstGLSLProfile profile,
                                                 gint n_strings,
                                                 const gchar ** str);
GST_GL_API
gboolean       gst_glsl_stage_compile           (GstGLSLStage * stage,
                                                 GError ** error);

G_END_DECLS

#endif /* __GST_GLSL_STAGE_H__ */
