/*
 * GStreamer
 * Copyright (C) 2008 Filippo Argiolas <filippo.argiolas@gmail.com>
 * Copyright (C) 2014 Julien Isorce <julien.isorce@collabora.co.uk>
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

#ifndef __GST_GL_SHADER_H__
#define __GST_GL_SHADER_H__

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_shader_get_type (void);
#define GST_TYPE_GL_SHADER         (gst_gl_shader_get_type())

#define GST_GL_SHADER(o)           (G_TYPE_CHECK_INSTANCE_CAST((o), GST_TYPE_GL_SHADER, GstGLShader))
#define GST_GL_SHADER_CLASS(k)     (G_TYPE_CHECK_CLASS((k), GST_TYPE_GL_SHADER, GstGLShaderClass))
#define GST_IS_GL_SHADER(o)        (G_TYPE_CHECK_INSTANCE_TYPE((o), GST_TYPE_GL_SHADER))
#define GST_IS_GL_SHADER_CLASS(k)  (G_TYPE_CHECK_CLASS_TYPE((k), GST_TYPE_GL_SHADER))
#define GST_GL_SHADER_GET_CLASS(o) (G_TYPE_INSTANCE_GET_CLASS((o), GST_TYPE_GL_SHADER, GstGLShaderClass))

struct _GstGLShader
{
  GstObject parent;

  GstGLContext *context;

  /*< private >*/
  GstGLShaderPrivate *priv;

  gpointer _padding[GST_PADDING];
};

struct _GstGLShaderClass {
  /*< private >*/
  GstObjectClass parent_class;
};

GST_GL_API
GstGLShader * gst_gl_shader_new                     (GstGLContext *context);
GST_GL_API
GstGLShader * gst_gl_shader_new_with_stages         (GstGLContext * context, GError ** error, ...);
GST_GL_API
GstGLShader * gst_gl_shader_new_link_with_stages    (GstGLContext * context, GError ** error, ...);
GST_GL_API
GstGLShader * gst_gl_shader_new_default             (GstGLContext * context, GError ** error);

GST_GL_API
gboolean gst_gl_shader_attach                       (GstGLShader * shader, GstGLSLStage * stage);
GST_GL_API
gboolean gst_gl_shader_attach_unlocked              (GstGLShader * shader, GstGLSLStage * stage);

GST_GL_API
void     gst_gl_shader_detach                       (GstGLShader * shader, GstGLSLStage * stage);
GST_GL_API
void     gst_gl_shader_detach_unlocked              (GstGLShader * shader, GstGLSLStage * stage);

GST_GL_API
gboolean gst_gl_shader_compile_attach_stage         (GstGLShader * shader,
                                                     GstGLSLStage *stage,
                                                     GError ** error);
GST_GL_API
gboolean gst_gl_shader_link                         (GstGLShader * shader, GError ** error);
GST_GL_API
gboolean gst_gl_shader_is_linked                    (GstGLShader *shader);

GST_GL_API
int gst_gl_shader_get_program_handle                (GstGLShader * shader);

GST_GL_API
void gst_gl_shader_release                          (GstGLShader *shader);
GST_GL_API
void gst_gl_shader_release_unlocked                 (GstGLShader * shader);
GST_GL_API
void gst_gl_shader_use                              (GstGLShader *shader);
GST_GL_API
void gst_gl_context_clear_shader                    (GstGLContext *context);

GST_GL_API
void gst_gl_shader_set_uniform_1i           (GstGLShader *shader, const gchar *name, gint value);
GST_GL_API
void gst_gl_shader_set_uniform_1iv          (GstGLShader *shader, const gchar *name, guint count, const gint *value);
GST_GL_API
void gst_gl_shader_set_uniform_1f           (GstGLShader *shader, const gchar *name, gfloat value);
GST_GL_API
void gst_gl_shader_set_uniform_1fv          (GstGLShader *shader, const gchar *name, guint count, const gfloat *value);
GST_GL_API
void gst_gl_shader_set_uniform_2i           (GstGLShader *shader, const gchar *name, gint v0,     gint v1);
GST_GL_API
void gst_gl_shader_set_uniform_2iv          (GstGLShader *shader, const gchar *name, guint count, const gint *value);
GST_GL_API
void gst_gl_shader_set_uniform_2f           (GstGLShader *shader, const gchar *name, gfloat v0,   gfloat v1);
GST_GL_API
void gst_gl_shader_set_uniform_2fv          (GstGLShader *shader, const gchar *name, guint count, const gfloat *value);
GST_GL_API
void gst_gl_shader_set_uniform_3i           (GstGLShader *shader, const gchar *name, gint v0,     gint v1,       gint v2);
GST_GL_API
void gst_gl_shader_set_uniform_3iv          (GstGLShader *shader, const gchar *name, guint count, const gint * value);
GST_GL_API
void gst_gl_shader_set_uniform_3f           (GstGLShader *shader, const gchar *name, gfloat v0,   gfloat v1,     gfloat v2);
GST_GL_API
void gst_gl_shader_set_uniform_3fv          (GstGLShader *shader, const gchar *name, guint count, const gfloat *value);
GST_GL_API
void gst_gl_shader_set_uniform_4i           (GstGLShader *shader, const gchar *name, gint v0,     gint v1,       gint v2,   gint v3);
GST_GL_API
void gst_gl_shader_set_uniform_4iv          (GstGLShader *shader, const gchar *name, guint count, const gint *value);
GST_GL_API
void gst_gl_shader_set_uniform_4f           (GstGLShader *shader, const gchar *name, gfloat v0,   gfloat v1,     gfloat v2, gfloat v3);
GST_GL_API
void gst_gl_shader_set_uniform_4fv          (GstGLShader *shader, const gchar *name, guint count, const gfloat *value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_2fv   (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_3fv   (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_4fv   (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_2x3fv (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_2x4fv (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_3x2fv (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_3x4fv (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_4x2fv (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);
GST_GL_API
void gst_gl_shader_set_uniform_matrix_4x3fv (GstGLShader *shader, const gchar *name, gint count, gboolean transpose, const gfloat* value);

GST_GL_API
gint gst_gl_shader_get_attribute_location  (GstGLShader *shader, const gchar *name);
GST_GL_API
void gst_gl_shader_bind_attribute_location (GstGLShader * shader, guint index, const gchar * name);
GST_GL_API
void gst_gl_shader_bind_frag_data_location (GstGLShader * shader, guint index, const gchar * name);

G_END_DECLS

#endif /* __GST_GL_SHADER_H__ */
