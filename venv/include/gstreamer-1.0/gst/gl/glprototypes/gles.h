/*
 * GStreamer
 * Copyright (C) 2012 Matthew Waters <ystreet00@gmail.com>
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
/*
 * Cogl
 *
 * An object oriented GL/GLES Abstraction/Utility Layer
 *
 * Copyright (C) 2009, 2011 Intel Corporation.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library. If not, see <http://www.gnu.org/licenses/>.
 */

GST_GL_EXT_BEGIN (only_in_both_gles,
                  GST_GL_API_GLES1 |
                  GST_GL_API_GLES2,
                  255, 255,
                  1, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, DepthRangef,
                     (GLfloat near_val, GLfloat far_val))
GST_GL_EXT_FUNCTION (void, ClearDepthf,
                     (GLclampf depth))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (only_in_gles1,
                GST_GL_API_GLES1,
                255, 255,
                1, 0,
                "\0",
                "\0")
GST_GL_EXT_FUNCTION (void, ClipPlanef, (GLenum plane, const GLfloat *equation))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (gles2_only_api,
                  GST_GL_API_GLES2,
                  255, 255,
                  2, 0,
                  "ARB:\0",
                  "ES2_compatibility\0")
GST_GL_EXT_FUNCTION (void, ReleaseShaderCompiler, (void))
GST_GL_EXT_FUNCTION (void, GetShaderPrecisionFormat,
                     (GLenum shadertype,
                      GLenum precisiontype,
                      GLint* range,
                      GLint* precision))
GST_GL_EXT_FUNCTION (void, ShaderBinary,
                     (GLsizei n,
                      const GLuint* shaders,
                      GLenum binaryformat,
                      const GLvoid* binary,
                      GLsizei length))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (IMG_multisampled_render_to_texture,
                  GST_GL_API_NONE,
                  255, 255,
                  255, 255, /* not in either GLES */
                  "\0",
                  "IMG_multisampled_render_to_texture\0")
GST_GL_EXT_FUNCTION (void, RenderbufferStorageMultisampleIMG,
                     (GLenum           target,
                      GLsizei          samples,
                      GLenum           internal_format,
                      GLsizei          width,
                      GLsizei          height))
GST_GL_EXT_FUNCTION (void, FramebufferTexture2DMultisampleIMG,
                     (GLenum           target,
                      GLenum           attachment,
                      GLenum           textarget,
                      GLuint           texture,
                      GLint            level,
                      GLsizei          samples))
GST_GL_EXT_END ()
