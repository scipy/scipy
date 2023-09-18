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

GST_GL_EXT_BEGIN (blending,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  1, 2,
                  2, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, BlendEquation,
                     (GLenum                mode))
GST_GL_EXT_FUNCTION (void, BlendColor,
                     (GLclampf              red,
                      GLclampf              green,
                      GLclampf              blue,
                      GLclampf              alpha))
GST_GL_EXT_END ()

/* Optional, declared in 1.4 or GLES 1.2 */
GST_GL_EXT_BEGIN (blend_func_separate,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  1, 4,
                  2, 0,
                  "EXT\0",
                  "blend_func_separate\0")
GST_GL_EXT_FUNCTION (void, BlendFuncSeparate,
                     (GLenum                srcRGB,
                      GLenum                dstRGB,
                      GLenum                srcAlpha,
                      GLenum                dstAlpha))
GST_GL_EXT_END ()

/* Optional, declared in 2.0 */
GST_GL_EXT_BEGIN (blend_equation_separate,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 0,
                  2, 0,
                  "EXT\0",
                  "blend_equation_separate\0")
GST_GL_EXT_FUNCTION (void, BlendEquationSeparate,
                     (GLenum                modeRGB,
                      GLenum                modeAlpha))
GST_GL_EXT_END ()

/* GL and GLES 2.0 apis */
GST_GL_EXT_BEGIN (two_point_zero_api,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 0,
                  2, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, StencilFuncSeparate,
                     (GLenum face, GLenum func, GLint ref, GLuint mask))
GST_GL_EXT_FUNCTION (void, StencilMaskSeparate,
                     (GLenum face, GLuint mask))
GST_GL_EXT_FUNCTION (void, StencilOpSeparate,
                     (GLenum face, GLenum fail, GLenum zfail, GLenum zpass))
GST_GL_EXT_END ()
