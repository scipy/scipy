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

GST_GL_EXT_BEGIN (multitexture_part1,
                  GST_GL_API_OPENGL |
                  GST_GL_API_GLES1,
                  1, 3,
                  1, 0,
                  "ARB\0",
                  "multitexture\0")
GST_GL_EXT_FUNCTION (void, ClientActiveTexture,
                     (GLenum                texture))
GST_GL_EXT_END ()

/* These are the core GL functions which are available when the API
   supports fixed-function (ie, GL and GLES1.1) */
GST_GL_EXT_BEGIN (fixed_function_core,
                  GST_GL_API_OPENGL |
                  GST_GL_API_GLES1,
                  0, 0,
                  1, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, AlphaFunc,
                     (GLenum func, GLclampf ref))
GST_GL_EXT_FUNCTION (void, Fogf,
                     (GLenum pname, GLfloat param))
GST_GL_EXT_FUNCTION (void, Fogfv,
                     (GLenum pname, const GLfloat *params))
GST_GL_EXT_FUNCTION (void, LoadMatrixf,
                     (const GLfloat *m))
GST_GL_EXT_FUNCTION (void, Materialfv,
                     (GLenum face, GLenum pname, const GLfloat *params))
GST_GL_EXT_FUNCTION (void, PointSize,
                     (GLfloat size))
GST_GL_EXT_FUNCTION (void, TexEnvfv,
                     (GLenum target, GLenum pname, const GLfloat *params))
GST_GL_EXT_FUNCTION (void, Color4ub,
                     (GLubyte red, GLubyte green, GLubyte blue, GLubyte alpha))
GST_GL_EXT_FUNCTION (void, ColorPointer,
                     (GLint size,
                      GLenum type,
                      GLsizei stride,
                      const GLvoid *pointer))
GST_GL_EXT_FUNCTION (void, DisableClientState,
                     (GLenum array))
GST_GL_EXT_FUNCTION (void, EnableClientState,
                     (GLenum array))
GST_GL_EXT_FUNCTION (void, LoadIdentity,
                     (void))
GST_GL_EXT_FUNCTION (void, MatrixMode,
                     (GLenum mode))
GST_GL_EXT_FUNCTION (void, NormalPointer,
                     (GLenum type, GLsizei stride, const GLvoid *pointer))
GST_GL_EXT_FUNCTION (void, TexCoordPointer,
                     (GLint size,
                      GLenum type,
                      GLsizei stride,
                      const GLvoid *pointer))
GST_GL_EXT_FUNCTION (void, TexEnvi,
                     (GLenum target,
                      GLenum pname,
                      GLint param))
GST_GL_EXT_FUNCTION (void, VertexPointer,
                     (GLint size,
                      GLenum type,
                      GLsizei stride,
                      const GLvoid *pointer))
GST_GL_EXT_FUNCTION (void, PushMatrix,
                     (void))
GST_GL_EXT_FUNCTION (void, PopMatrix,
                     (void))
GST_GL_EXT_END ()

/* Eventually we want to remove this category */
GST_GL_EXT_BEGIN (fixed_function_gl_only,
                  GST_GL_API_OPENGL,
                  0, 0,
                  0, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, PushAttrib,
                     (GLbitfield            mask))
GST_GL_EXT_FUNCTION (void, PopAttrib,
                     (void))
GST_GL_EXT_FUNCTION (void, TexImage1D,
                     (GLenum                target,
                      GLint                 level,
                      GLint                 internalFormat,
                      GLsizei               width,
                      GLint                 border,
                      GLenum                format,
                      GLenum                type,
                      const GLvoid         *data))
GST_GL_EXT_FUNCTION (void, Rotatef,
                     (GLfloat angle, GLfloat x, GLfloat y, GLfloat z))
GST_GL_EXT_FUNCTION (void, Translatef,
                     (GLfloat x, GLfloat y, GLfloat z))
GST_GL_EXT_FUNCTION (void, Scalef,
                     (GLfloat x, GLfloat y, GLfloat z))
GST_GL_EXT_FUNCTION (void, Lightfv,
                     (GLenum light, GLenum pname, const GLfloat *params))
GST_GL_EXT_FUNCTION (void, ColorMaterial,
                     (GLenum face, GLenum pname))
GST_GL_EXT_FUNCTION (void, ShadeModel,
                     (GLenum value))
GST_GL_EXT_END ()
