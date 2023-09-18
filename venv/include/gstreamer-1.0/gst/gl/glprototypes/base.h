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


/* These are the core GL functions which we assume will always be
   available */
GST_GL_EXT_BEGIN (core,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES1 | GST_GL_API_GLES2,
                  1, 0,
                  1, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, BindTexture,
                     (GLenum target, GLuint texture))
GST_GL_EXT_FUNCTION (void, BlendFunc,
                     (GLenum sfactor, GLenum dfactor))
GST_GL_EXT_FUNCTION (void, Clear,
                     (GLbitfield mask))
GST_GL_EXT_FUNCTION (void, ClearColor,
                     (GLclampf red,
                      GLclampf green,
                      GLclampf blue,
                      GLclampf alpha))
GST_GL_EXT_FUNCTION (void, ClearStencil,
                     (GLint s))
GST_GL_EXT_FUNCTION (void, ColorMask,
                     (GLboolean red,
                      GLboolean green,
                      GLboolean blue,
                      GLboolean alpha))
GST_GL_EXT_FUNCTION (void, CopyTexSubImage2D,
                     (GLenum target,
                      GLint level,
                      GLint xoffset,
                      GLint yoffset,
                      GLint x,
                      GLint y,
                      GLsizei width,
                      GLsizei height))
GST_GL_EXT_FUNCTION (void, DeleteTextures,
                     (GLsizei n, const GLuint* textures))
GST_GL_EXT_FUNCTION (void, DepthFunc,
                     (GLenum func))
GST_GL_EXT_FUNCTION (void, DepthMask,
                     (GLboolean flag))
GST_GL_EXT_FUNCTION (void, Disable,
                     (GLenum cap))
GST_GL_EXT_FUNCTION (void, DrawArrays,
                     (GLenum mode, GLint first, GLsizei count))
GST_GL_EXT_FUNCTION (void, DrawElements,
                     (GLenum mode,
                      GLsizei count,
                      GLenum type,
                      const GLvoid* indices))
GST_GL_EXT_FUNCTION (void, Enable,
                     (GLenum cap))
GST_GL_EXT_FUNCTION (void, Finish,
                     (void))
GST_GL_EXT_FUNCTION (void, Flush,
                     (void))
GST_GL_EXT_FUNCTION (void, FrontFace,
                     (GLenum mode))
GST_GL_EXT_FUNCTION (void, CullFace,
                     (GLenum mode))
GST_GL_EXT_FUNCTION (void, GenTextures,
                     (GLsizei n, GLuint* textures))
GST_GL_EXT_FUNCTION (GLenum, GetError,
                     (void))
GST_GL_EXT_FUNCTION (void, GetIntegerv,
                     (GLenum pname, GLint* params))
GST_GL_EXT_FUNCTION (void, GetBooleanv,
                     (GLenum pname, GLboolean* params))
GST_GL_EXT_FUNCTION (void, GetFloatv,
                     (GLenum pname, GLfloat* params))
GST_GL_EXT_FUNCTION (const GLubyte*, GetString,
                     (GLenum name))
GST_GL_EXT_FUNCTION (void, Hint,
                     (GLenum target, GLenum mode))
GST_GL_EXT_FUNCTION (GLboolean, IsTexture,
                     (GLuint texture))
GST_GL_EXT_FUNCTION (void, PixelStorei,
                     (GLenum pname, GLint param))
GST_GL_EXT_FUNCTION (void, ReadPixels,
                     (GLint x,
                      GLint y,
                      GLsizei width,
                      GLsizei height,
                      GLenum format,
                      GLenum type,
                      GLvoid* pixels))
GST_GL_EXT_FUNCTION (void, Scissor,
                     (GLint x, GLint y, GLsizei width, GLsizei height))
GST_GL_EXT_FUNCTION (void, StencilFunc,
                     (GLenum func, GLint ref, GLuint mask))
GST_GL_EXT_FUNCTION (void, StencilMask,
                     (GLuint mask))
GST_GL_EXT_FUNCTION (void, StencilOp,
                     (GLenum fail, GLenum zfail, GLenum zpass))
GST_GL_EXT_FUNCTION (void, TexImage2D,
                     (GLenum target,
                      GLint level,
                      GLint internalformat,
                      GLsizei width,
                      GLsizei height,
                      GLint border,
                      GLenum format,
                      GLenum type,
                      const GLvoid* pixels))
GST_GL_EXT_FUNCTION (void, TexParameterfv,
                     (GLenum target, GLenum pname, const GLfloat* params))
GST_GL_EXT_FUNCTION (void, TexParameteri,
                     (GLenum target, GLenum pname, GLint param))
GST_GL_EXT_FUNCTION (void, TexParameteriv,
                     (GLenum target, GLenum pname, const GLint* params))
GST_GL_EXT_FUNCTION (void, GetTexParameterfv,
                     (GLenum target, GLenum pname, GLfloat* params))
GST_GL_EXT_FUNCTION (void, GetTexParameteriv,
                     (GLenum target, GLenum pname, GLint* params))
GST_GL_EXT_FUNCTION (void, TexSubImage2D,
                     (GLenum target,
                      GLint level,
                      GLint xoffset,
                      GLint yoffset,
                      GLsizei width,
                      GLsizei height,
                      GLenum format,
                      GLenum type,
                      const GLvoid* pixels))
GST_GL_EXT_FUNCTION (void, CopyTexImage2D,
                     (GLenum target,
                      GLint level,
                      GLenum internalformat,
                      GLint x,
                      GLint y,
                      GLsizei width,
                      GLsizei height,
                      GLint border))
GST_GL_EXT_FUNCTION (void, Viewport,
                     (GLint x, GLint y, GLsizei width, GLsizei height))
GST_GL_EXT_FUNCTION (GLboolean, IsEnabled, (GLenum cap))
GST_GL_EXT_FUNCTION (void, LineWidth, (GLfloat width))
GST_GL_EXT_FUNCTION (void, PolygonOffset, (GLfloat factor, GLfloat units))
GST_GL_EXT_FUNCTION (void, TexParameterf,
                     (GLenum target,
                      GLenum pname,
                      GLfloat param))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (texture_3d,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3,
                  1, 2,
                  255, 255, /* not in either GLES */
                  "OES\0",
                  "texture_3D\0")
GST_GL_EXT_FUNCTION (void, TexImage3D,
                     (GLenum target, GLint level,
                      GLint internalFormat,
                      GLsizei width, GLsizei height,
                      GLsizei depth, GLint border,
                      GLenum format, GLenum type,
                      const GLvoid *pixels))
GST_GL_EXT_FUNCTION (void, TexSubImage3D,
                     (GLenum target, GLint level,
                      GLint xoffset, GLint yoffset,
                      GLint zoffset, GLsizei width,
                      GLsizei height, GLsizei depth,
                      GLenum format,
                      GLenum type, const GLvoid *pixels))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (only_in_both_gles_and_gl_1_3,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES1 | GST_GL_API_GLES2,
                  1, 3,
                  1, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, CompressedTexImage2D,
                     (GLenum target,
                      GLint level,
                      GLenum internalformat,
                      GLsizei width,
                      GLsizei height,
                      GLint border,
                      GLsizei imageSize,
                      const GLvoid* data))
GST_GL_EXT_FUNCTION (void, CompressedTexSubImage2D,
                     (GLenum target,
                      GLint level,
                      GLint xoffset,
                      GLint yoffset,
                      GLsizei width,
                      GLsizei height,
                      GLenum format,
                      GLsizei imageSize,
                      const GLvoid* data))
GST_GL_EXT_FUNCTION (void, SampleCoverage,
                     (GLclampf value, GLboolean invert))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (only_in_both_gles_and_gl_1_5,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES1 | GST_GL_API_GLES2,
                  1, 5,
                  1, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, GetBufferParameteriv,
                     (GLenum target, GLenum pname, GLint* params))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (vbos,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES1 | GST_GL_API_GLES2,
                  1, 5,
                  1, 0,
                  "ARB\0",
                  "vertex_buffer_object\0")
GST_GL_EXT_FUNCTION (void, GenBuffers,
                     (GLsizei		 n,
                      GLuint		*buffers))
GST_GL_EXT_FUNCTION (void, BindBuffer,
                     (GLenum		 target,
                      GLuint		 buffer))
GST_GL_EXT_FUNCTION (void, BufferData,
                     (GLenum		 target,
                      GLsizeiptr		 size,
                      const GLvoid		*data,
                      GLenum		 usage))
GST_GL_EXT_FUNCTION (void, BufferSubData,
                     (GLenum		 target,
                      GLintptr		 offset,
                      GLsizeiptr		 size,
                      const GLvoid		*data))
GST_GL_EXT_FUNCTION (void, DeleteBuffers,
                     (GLsizei		 n,
                      const GLuint		*buffers))
GST_GL_EXT_FUNCTION (GLboolean, IsBuffer,
                     (GLuint               buffer))
GST_GL_EXT_END ()

/* Available in GL 1.3, the multitexture extension or GLES. These are
   required */
GST_GL_EXT_BEGIN (multitexture_part0,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES1 | GST_GL_API_GLES2,
                  1, 3,
                  1, 0,
                  "ARB\0",
                  "multitexture\0")
GST_GL_EXT_FUNCTION (void, ActiveTexture,
                   (GLenum                texture))
GST_GL_EXT_END ()


 /* GLES doesn't support mapping buffers in core so this has to be a
   separate check */
GST_GL_EXT_BEGIN (map_vbos, 
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3,
                  1, 5,
                  255, 255, /* not in GLES core */
                  "ARB\0OES\0",
                  "vertex_buffer_object\0mapbuffer\0")
GST_GL_EXT_FUNCTION (void *, MapBuffer,
                   (GLenum		 target,
                    GLenum		 access))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (unmap_buffer,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 | GST_GL_API_GLES2,
                  1, 5,
                  3, 0,
                  "ARB\0OES\0",
                  "vertex_buffer_object\0mapbuffer\0")
GST_GL_EXT_FUNCTION (GLboolean, UnmapBuffer,
                   (GLenum		 target))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (gl3,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  3, 0,
                  3, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (const GLubyte*, GetStringi,
                     (GLenum name, GLint index))
GST_GL_EXT_FUNCTION (void *, MapBufferRange,
                     (GLenum target, GLintptr offset, GLsizeiptr length, GLbitfield access))
GST_GL_EXT_END ()
