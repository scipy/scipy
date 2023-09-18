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

/* This lists functions that are unique to GL 2.0 or GLES 2.0 and are
 * not in the old GLSL extensions */
GST_GL_EXT_BEGIN (shaders_glsl_2_only,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 0,
                  2, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (GLuint, CreateProgram,
                     (void))
GST_GL_EXT_FUNCTION (GLuint, CreateShader,
                     (GLenum                shaderType))
GST_GL_EXT_FUNCTION (void, DeleteShader,
                     (GLuint                shader))
GST_GL_EXT_FUNCTION (void, AttachShader,
                     (GLuint                program,
                      GLuint                shader))
GST_GL_EXT_FUNCTION (void, UseProgram,
                     (GLuint                program))
GST_GL_EXT_FUNCTION (void, DeleteProgram,
                     (GLuint                program))
GST_GL_EXT_FUNCTION (void, GetShaderInfoLog,
                     (GLuint                shader,
                      GLsizei               maxLength,
                      GLsizei              *length,
                      char                 *infoLog))
GST_GL_EXT_FUNCTION (void, GetProgramInfoLog,
                     (GLuint                program,
                      GLsizei               bufSize,
                      GLsizei              *length,
                      char                 *infoLog))
GST_GL_EXT_FUNCTION (void, GetShaderiv,
                     (GLuint                shader,
                      GLenum                pname,
                      GLint                *params))
GST_GL_EXT_FUNCTION (void, GetProgramiv,
                     (GLuint                program,
                      GLenum                pname,
                      GLint                *params))
GST_GL_EXT_FUNCTION (void, DetachShader,
                     (GLuint program, GLuint shader))
GST_GL_EXT_FUNCTION (void, GetAttachedShaders,
                     (GLuint program,
                      GLsizei maxcount,
                      GLsizei* count,
                      GLuint* shaders))
GST_GL_EXT_FUNCTION (GLboolean, IsShader,
                     (GLuint shader))
GST_GL_EXT_FUNCTION (GLboolean, IsProgram,
                     (GLuint program))
GST_GL_EXT_END ()

/* These functions are provided by GL_ARB_shader_objects or are in GL
 * 2.0 core */
GST_GL_EXT_BEGIN (shader_objects_or_gl2,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 0,
                  2, 0,
                  "ARB\0",
                  "shader_objects\0")
GST_GL_EXT_FUNCTION (void, ShaderSource,
                     (GLuint                shader,
                      GLsizei               count,
                      const char          **string,
                      const GLint          *length))
GST_GL_EXT_FUNCTION (void, CompileShader,
                     (GLuint                shader))
GST_GL_EXT_FUNCTION (void, LinkProgram,
                     (GLuint                program))
GST_GL_EXT_FUNCTION (GLint, GetUniformLocation,
                     (GLuint                program,
                      const char           *name))
GST_GL_EXT_FUNCTION (void, Uniform1f,
                     (GLint                 location,
                      GLfloat               v0))
GST_GL_EXT_FUNCTION (void, Uniform2f,
                     (GLint                 location,
                      GLfloat               v0,
                      GLfloat               v1))
GST_GL_EXT_FUNCTION (void, Uniform3f,
                     (GLint                 location,
                      GLfloat               v0,
                      GLfloat               v1,
                      GLfloat               v2))
GST_GL_EXT_FUNCTION (void, Uniform4f,
                     (GLint                 location,
                      GLfloat               v0,
                      GLfloat               v1,
                      GLfloat               v2,
                      GLfloat               v3))
GST_GL_EXT_FUNCTION (void, Uniform1fv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLfloat *       value))
GST_GL_EXT_FUNCTION (void, Uniform2fv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLfloat *       value))
GST_GL_EXT_FUNCTION (void, Uniform3fv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLfloat *       value))
GST_GL_EXT_FUNCTION (void, Uniform4fv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLfloat *       value))
GST_GL_EXT_FUNCTION (void, Uniform1i,
                     (GLint                 location,
                      GLint                 v0))
GST_GL_EXT_FUNCTION (void, Uniform2i,
                     (GLint                 location,
                      GLint                 v0,
                      GLint                 v1))
GST_GL_EXT_FUNCTION (void, Uniform3i,
                     (GLint                 location,
                      GLint                 v0,
                      GLint                 v1,
                      GLint                 v2))
GST_GL_EXT_FUNCTION (void, Uniform4i,
                     (GLint                 location,
                      GLint                 v0,
                      GLint                 v1,
                      GLint                 v2,
                      GLint                 v3))
GST_GL_EXT_FUNCTION (void, Uniform1iv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLint *         value))
GST_GL_EXT_FUNCTION (void, Uniform2iv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLint *         value))
GST_GL_EXT_FUNCTION (void, Uniform3iv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLint *         value))
GST_GL_EXT_FUNCTION (void, Uniform4iv,
                     (GLint                 location,
                      GLsizei               count,
                      const GLint *         value))
GST_GL_EXT_FUNCTION (void, UniformMatrix2fv,
                     (GLint                 location,
                      GLsizei               count,
                      GLboolean             transpose,
                      const GLfloat        *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix3fv,
                     (GLint                 location,
                      GLsizei               count,
                      GLboolean             transpose,
                      const GLfloat        *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix4fv,
                     (GLint                 location,
                      GLsizei               count,
                      GLboolean             transpose,
                      const GLfloat        *value))

GST_GL_EXT_FUNCTION (void, GetUniformfv,
                     (GLuint                program,
                      GLint                 location,
                      GLfloat              *params))
GST_GL_EXT_FUNCTION (void, GetUniformiv,
                     (GLuint                program,
                      GLint                 location,
                      GLint                *params))
GST_GL_EXT_FUNCTION (void, GetActiveUniform,
                     (GLuint program,
                      GLuint index,
                      GLsizei bufsize,
                      GLsizei* length,
                      GLint* size,
                      GLenum* type,
                      GLchar* name))
GST_GL_EXT_FUNCTION (void, GetShaderSource,
                     (GLuint shader,
                      GLsizei bufsize,
                      GLsizei* length,
                      GLchar* source))
GST_GL_EXT_FUNCTION (void, ValidateProgram, (GLuint program))
GST_GL_EXT_END ()

/* These functions are provided by GL_ARB_vertex_shader or are in GL
 * 2.0 core */
GST_GL_EXT_BEGIN (vertex_shaders,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 0,
                  2, 0,
                  "ARB\0",
                  "vertex_shader\0")
GST_GL_EXT_FUNCTION (void, VertexAttribPointer,
                     (GLuint		 index,
                      GLint		 size,
                      GLenum		 type,
                      GLboolean		 normalized,
                      GLsizei		 stride,
                      const GLvoid        *pointer))
GST_GL_EXT_FUNCTION (void, EnableVertexAttribArray,
                     (GLuint		 index))
GST_GL_EXT_FUNCTION (void, DisableVertexAttribArray,
                     (GLuint		 index))
GST_GL_EXT_FUNCTION (void, VertexAttrib1f, (GLuint indx, GLfloat x))
GST_GL_EXT_FUNCTION (void, VertexAttrib1fv,
                     (GLuint indx, const GLfloat* values))
GST_GL_EXT_FUNCTION (void, VertexAttrib2f, (GLuint indx, GLfloat x, GLfloat y))
GST_GL_EXT_FUNCTION (void, VertexAttrib2fv,
                     (GLuint indx, const GLfloat* values))
GST_GL_EXT_FUNCTION (void, VertexAttrib3f,
                     (GLuint indx, GLfloat x, GLfloat y, GLfloat z))
GST_GL_EXT_FUNCTION (void, VertexAttrib3fv,
                     (GLuint indx, const GLfloat* values))
GST_GL_EXT_FUNCTION (void, VertexAttrib4f,
                     (GLuint index, GLfloat x, GLfloat y, GLfloat z, GLfloat w))
GST_GL_EXT_FUNCTION (void, VertexAttrib4fv,
                     (GLuint indx, const GLfloat* values))
GST_GL_EXT_FUNCTION (void, GetVertexAttribfv,
                     (GLuint index, GLenum pname, GLfloat* params))
GST_GL_EXT_FUNCTION (void, GetVertexAttribiv,
                     (GLuint index, GLenum pname, GLint* params))
GST_GL_EXT_FUNCTION (void, GetVertexAttribPointerv,
                     (GLuint index, GLenum pname, GLvoid** pointer))
GST_GL_EXT_FUNCTION (GLint, GetAttribLocation,
                     (GLuint program, const char *name))
GST_GL_EXT_FUNCTION (void, BindAttribLocation,
                     (GLuint program,
                      GLuint index,
                      const GLchar* name))
GST_GL_EXT_FUNCTION (void, GetActiveAttrib,
                     (GLuint program,
                      GLuint index,
                      GLsizei bufsize,
                      GLsizei* length,
                      GLint* size,
                      GLenum* type,
                      GLchar* name))
GST_GL_EXT_END ()

/* These only list functions that come from the old GLSL extensions.
 * Functions that are common to the extensions and GLSL 2.0 should
 * instead be listed in cogl-glsl-functions.h */
GST_GL_EXT_BEGIN (shader_objects,
                  GST_GL_API_NONE,
                  255, 255,
                  255, 255, /* not in either GLES */
                  "ARB\0",
                  "shader_objects\0")
GST_GL_EXT_FUNCTION (GLuint, CreateProgramObject,
                     (void))
GST_GL_EXT_FUNCTION (GLuint, CreateShaderObject,
                     (GLenum shaderType))
GST_GL_EXT_FUNCTION (void, DeleteObject,
                     (GLuint obj))
GST_GL_EXT_FUNCTION (void, AttachObject,
                     (GLuint container, GLuint obj))
GST_GL_EXT_FUNCTION (void, UseProgramObject,
                     (GLuint programObj))
GST_GL_EXT_FUNCTION (void, GetInfoLog,
                     (GLuint                obj,
                      GLsizei               maxLength,
                      GLsizei              *length,
                      char                 *infoLog))
GST_GL_EXT_FUNCTION (void, GetObjectParameteriv,
                     (GLuint                obj,
                      GLenum                pname,
                      GLint                *params))
GST_GL_EXT_FUNCTION (void, DetachObject,
                     (GLuint container, GLuint obj))
GST_GL_EXT_FUNCTION (void, GetAttachedObjects,
                     (GLuint program,
                      GLsizei maxcount,
                      GLsizei* count,
                      GLuint* shaders))
GST_GL_EXT_END ()

/* ARB_fragment_program */
GST_GL_EXT_BEGIN (arbfp,
                  GST_GL_API_NONE,
                  255, 255,
                  255, 255, /* not in either GLES */
                  "ARB\0",
                  "fragment_program\0")
GST_GL_EXT_FUNCTION (void, GenPrograms,
                     (GLsizei               n,
                      GLuint               *programs))
GST_GL_EXT_FUNCTION (void, DeletePrograms,
                     (GLsizei               n,
                      GLuint               *programs))
GST_GL_EXT_FUNCTION (void, BindProgram,
                     (GLenum                target,
                      GLuint                program))
GST_GL_EXT_FUNCTION (void, ProgramString,
                     (GLenum                target,
                      GLenum                format,
                      GLsizei               len,
                      const void           *program))
GST_GL_EXT_FUNCTION (void, ProgramLocalParameter4fv,
                     (GLenum                target,
                      GLuint                index,
                      GLfloat              *params))
GST_GL_EXT_END ()

/* This lists functions that are unique to GL 2.1 or GLES 3.0 and are
 * not in the old GLSL extensions */
GST_GL_EXT_BEGIN (shaders_2_1,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 1,
                  3, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, UniformMatrix2x3fv,
                     (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix3x2fv,
                     (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix2x4fv,
                     (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix4x2fv,
                     (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix3x4fv,
                     (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value))
GST_GL_EXT_FUNCTION (void, UniformMatrix4x3fv,
                     (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (bind_frag_data,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3,
                  3, 0,
                  255, 255,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, BindFragDataLocation,
                     (GLuint program, GLuint index, const GLchar * name))
GST_GL_EXT_END ()
