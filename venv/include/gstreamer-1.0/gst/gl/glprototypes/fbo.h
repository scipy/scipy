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

GST_GL_EXT_BEGIN (offscreen,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  3, 0,
                  2, 0,
                  /* for some reason the ARB version of this
                     extension doesn't have an ARB suffix for the
                     functions */
                  "ARB:\0EXT\0OES\0",
                  "framebuffer_object\0")
GST_GL_EXT_FUNCTION (void, GenRenderbuffers,
                     (GLsizei               n,
                      GLuint               *renderbuffers))
GST_GL_EXT_FUNCTION (void, DeleteRenderbuffers,
                     (GLsizei               n,
                      const GLuint         *renderbuffers))
GST_GL_EXT_FUNCTION (void, BindRenderbuffer,
                     (GLenum                target,
                      GLuint                renderbuffer))
GST_GL_EXT_FUNCTION (void, RenderbufferStorage,
                     (GLenum                target,
                      GLenum                internalformat,
                      GLsizei               width,
                      GLsizei               height))
GST_GL_EXT_FUNCTION (void, GenFramebuffers,
                     (GLsizei               n,
                      GLuint               *framebuffers))
GST_GL_EXT_FUNCTION (void, BindFramebuffer,
                     (GLenum                target,
                      GLuint                framebuffer))
GST_GL_EXT_FUNCTION (void, FramebufferTexture2D,
                     (GLenum                target,
                      GLenum                attachment,
                      GLenum                textarget,
                      GLuint                texture,
                      GLint                 level))
GST_GL_EXT_FUNCTION (void, FramebufferRenderbuffer,
                     (GLenum                target,
                      GLenum                attachment,
                      GLenum                renderbuffertarget,
                      GLuint                renderbuffer))
GST_GL_EXT_FUNCTION (GLboolean, IsRenderbuffer,
                     (GLuint                renderbuffer))
GST_GL_EXT_FUNCTION (GLenum, CheckFramebufferStatus,
                     (GLenum                target))
GST_GL_EXT_FUNCTION (void, DeleteFramebuffers,
                     (GLsizei               n,
                      const                 GLuint *framebuffers))
GST_GL_EXT_FUNCTION (void, GenerateMipmap,
                     (GLenum                target))
GST_GL_EXT_FUNCTION (void, GetFramebufferAttachmentParameteriv,
                     (GLenum                target,
                      GLenum                attachment,
                      GLenum                pname,
                      GLint                *params))
GST_GL_EXT_FUNCTION (void, GetRenderbufferParameteriv,
                     (GLenum                target,
                      GLenum                pname,
                      GLint                *params))
GST_GL_EXT_FUNCTION (GLboolean, IsFramebuffer,
                     (GLuint                framebuffer))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (offscreen_blit,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  3, 0,
                  3, 0,
                  "EXT\0ANGLE\0",
                  "framebuffer_blit\0")
GST_GL_EXT_FUNCTION (void, BlitFramebuffer,
                     (GLint                 srcX0,
                      GLint                 srcY0,
                      GLint                 srcX1,
                      GLint                 srcY1,
                      GLint                 dstX0,
                      GLint                 dstY0,
                      GLint                 dstX1,
                      GLint                 dstY1,
                      GLbitfield            mask,
                      GLenum                filter))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (framebuffer_discard, 
                  GST_GL_API_NONE,
                  255, 255,
                  255, 255, /* not in either GLES */
                  "EXT\0",
                  "framebuffer_discard\0")
GST_GL_EXT_FUNCTION (void, DiscardFramebuffer,
                     (GLenum           target,
                      GLsizei          numAttachments,
                      const GLenum    *attachments))
GST_GL_EXT_END ()


GST_GL_EXT_BEGIN (read_buffer,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  1, 0,
                  3, 0,
                  "NV\0",
                  "read_buffer\0")
GST_GL_EXT_FUNCTION (void, ReadBuffer,
                     (GLenum mode))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (draw_buffers,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  2, 1,
                  3, 0,
                  "ARB\0ATI\0NV\0",
                  "draw_buffers\0")
GST_GL_EXT_FUNCTION (void, DrawBuffers,
                     (GLsizei n, const GLenum *bufs))
GST_GL_EXT_END ()
