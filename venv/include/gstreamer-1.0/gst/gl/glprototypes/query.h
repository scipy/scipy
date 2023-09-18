/*
 * GStreamer
 * Copyright (C) 2016 Matthew Waters <matthew@centricular.com>
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

GST_GL_EXT_BEGIN (timer_query,
                  GST_GL_API_OPENGL3,
                  3, 3,
                  3, 0,
                  "ARB:\0ANGLE\0EXT\0",
                  "timer_query\0disjoint_timer_query\0")
GST_GL_EXT_FUNCTION (void, GenQueries,
                     (GLsizei n,
                     GLuint *ids))
GST_GL_EXT_FUNCTION (void, DeleteQueries,
                     (GLsizei n,
                      GLuint *ids))
GST_GL_EXT_FUNCTION (GLboolean, IsQuery,
                     (GLuint id))
GST_GL_EXT_FUNCTION (void, BeginQuery,
                     (GLenum target,
                      GLuint id))
GST_GL_EXT_FUNCTION (void, EndQuery,
                     (GLenum target))
GST_GL_EXT_FUNCTION (void, QueryCounter,
                     (GLuint id,
                     GLenum target))
GST_GL_EXT_FUNCTION (void, GetQueryiv,
                     (GLenum target,
                      GLenum pname,
                      GLint *params))
GST_GL_EXT_FUNCTION (void, GetQueryObjectiv,
                     (GLuint id,
                      GLenum pname,
                      GLint *params))
GST_GL_EXT_FUNCTION (void, GetQueryObjectuiv,
                     (GLuint id,
                      GLenum pname,
                      GLuint *params))
GST_GL_EXT_FUNCTION (void, GetQueryObjecti64v,
                     (GLuint id,
                      GLenum pname,
                      GLint64 *params))
GST_GL_EXT_FUNCTION (void, GetQueryObjectui64v,
                     (GLuint id,
                      GLenum pname,
                      GLuint64 *params))
GST_GL_EXT_END ()
