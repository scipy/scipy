/*
 * GStreamer
 * Copyright (C) 2014-2015 Matthew Waters <matthew@centricular.com>
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

GST_GL_EXT_BEGIN (sync,
                  GST_GL_API_OPENGL3 | GST_GL_API_OPENGL |
                  GST_GL_API_GLES2,
                  3, 2,
                  3, 0,
                  "APPLE\0",
                  "sync\0")
GST_GL_EXT_FUNCTION (GLsync, FenceSync,
                     (GLenum condition,
                      GLbitfield flags))
GST_GL_EXT_FUNCTION (GLboolean, IsSync,
                     (GLsync sync))
GST_GL_EXT_FUNCTION (void, DeleteSync,
                     (GLsync sync))
GST_GL_EXT_FUNCTION (GLenum, ClientWaitSync,
                     (GLsync sync,
                      GLbitfield flags,
                      GLuint64 timeout))
GST_GL_EXT_FUNCTION (void, WaitSync,
                     (GLsync sync,
                      GLbitfield flags,
                      GLuint64 timeout))
GST_GL_EXT_FUNCTION (void, GetSynciv,
                     (GLsync sync,
                      GLenum name,
                      GLsizei buf_size,
                      GLsizei *length,
                      GLint *values))
GST_GL_EXT_END ()
