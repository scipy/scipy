/*
 * GStreamer
 * Copyright (C) 2021 Matthew Waters <matthew@centricular.com>
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

GST_GL_EXT_BEGIN (buffer_storage,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  4, 4,
                  255, 255,
                  "EXT\0ARB:\0", /* ARB version doesn't have function suffixes */
                  "buffer_storage\0")
GST_GL_EXT_FUNCTION (void, BufferStorage,
                     (GLenum target, GLsizeiptr, const void * data, GLbitfield flags))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (flush_mapped,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  3, 0,
                  3, 0,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, FlushMappedBufferRange,
                     (GLenum target, GLintptr offset, GLsizeiptr length))
GST_GL_EXT_END ()

