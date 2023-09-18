/*
 * GStreamer
 * Copyright (C) 2015 Matthew Waters <matthew@centricular.com>
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

GST_GL_EXT_BEGIN (buffer_copy_sub_data,
                  GST_GL_API_OPENGL3 |
                  GST_GL_API_GLES2,
                  3, 1,
                  3, 0,
                  /* extension is available in GL 3.0 */
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, CopyBufferSubData,
                     (GLenum                readTarget,
                      GLenum                writeTarget,
                      GLintptr              readOffset,
                      GLintptr              writeOffset,
                      GLsizeiptr            size))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (get_buffer_sub_data,
                  GST_GL_API_OPENGL3,
                  1, 5,
                  255, 255,
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, GetBufferSubData,
                     (GLenum                target,
                      GLintptr              offset,
                      GLsizeiptr            size,
                      void *                data))
GST_GL_EXT_END ()
