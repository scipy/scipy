/*
 * GStreamer
 * Copyright (C) 2014 Matthew Waters <matthew@centricular.com>
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

GST_GL_EXT_BEGIN (vao,
                  GST_GL_API_OPENGL3 | GST_GL_API_GLES2,
                  3, 0,
                  3, 0,
                  "ARB:\0OES\0",
                  "vertex_array_object\0")
GST_GL_EXT_FUNCTION (void, GenVertexArrays,
                     (GLsizei               n,
                      GLuint               *arrays))
GST_GL_EXT_FUNCTION (void, DeleteVertexArrays,
                     (GLsizei               n,
                      GLuint               *arrays))
GST_GL_EXT_FUNCTION (void, BindVertexArray,
                     (GLuint                array))
GST_GL_EXT_FUNCTION (GLboolean, IsVertexArray,
                     (GLuint                array))
GST_GL_EXT_END ()
