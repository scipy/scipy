/*
 * GStreamer
 * Copyright (C) 2014 Matthew Waters <ystreet00@gmail.com>
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

GST_GL_EXT_BEGIN (debug,
                  GST_GL_API_OPENGL3,
                  4, 3,
                  255, 255,
                  "KHR:\0KHR\0ARB\0",
                  "debug\0debug_output\0")
GST_GL_EXT_FUNCTION (void, DebugMessageControl,
                     (GLenum source,
                      GLenum type,
                      GLenum severity,
                      GLsizei count,
                      const GLuint* ids,
                      gboolean enabled))
GST_GL_EXT_FUNCTION (void, DebugMessageInsert,
                     (GLenum source,
                      GLenum type,
                      GLuint id,
                      GLenum severity,
                      GLsizei length,
                      const gchar *message))
GST_GL_EXT_FUNCTION (void, DebugMessageCallback,
                     (GST_GL_DEBUG_PROC callback,
                      gpointer user_data))
GST_GL_EXT_FUNCTION (GLuint, GetDebugMessageLog,
                     (GLuint count,
                      GLsizei bufSize,
                      GLenum* sources,
                      GLenum* types,
                      GLuint* ids,
                      GLenum* severities,
                      GLsizei* lengths,
                      gchar* messageLog))
GST_GL_EXT_FUNCTION (void, GetPointerv,
                     (GLenum pname,
                      gpointer * params))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (khr_debug,
                  GST_GL_API_OPENGL3,
                  4, 3,
                  255, 255,
                  "KHR:\0KHR\0",
                  "debug\0")
GST_GL_EXT_FUNCTION (void, PushDebugGroup,
                     (GLenum source,
                      GLuint id,
                      GLsizei length,
                      const gchar * message))
GST_GL_EXT_FUNCTION (void, PopDebugGroup, (void))
GST_GL_EXT_FUNCTION (void, ObjectLabel,
                     (GLenum identifier,
                      GLuint name,
                      GLsizei length,
                      const gchar *label))
GST_GL_EXT_FUNCTION (void, GetObjectLabel,
                     (GLenum identifier,
                      GLuint name,
                      GLsizei bufSize, 
                      GLsizei *length,
                      gchar *label))
GST_GL_EXT_FUNCTION (void, ObjectPtrLabel,
                     (gpointer ptr,
                      GLsizei length,
                      const gchar *label))
GST_GL_EXT_FUNCTION (void, GetObjectPtrLabel,
                     (gpointer ptr,
                      GLsizei bufSize,
                      GLsizei *length,
                      gchar *label))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (ext_debug_marker,
                  GST_GL_API_NONE,
                  255, 255,
                  255, 255,
                  "EXT\0",
                  "debug_marker\0")
GST_GL_EXT_FUNCTION (void, InsertEventMarker,
                     (GLsizei length,
                      const gchar * message))
GST_GL_EXT_FUNCTION (void, PushGroupMarker,
                     (GLsizei length,
                      const gchar * message))
GST_GL_EXT_FUNCTION (void, PopGroupMarker,
                     (void))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (gremedy_string_marker,
                  GST_GL_API_NONE,
                  255, 255,
                  255, 255,
                  "GREMEDY\0",
                  "string_marker\0")
GST_GL_EXT_FUNCTION (void, StringMarker,
                     (GLsizei length,
                      const gchar * message))
GST_GL_EXT_END ()
