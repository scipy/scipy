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

/* These are the core GL functions which are only available in big
   GL */
GST_GL_EXT_BEGIN (only_in_big_gl,
                  GST_GL_API_OPENGL | GST_GL_API_OPENGL3,
                  1, 0,
                  255, 255, /* not in GLES */
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, GetTexLevelParameteriv,
                     (GLenum target, GLint level,
                      GLenum pname, GLint *params))
GST_GL_EXT_FUNCTION (void, GetTexImage,
                     (GLenum target, GLint level,
                      GLenum format, GLenum type,
                      GLvoid *pixels))
GST_GL_EXT_FUNCTION (void, DepthRange,
                     (double near_val, double far_val))
GST_GL_EXT_FUNCTION (void, DrawBuffer,
                     (GLenum mode))
GST_GL_EXT_FUNCTION (void, ClearDepth,
                     (double depth))
GST_GL_EXT_END ()

GST_GL_EXT_BEGIN (only_in_big_gl_compat,
                  GST_GL_API_OPENGL,
                  1, 0,
                  255, 255, /* not in GLES */
                  "\0",
                  "\0")
GST_GL_EXT_FUNCTION (void, ClipPlane,
                     (GLenum plane, const double *equation))
GST_GL_EXT_END ()
