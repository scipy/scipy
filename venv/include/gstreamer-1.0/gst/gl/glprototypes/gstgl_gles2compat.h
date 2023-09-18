/*
 * GStreamer
 * Copyright (C) 2009 Julien Isorce <julien.isorce@gmail.com>
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

/* Compatibility for OpenGL ES 2.0 */

#ifndef __GST_GL_ES2__
#define __GST_GL_ES2__

#include <glib.h>

G_BEGIN_DECLS

#ifndef GL_RGB16
#define GL_RGB16 GL_RGB565
#endif
#ifndef GL_RGB8
#define GL_RGB8 GL_RGB
#endif

/* UNSUPPORTED */

#ifndef GL_COLOR_ATTACHMENT1
#define GL_COLOR_ATTACHMENT1 0x8CE1
#endif
#ifndef GL_COLOR_ATTACHMENT2
#define GL_COLOR_ATTACHMENT2 0x8CE2
#endif
#ifndef GL_COLOR_ATTACHMENT3
#define GL_COLOR_ATTACHMENT3 0x8CE3
#endif
#ifndef GL_TEXTURE_ENV
#define GL_TEXTURE_ENV 0
#endif
#ifndef GL_TEXTURE_ENV_MODE
#define GL_TEXTURE_ENV_MODE 0
#endif
#ifndef GL_DEPTH24_STENCIL8
#define GL_DEPTH24_STENCIL8 0x88F0
#endif

G_END_DECLS

#endif /* __GST_GL_ES2__ */
