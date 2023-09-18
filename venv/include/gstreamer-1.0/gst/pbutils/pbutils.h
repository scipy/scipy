/* GStreamer base utils library
 * Copyright (C) 2006 Tim-Philipp Müller <tim centricular net>
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

#ifndef __GST_PB_UTILS_BASE_UTILS_H__
#define __GST_PB_UTILS_BASE_UTILS_H__

#include <gst/gst.h>

#include <gst/pbutils/gstpluginsbaseversion.h>
#include <gst/pbutils/descriptions.h>
#include <gst/pbutils/missing-plugins.h>
#include <gst/pbutils/install-plugins.h>
#include <gst/pbutils/codec-utils.h>
#include <gst/pbutils/pbutils-enumtypes.h>
#include <gst/pbutils/gstdiscoverer.h>
#include <gst/pbutils/encoding-profile.h>
#include <gst/pbutils/encoding-target.h>
#include <gst/pbutils/gstaudiovisualizer.h>

G_BEGIN_DECLS

GST_PBUTILS_API
void    gst_pb_utils_init (void);

G_END_DECLS

#endif /* __GST_PB_UTILS_BASE_UTILS_H__ */

