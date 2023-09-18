/* GStreamer
 * Copyright (C) <2005> Philippe Khalaf <burger@speedy.org>
 *               <2005> Wim Taymans <wim@fluendo.com>
 *
 * gstrtpbuffer.h: various helper functions to manipulate buffers
 *     with RTP payload.
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

#ifndef __GST_RTPDEFS_H__
#define __GST_RTPDEFS_H__

#include <gst/gst.h>
#include <gst/rtp/rtp-prelude.h>

/**
 * SECTION:gstrtpdefs
 * @title: GstRTPdefs
 * @short_description: common RTP defines
 *
 * Provides common defines for the RTP library.
 */

/**
 * GstRTPProfile:
 * @GST_RTP_PROFILE_UNKNOWN: invalid profile
 * @GST_RTP_PROFILE_AVP: the Audio/Visual profile (RFC 3551)
 * @GST_RTP_PROFILE_SAVP: the secure Audio/Visual profile (RFC 3711)
 * @GST_RTP_PROFILE_AVPF: the Audio/Visual profile with feedback (RFC 4585)
 * @GST_RTP_PROFILE_SAVPF: the secure Audio/Visual profile with feedback (RFC 5124)
 *
 * The transfer profile to use.
 *
 * Since: 1.6
 */
typedef enum {
  GST_RTP_PROFILE_UNKNOWN = 0,
  GST_RTP_PROFILE_AVP,
  GST_RTP_PROFILE_SAVP,
  GST_RTP_PROFILE_AVPF,
  GST_RTP_PROFILE_SAVPF
} GstRTPProfile;

#endif /* __GST_RTPDEFS_H__ */
