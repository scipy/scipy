/* GStreamer RIFF I/O
 * Copyright (C) 2003 Ronald Bultje <rbultje@ronald.bitfreak.net>
 *
 * riff-media.h: RIFF-id to/from caps routines
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

#ifndef __GST_RIFF_MEDIA_H__
#define __GST_RIFF_MEDIA_H__

#include <glib.h>
#include <gst/gst.h>
#include "riff-ids.h"

G_BEGIN_DECLS

/*
 * Create caos. strh/strf, strf/strd_data and codec_name can be NULL.
 */

GST_RIFF_API
GstCaps * gst_riff_create_video_caps (guint32              codec_fcc,
                                      gst_riff_strh      * strh,
                                      gst_riff_strf_vids * strf,
                                      GstBuffer          * strf_data,
                                      GstBuffer          * strd_data,
                                      char              ** codec_name);

GST_RIFF_API
GstCaps * gst_riff_create_audio_caps (guint16              codec_id,
                                      gst_riff_strh      * strh,
                                      gst_riff_strf_auds * strf,
                                      GstBuffer          * strf_data,
                                      GstBuffer          * strd_data,
                                      char              ** codec_name,
                                      gint                 channel_reorder_map[18]);

GST_RIFF_API
GstCaps * gst_riff_create_iavs_caps  (guint32              codec_fcc,
                                      gst_riff_strh      * strh,
                                      gst_riff_strf_iavs * strf,
                                      GstBuffer          * strf_data,
                                      GstBuffer          * strd_data,
                                      char              ** codec_name);
/*
 * Create template caps (includes all known types).
 */

GST_RIFF_API
GstCaps * gst_riff_create_video_template_caps (void);

GST_RIFF_API
GstCaps * gst_riff_create_audio_template_caps (void);

GST_RIFF_API
GstCaps * gst_riff_create_iavs_template_caps  (void);

G_END_DECLS

#endif /* __GST_RIFF_READ_H__ */
