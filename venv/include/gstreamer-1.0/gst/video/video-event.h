/* GStreamer
 * Copyright (C) <2011> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_VIDEO_EVENT_H__
#define __GST_VIDEO_EVENT_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

/* video still frame event creation and parsing */

GST_VIDEO_API
GstEvent *     gst_video_event_new_still_frame   (gboolean in_still);

GST_VIDEO_API
gboolean       gst_video_event_parse_still_frame (GstEvent * event, gboolean * in_still);

/* video force key unit event creation and parsing */

GST_VIDEO_API
GstEvent * gst_video_event_new_downstream_force_key_unit (GstClockTime timestamp,
                                                          GstClockTime stream_time,
                                                          GstClockTime running_time,
                                                          gboolean all_headers,
                                                          guint count);

GST_VIDEO_API
gboolean gst_video_event_parse_downstream_force_key_unit (GstEvent * event,
                                                          GstClockTime * timestamp,
                                                          GstClockTime * stream_time,
                                                          GstClockTime * running_time,
                                                          gboolean * all_headers,
                                                          guint * count);

GST_VIDEO_API
GstEvent * gst_video_event_new_upstream_force_key_unit (GstClockTime running_time,
                                                        gboolean all_headers,
                                                        guint count);

GST_VIDEO_API
gboolean gst_video_event_parse_upstream_force_key_unit (GstEvent * event,
                                                        GstClockTime * running_time,
                                                        gboolean * all_headers,
                                                        guint * count);

GST_VIDEO_API
gboolean gst_video_event_is_force_key_unit(GstEvent *event);

G_END_DECLS

#endif /* __GST_VIDEO_EVENT_H__ */
