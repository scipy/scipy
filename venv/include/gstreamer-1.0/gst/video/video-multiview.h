/* GStreamer
 * Copyright (C) <2015> Jan Schmidt <jan@centricular.com>
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

#ifndef __GST_VIDEO_MULTIVIEW_H__
#define __GST_VIDEO_MULTIVIEW_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

/**
 * GstVideoMultiviewFlagsSet:
 *
 * See #GstVideoMultiviewFlags.
 */
#define GST_TYPE_VIDEO_MULTIVIEW_FLAGSET (gst_video_multiview_flagset_get_type())
GST_VIDEO_API
GType gst_video_multiview_flagset_get_type (void);

GST_VIDEO_API
const gchar * gst_video_multiview_mode_to_caps_string (GstVideoMultiviewMode mview_mode);

GST_VIDEO_API
GstVideoMultiviewMode gst_video_multiview_mode_from_caps_string (const gchar * caps_mview_mode);

GST_VIDEO_API
const GValue *gst_video_multiview_get_mono_modes(void);

GST_VIDEO_API
const GValue *gst_video_multiview_get_unpacked_modes(void);

GST_VIDEO_API
const GValue *gst_video_multiview_get_doubled_height_modes(void);

GST_VIDEO_API
const GValue *gst_video_multiview_get_doubled_width_modes(void);

GST_VIDEO_API
const GValue *gst_video_multiview_get_doubled_size_modes(void);

GST_VIDEO_API
void gst_video_multiview_video_info_change_mode (GstVideoInfo *info,
    GstVideoMultiviewMode out_mview_mode, GstVideoMultiviewFlags out_mview_flags);

GST_VIDEO_API
gboolean gst_video_multiview_guess_half_aspect (GstVideoMultiviewMode mv_mode,
    guint width, guint height, guint par_n, guint par_d);


#if 0 /* Place-holder for later MVC support */
#define GST_VIDEO_MULTIVIEW_META_API_TYPE (gst_video_multiview_meta_api_get_type())
#define GST_VIDEO_MULTIVIEW_META_INFO  (gst_video_multiview_meta_get_info())

typedef struct _GstVideoMultiviewMeta GstVideoMultiviewMeta;
typedef struct _GstVideoMultiviewViewInfo GstVideoMultiviewViewInfo;

GType gst_video_multiview_meta_api_get_type (void);
const GstMetaInfo * gst_video_multiview_meta_get_info (void);

GstVideoMultiviewMeta * gst_buffer_add_video_multiview_meta (GstBuffer *buffer, guint n_views);
#define gst_buffer_get_video_multiview_meta(b) ((GstVideoMultiviewMeta *)gst_buffer_get_meta((b),GST_VIDEO_MULTIVIEW_META_API_TYPE))

void gst_video_multiview_meta_set_n_views (GstVideoMultiviewMeta *mview_meta, guint n_views);

typedef enum {
  GST_VIDEO_MULTIVIEW_VIEW_UNKNOWN = 0,
  GST_VIDEO_MULTIVIEW_VIEW_MONO = 1,
  GST_VIDEO_MULTIVIEW_VIEW_LEFT = 2,
  GST_VIDEO_MULTIVIEW_VIEW_RIGHT = 3
} GstVideoMultiviewViewLabel;

struct _GstVideoMultiviewViewInfo {
  GstVideoMultiviewViewLabel view_label;

  guint meta_id; /* id of the GstVideoMeta for this view */

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

struct _GstVideoMultiviewMeta {
  GstMeta       meta;

  guint         n_views;
  GstVideoMultiviewViewInfo *view_info;
};
#endif

G_END_DECLS

#endif
