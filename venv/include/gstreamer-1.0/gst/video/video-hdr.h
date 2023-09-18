/* GStreamer
 * Copyright (C) <2018-2019> Seungha Yang <seungha.yang@navercorp.com>
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

#ifndef __GST_VIDEO_HDR_H__
#define __GST_VIDEO_HDR_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

typedef struct _GstVideoMasteringDisplayInfoCoordinates GstVideoMasteringDisplayInfoCoordinates;
typedef struct _GstVideoMasteringDisplayInfo GstVideoMasteringDisplayInfo;
typedef struct _GstVideoContentLightLevel GstVideoContentLightLevel;

/**
 * GstVideoMasteringDisplayInfoCoordinates:
 * @x: the x coordinate of CIE 1931 color space in unit of 0.00002.
 * @y: the y coordinate of CIE 1931 color space in unit of 0.00002.
 *
 * Used to represent display_primaries and white_point of
 * #GstVideoMasteringDisplayInfo struct. See #GstVideoMasteringDisplayInfo
 *
 * Since: 1.18
 */
struct _GstVideoMasteringDisplayInfoCoordinates
{
  guint16 x;
  guint16 y;
};

/**
 * GstVideoMasteringDisplayInfo:
 * @display_primaries: the xy coordinates of primaries in the CIE 1931 color space.
 *   the index 0 contains red, 1 is for green and 2 is for blue.
 *   each value is normalized to 50000 (meaning that in unit of 0.00002)
 * @white_point: the xy coordinates of white point in the CIE 1931 color space.
 *   each value is normalized to 50000 (meaning that in unit of 0.00002)
 * @max_display_mastering_luminance: the maximum value of display luminance
 *   in unit of 0.0001 candelas per square metre (cd/m^2 and nit)
 * @min_display_mastering_luminance: the minimum value of display luminance
 *   in unit of 0.0001 candelas per square metre (cd/m^2 and nit)
 *
 * Mastering display color volume information defined by SMPTE ST 2086
 * (a.k.a static HDR metadata).
 *
 * Since: 1.18
 */
struct _GstVideoMasteringDisplayInfo
{
  GstVideoMasteringDisplayInfoCoordinates display_primaries[3];
  GstVideoMasteringDisplayInfoCoordinates white_point;
  guint32 max_display_mastering_luminance;
  guint32 min_display_mastering_luminance;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_VIDEO_API
void      gst_video_mastering_display_info_init         (GstVideoMasteringDisplayInfo * minfo);

GST_VIDEO_API
gboolean  gst_video_mastering_display_info_from_string  (GstVideoMasteringDisplayInfo * minfo,
                                                         const gchar * mastering);

GST_VIDEO_API
gchar *   gst_video_mastering_display_info_to_string    (const GstVideoMasteringDisplayInfo * minfo);

GST_VIDEO_API
gboolean  gst_video_mastering_display_info_is_equal     (const GstVideoMasteringDisplayInfo * minfo,
                                                         const GstVideoMasteringDisplayInfo * other);

GST_VIDEO_API
gboolean  gst_video_mastering_display_info_from_caps    (GstVideoMasteringDisplayInfo * minfo,
                                                         const GstCaps * caps);

GST_VIDEO_API
gboolean  gst_video_mastering_display_info_add_to_caps  (const GstVideoMasteringDisplayInfo * minfo,
                                                         GstCaps * caps);

/**
 * GstVideoContentLightLevel:
 * @max_content_light_level: the maximum content light level
 *   (abbreviated to MaxCLL) in candelas per square meter (cd/m^2 and nit)
 * @max_frame_average_light_level: the maximum frame average light level
 *   (abbreviated to MaxFLL) in candelas per square meter (cd/m^2 and nit)
 *
 * Content light level information specified in CEA-861.3, Appendix A.
 *
 * Since: 1.18
 */
struct _GstVideoContentLightLevel
{
  guint16 max_content_light_level;
  guint16 max_frame_average_light_level;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_VIDEO_API
void      gst_video_content_light_level_init         (GstVideoContentLightLevel * linfo);

GST_VIDEO_API
gboolean  gst_video_content_light_level_from_string  (GstVideoContentLightLevel * linfo,
                                                      const gchar * level);

GST_VIDEO_API
gchar *   gst_video_content_light_level_to_string    (const GstVideoContentLightLevel * linfo);

GST_VIDEO_API
gboolean  gst_video_content_light_level_is_equal     (const GstVideoContentLightLevel * linfo,
                                                      const GstVideoContentLightLevel * other);

GST_VIDEO_API
gboolean  gst_video_content_light_level_from_caps    (GstVideoContentLightLevel * linfo,
                                                      const GstCaps * caps);

GST_VIDEO_API
gboolean  gst_video_content_light_level_add_to_caps  (const GstVideoContentLightLevel * linfo,
                                                      GstCaps * caps);


G_END_DECLS

#endif /* __GST_VIDEO_HDR_H__ */
