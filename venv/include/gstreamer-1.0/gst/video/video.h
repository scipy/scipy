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

#ifndef __GST_VIDEO_H__
#define __GST_VIDEO_H__

#include <gst/gst.h>

#include <gst/video/video-prelude.h>

typedef struct _GstVideoAlignment GstVideoAlignment;

#include <gst/video/video-format.h>
#include <gst/video/video-color.h>
#include <gst/video/video-dither.h>
#include <gst/video/video-info.h>
#include <gst/video/video-frame.h>
#include <gst/video/video-enumtypes.h>
#include <gst/video/video-converter.h>
#include <gst/video/video-scaler.h>
#include <gst/video/video-multiview.h>

G_BEGIN_DECLS

/**
 * GstVideoAlignment:
 * @padding_left: extra pixels on the left side
 * @padding_right: extra pixels on the right side
 * @padding_top: extra pixels on the top
 * @padding_bottom: extra pixels on the bottom
 * @stride_align: array with extra alignment requirements for the strides
 *
 * Extra alignment parameters for the memory of video buffers. This
 * structure is usually used to configure the bufferpool if it supports the
 * #GST_BUFFER_POOL_OPTION_VIDEO_ALIGNMENT.
 */
struct _GstVideoAlignment
{
  guint padding_top;
  guint padding_bottom;
  guint padding_left;
  guint padding_right;
  guint stride_align[GST_VIDEO_MAX_PLANES];
};

/**
 * GstVideoOrientationMethod:
 * @GST_VIDEO_ORIENTATION_IDENTITY: Identity (no rotation)
 * @GST_VIDEO_ORIENTATION_90R: Rotate clockwise 90 degrees
 * @GST_VIDEO_ORIENTATION_180: Rotate 180 degrees
 * @GST_VIDEO_ORIENTATION_90L: Rotate counter-clockwise 90 degrees
 * @GST_VIDEO_ORIENTATION_HORIZ: Flip horizontally
 * @GST_VIDEO_ORIENTATION_VERT: Flip vertically
 * @GST_VIDEO_ORIENTATION_UL_LR: Flip across upper left/lower right diagonal
 * @GST_VIDEO_ORIENTATION_UR_LL: Flip across upper right/lower left diagonal
 * @GST_VIDEO_ORIENTATION_AUTO: Select flip method based on image-orientation tag
 * @GST_VIDEO_ORIENTATION_CUSTOM: Current status depends on plugin internal setup
 *
 * The different video orientation methods.
 *
 * Since: 1.10
 */
typedef enum {
  GST_VIDEO_ORIENTATION_IDENTITY,
  GST_VIDEO_ORIENTATION_90R,
  GST_VIDEO_ORIENTATION_180,
  GST_VIDEO_ORIENTATION_90L,
  GST_VIDEO_ORIENTATION_HORIZ,
  GST_VIDEO_ORIENTATION_VERT,
  GST_VIDEO_ORIENTATION_UL_LR,
  GST_VIDEO_ORIENTATION_UR_LL,
  GST_VIDEO_ORIENTATION_AUTO,
  GST_VIDEO_ORIENTATION_CUSTOM,
} GstVideoOrientationMethod;

/**
 * GST_TYPE_VIDEO_ORIENTATION_METHOD:
 *
 * Since: 1.20
 */

/* metadata macros */
/**
 * GST_META_TAG_VIDEO_STR:
 *
 * This metadata is relevant for video streams.
 *
 * Since: 1.2
 */
#define GST_META_TAG_VIDEO_STR "video"
/**
 * GST_META_TAG_VIDEO_ORIENTATION_STR:
 *
 * This metadata stays relevant as long as video orientation is unchanged.
 *
 * Since: 1.2
 */
#define GST_META_TAG_VIDEO_ORIENTATION_STR "orientation"
/**
 * GST_META_TAG_VIDEO_SIZE_STR:
 *
 * This metadata stays relevant as long as video size is unchanged.
 *
 * Since: 1.2
 */
#define GST_META_TAG_VIDEO_SIZE_STR "size"
/**
 * GST_META_TAG_VIDEO_COLORSPACE_STR:
 *
 * This metadata stays relevant as long as video colorspace is unchanged.
 *
 * Since: 1.2
 */
#define GST_META_TAG_VIDEO_COLORSPACE_STR "colorspace"

GST_VIDEO_API
void           gst_video_alignment_reset         (GstVideoAlignment *align);


/* some helper functions */

GST_VIDEO_API
gboolean       gst_video_calculate_display_ratio (guint * dar_n,
                                                  guint * dar_d,
                                                  guint   video_width,
                                                  guint   video_height,
                                                  guint   video_par_n,
                                                  guint   video_par_d,
                                                  guint   display_par_n,
                                                  guint   display_par_d);

GST_VIDEO_API
gboolean       gst_video_guess_framerate (GstClockTime duration,
                                          gint * dest_n, gint * dest_d);

GST_VIDEO_API
gboolean       gst_video_is_common_aspect_ratio (gint width, gint height, gint par_n, gint par_d);

/* convert/encode video sample from one format to another */

typedef void (*GstVideoConvertSampleCallback) (GstSample * sample, GError *error, gpointer user_data);

GST_VIDEO_API
void          gst_video_convert_sample_async (GstSample                    * sample,
                                              const GstCaps                * to_caps,
                                              GstClockTime                   timeout,
                                              GstVideoConvertSampleCallback  callback,
                                              gpointer                       user_data,
                                              GDestroyNotify                 destroy_notify);

GST_VIDEO_API
GstSample *   gst_video_convert_sample       (GstSample     * sample,
                                              const GstCaps * to_caps,
                                              GstClockTime    timeout,
                                              GError       ** error);


GST_VIDEO_API
gboolean gst_video_orientation_from_tag (GstTagList * taglist,
                                         GstVideoOrientationMethod * method);

G_END_DECLS

#include <gst/video/colorbalancechannel.h>
#include <gst/video/colorbalance.h>
#include <gst/video/gstvideoaffinetransformationmeta.h>
#include <gst/video/gstvideoaggregator.h>
#include <gst/video/gstvideocodecalphameta.h>
#include <gst/video/gstvideodecoder.h>
#include <gst/video/gstvideoencoder.h>
#include <gst/video/gstvideofilter.h>
#include <gst/video/gstvideometa.h>
#include <gst/video/gstvideopool.h>
#include <gst/video/gstvideosink.h>
#include <gst/video/gstvideotimecode.h>
#include <gst/video/gstvideoutils.h>
#include <gst/video/navigation.h>
#include <gst/video/video-anc.h>
#include <gst/video/video-blend.h>
#include <gst/video/videodirection.h>
#include <gst/video/video-event.h>
#include <gst/video/video-hdr.h>
#include <gst/video/videoorientation.h>
#include <gst/video/video-overlay-composition.h>
#include <gst/video/videooverlay.h>
#include <gst/video/video-sei.h>

#endif /* __GST_VIDEO_H__ */
