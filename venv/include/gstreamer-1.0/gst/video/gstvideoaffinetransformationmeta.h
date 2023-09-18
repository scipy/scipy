/* GStreamer
 * Copyright (C) Collabora Ltd.
 *   Author: Matthieu Bouron <matthieu.bouron@collabora.com>
 * Copyright (C) 2015, Matthew Waters <matthew@centricular.com>
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

#ifndef __GST_VIDEO_AFFINE_TRANSFORMATION_META_H__
#define __GST_VIDEO_AFFINE_TRANSFORMATION_META_H__

#include <gst/gst.h>
#include <gst/video/video.h>

G_BEGIN_DECLS

#define GST_VIDEO_AFFINE_TRANSFORMATION_META_API_TYPE (gst_video_affine_transformation_meta_api_get_type())
#define GST_VIDEO_AFFINE_TRANSFORMATION_META_INFO  (gst_video_affine_transformation_meta_get_info())

typedef struct _GstVideoAffineTransformationMeta GstVideoAffineTransformationMeta;
typedef gboolean (*GstVideoAffineTransformationGetMatrix) (GstVideoAffineTransformationMeta * meta, gfloat * matrix);

#define GST_CAPS_FEATURE_META_GST_VIDEO_AFFINE_TRANSFORMATION_META "meta:GstVideoAffineTransformation"
#define GST_BUFFER_POOL_OPTION_VIDEO_AFFINE_TRANSFORMATION_META "GstBufferPoolOptionVideoAffineTransformation"

/**
 * GstVideoAffineTransformationMeta:
 * @meta: parent #GstMeta
 * @matrix: the column-major 4x4 transformation matrix
 *
 * Extra buffer metadata for performing an affine transformation using a 4x4
 * matrix. The transformation matrix can be composed with
 * gst_video_affine_transformation_meta_apply_matrix().
 *
 * The vertices operated on are all in the range 0 to 1, not in
 * Normalized Device Coordinates (-1 to +1). Transforming points in this space
 * are assumed to have an origin at (0.5, 0.5, 0.5) in a left-handed coordinate
 * system with the x-axis moving horizontally (positive values to the right),
 * the y-axis moving vertically (positive values up the screen) and the z-axis
 * perpendicular to the screen (positive values into the screen).
 *
 * Since: 1.8
 */
struct _GstVideoAffineTransformationMeta
{
  GstMeta meta;

  gfloat matrix[16];
};

GST_VIDEO_API
GType gst_video_affine_transformation_meta_api_get_type          (void);

GST_VIDEO_API
const GstMetaInfo *gst_video_affine_transformation_meta_get_info (void);

#define gst_buffer_get_video_affine_transformation_meta(b) \
    ((GstVideoAffineTransformationMeta *)gst_buffer_get_meta((b),GST_VIDEO_AFFINE_TRANSFORMATION_META_API_TYPE))
GST_VIDEO_API
GstVideoAffineTransformationMeta *gst_buffer_add_video_affine_transformation_meta (GstBuffer * buffer);

GST_VIDEO_API
void gst_video_affine_transformation_meta_apply_matrix                           (GstVideoAffineTransformationMeta * meta,
                                                                                  const gfloat matrix[16]);

G_END_DECLS

#endif /* __GST_VIDEO_AFFINE_TRANSFORMATION_META_H__ */
