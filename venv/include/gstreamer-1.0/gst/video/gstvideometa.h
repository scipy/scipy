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

#ifndef __GST_VIDEO_META_H__
#define __GST_VIDEO_META_H__

#include <gst/gst.h>

#include <gst/video/video.h>
#include <gst/video/gstvideotimecode.h>

G_BEGIN_DECLS

#define GST_VIDEO_META_API_TYPE (gst_video_meta_api_get_type())
#define GST_VIDEO_META_INFO  (gst_video_meta_get_info())
typedef struct _GstVideoMeta GstVideoMeta;

#define GST_CAPS_FEATURE_META_GST_VIDEO_META "meta:GstVideoMeta"

#define GST_VIDEO_CROP_META_API_TYPE  (gst_video_crop_meta_api_get_type())
#define GST_VIDEO_CROP_META_INFO  (gst_video_crop_meta_get_info())
typedef struct _GstVideoCropMeta GstVideoCropMeta;

/**
 * GstVideoMeta:
 * @meta: parent #GstMeta
 * @buffer: the buffer this metadata belongs to
 * @flags: additional video flags
 * @format: the video format
 * @id: identifier of the frame
 * @width: the video width
 * @height: the video height
 * @n_planes: the number of planes in the image
 * @offset: array of offsets for the planes. This field might not always be
 *          valid, it is used by the default implementation of @map.
 * @stride: array of strides for the planes. This field might not always be
 *          valid, it is used by the default implementation of @map.
 * @map: map the memory of a plane
 * @unmap: unmap the memory of a plane
 * @alignment: the paddings and alignment constraints of the video buffer.
 * It is up to the caller of `gst_buffer_add_video_meta_full()` to set it
 * using gst_video_meta_set_alignment(), if they did not it defaults
 * to no padding and no alignment. Since: 1.18
 *
 * Extra buffer metadata describing image properties
 *
 * This meta can also be used by downstream elements to specifiy their
 * buffer layout requirements for upstream. Upstream should try to
 * fit those requirements, if possible, in order to prevent buffer copies.
 *
 * This is done by passing a custom #GstStructure to
 * gst_query_add_allocation_meta() when handling the ALLOCATION query.
 * This structure should be named 'video-meta' and can have the following
 * fields:
 * - padding-top (uint): extra pixels on the top
 * - padding-bottom (uint): extra pixels on the bottom
 * - padding-left (uint): extra pixels on the left side
 * - padding-right (uint): extra pixels on the right side
 * The padding fields have the same semantic as #GstVideoMeta.alignment
 * and so represent the paddings requested on produced video buffers.
 */
struct _GstVideoMeta {
  GstMeta            meta;

  GstBuffer         *buffer;

  GstVideoFrameFlags flags;
  GstVideoFormat     format;
  gint               id;
  guint              width;
  guint              height;

  guint              n_planes;
  gsize              offset[GST_VIDEO_MAX_PLANES];
  gint               stride[GST_VIDEO_MAX_PLANES];

  gboolean (*map)    (GstVideoMeta *meta, guint plane, GstMapInfo *info,
                      gpointer *data, gint * stride, GstMapFlags flags);
  gboolean (*unmap)  (GstVideoMeta *meta, guint plane, GstMapInfo *info);

  GstVideoAlignment  alignment;
};

GST_VIDEO_API
GType gst_video_meta_api_get_type (void);

GST_VIDEO_API
const GstMetaInfo * gst_video_meta_get_info (void);

GST_VIDEO_API
GstVideoMeta * gst_buffer_get_video_meta (GstBuffer *buffer);

GST_VIDEO_API
GstVideoMeta * gst_buffer_get_video_meta_id    (GstBuffer *buffer, gint id);

GST_VIDEO_API
GstVideoMeta * gst_buffer_add_video_meta       (GstBuffer *buffer, GstVideoFrameFlags flags,
                                                GstVideoFormat format, guint width, guint height);

GST_VIDEO_API
GstVideoMeta * gst_buffer_add_video_meta_full  (GstBuffer *buffer, GstVideoFrameFlags flags,
                                                GstVideoFormat format, guint width, guint height,
                                                guint n_planes, gsize offset[GST_VIDEO_MAX_PLANES],
                                                gint stride[GST_VIDEO_MAX_PLANES]);

GST_VIDEO_API
gboolean       gst_video_meta_map        (GstVideoMeta *meta, guint plane, GstMapInfo *info,
                                          gpointer *data, gint *stride, GstMapFlags flags);

GST_VIDEO_API
gboolean       gst_video_meta_unmap      (GstVideoMeta *meta, guint plane, GstMapInfo *info);

GST_VIDEO_API
gboolean       gst_video_meta_set_alignment (GstVideoMeta * meta, GstVideoAlignment alignment);

GST_VIDEO_API
gboolean       gst_video_meta_get_plane_size (GstVideoMeta * meta, gsize plane_size[GST_VIDEO_MAX_PLANES]);

GST_VIDEO_API
gboolean       gst_video_meta_get_plane_height (GstVideoMeta * meta, guint plane_height[GST_VIDEO_MAX_PLANES]);

/**
 * GstVideoCropMeta:
 * @meta: parent #GstMeta
 * @x: the horizontal offset
 * @y: the vertical offset
 * @width: the cropped width
 * @height: the cropped height
 *
 * Extra buffer metadata describing image cropping.
 */
struct _GstVideoCropMeta {
  GstMeta       meta;

  guint         x;
  guint         y;
  guint         width;
  guint         height;
};

GST_VIDEO_API
GType gst_video_crop_meta_api_get_type (void);

GST_VIDEO_API
const GstMetaInfo * gst_video_crop_meta_get_info (void);

#define gst_buffer_get_video_crop_meta(b) ((GstVideoCropMeta*)gst_buffer_get_meta((b),GST_VIDEO_CROP_META_API_TYPE))
#define gst_buffer_add_video_crop_meta(b) ((GstVideoCropMeta*)gst_buffer_add_meta((b),GST_VIDEO_CROP_META_INFO, NULL))

/* video metadata transforms */

GST_VIDEO_API
GQuark gst_video_meta_transform_scale_get_quark (void);
/**
 * gst_video_meta_transform_scale:
 *
 * GQuark for the video "gst-video-scale" transform.
 */
#define GST_VIDEO_META_TRANSFORM_IS_SCALE(type) ((type) == gst_video_meta_transform_scale_get_quark())

/**
 * GstVideoMetaTransform:
 * @in_info: the input #GstVideoInfo
 * @out_info: the output #GstVideoInfo
 *
 * Extra data passed to a video transform #GstMetaTransformFunction such as:
 * "gst-video-scale".
 */
typedef struct {
  GstVideoInfo *in_info;
  GstVideoInfo *out_info;
} GstVideoMetaTransform;

/**
 * GstVideoGLTextureType:
 * @GST_VIDEO_GL_TEXTURE_TYPE_LUMINANCE: Luminance texture, GL_LUMINANCE
 * @GST_VIDEO_GL_TEXTURE_TYPE_LUMINANCE_ALPHA: Luminance-alpha texture, GL_LUMINANCE_ALPHA
 * @GST_VIDEO_GL_TEXTURE_TYPE_RGB16: RGB 565 texture, GL_RGB
 * @GST_VIDEO_GL_TEXTURE_TYPE_RGB: RGB texture, GL_RGB
 * @GST_VIDEO_GL_TEXTURE_TYPE_RGBA: RGBA texture, GL_RGBA
 * @GST_VIDEO_GL_TEXTURE_TYPE_R: R texture, GL_RED_EXT
 * @GST_VIDEO_GL_TEXTURE_TYPE_RG: RG texture, GL_RG_EXT
 *
 * The GL texture type.
 */
typedef enum
{
  GST_VIDEO_GL_TEXTURE_TYPE_LUMINANCE,
  GST_VIDEO_GL_TEXTURE_TYPE_LUMINANCE_ALPHA,
  GST_VIDEO_GL_TEXTURE_TYPE_RGB16,
  GST_VIDEO_GL_TEXTURE_TYPE_RGB,
  GST_VIDEO_GL_TEXTURE_TYPE_RGBA,
  GST_VIDEO_GL_TEXTURE_TYPE_R,
  GST_VIDEO_GL_TEXTURE_TYPE_RG
} GstVideoGLTextureType;

/**
 * GstVideoGLTextureOrientation:
 * @GST_VIDEO_GL_TEXTURE_ORIENTATION_X_NORMAL_Y_NORMAL: Top line first in memory, left row first
 * @GST_VIDEO_GL_TEXTURE_ORIENTATION_X_NORMAL_Y_FLIP: Bottom line first in memory, left row first
 * @GST_VIDEO_GL_TEXTURE_ORIENTATION_X_FLIP_Y_NORMAL: Top line first in memory, right row first
 * @GST_VIDEO_GL_TEXTURE_ORIENTATION_X_FLIP_Y_FLIP: Bottom line first in memory, right row first
 *
 * The orientation of the GL texture.
 */
typedef enum
{
  GST_VIDEO_GL_TEXTURE_ORIENTATION_X_NORMAL_Y_NORMAL,
  GST_VIDEO_GL_TEXTURE_ORIENTATION_X_NORMAL_Y_FLIP,
  GST_VIDEO_GL_TEXTURE_ORIENTATION_X_FLIP_Y_NORMAL,
  GST_VIDEO_GL_TEXTURE_ORIENTATION_X_FLIP_Y_FLIP
} GstVideoGLTextureOrientation;

#define GST_VIDEO_GL_TEXTURE_UPLOAD_META_API_TYPE (gst_video_gl_texture_upload_meta_api_get_type())
#define GST_VIDEO_GL_TEXTURE_UPLOAD_META_INFO  (gst_video_gl_texture_upload_meta_get_info())

typedef struct _GstVideoGLTextureUploadMeta GstVideoGLTextureUploadMeta;
typedef gboolean (*GstVideoGLTextureUpload) (GstVideoGLTextureUploadMeta *meta, guint texture_id[4]);

#define GST_CAPS_FEATURE_META_GST_VIDEO_GL_TEXTURE_UPLOAD_META "meta:GstVideoGLTextureUploadMeta"

/**
 * GST_BUFFER_POOL_OPTION_VIDEO_GL_TEXTURE_UPLOAD_META:
 *
 * An option that can be activated on a bufferpool to request gl texture upload
 * meta on buffers from the pool.
 *
 * When this option is enabled on the bufferpool,
 * @GST_BUFFER_POOL_OPTION_VIDEO_META should also be enabled.
 *
 * Since: 1.2.2
 */
#define GST_BUFFER_POOL_OPTION_VIDEO_GL_TEXTURE_UPLOAD_META "GstBufferPoolOptionVideoGLTextureUploadMeta"

/**
 * GstVideoGLTextureUploadMeta:
 * @meta: parent #GstMeta
 * @texture_orientation: Orientation of the textures
 * @n_textures: Number of textures that are generated
 * @texture_type: Type of each texture
 *
 * Extra buffer metadata for uploading a buffer to an OpenGL texture
 * ID. The caller of gst_video_gl_texture_upload_meta_upload() must
 * have OpenGL set up and call this from a thread where it is valid
 * to upload something to an OpenGL texture.
 */

struct _GstVideoGLTextureUploadMeta {
  GstMeta       meta;

  GstVideoGLTextureOrientation texture_orientation;
  guint n_textures;
  GstVideoGLTextureType texture_type[4];

  /* <private> */
  GstBuffer *buffer;
  GstVideoGLTextureUpload upload;

  gpointer      user_data;
  GBoxedCopyFunc user_data_copy;
  GBoxedFreeFunc user_data_free;
};

GST_VIDEO_API
GType gst_video_gl_texture_upload_meta_api_get_type (void);

GST_VIDEO_API
const GstMetaInfo * gst_video_gl_texture_upload_meta_get_info (void);

#define gst_buffer_get_video_gl_texture_upload_meta(b) ((GstVideoGLTextureUploadMeta*)gst_buffer_get_meta((b),GST_VIDEO_GL_TEXTURE_UPLOAD_META_API_TYPE))

GST_VIDEO_API
GstVideoGLTextureUploadMeta *
          gst_buffer_add_video_gl_texture_upload_meta (GstBuffer *buffer,
                                                       GstVideoGLTextureOrientation texture_orientation,
                                                       guint n_textures,
                                                       GstVideoGLTextureType texture_type[4],
                                                       GstVideoGLTextureUpload upload,
                                                       gpointer user_data,
                                                       GBoxedCopyFunc user_data_copy,
                                                       GBoxedFreeFunc user_data_free);

GST_VIDEO_API
gboolean  gst_video_gl_texture_upload_meta_upload     (GstVideoGLTextureUploadMeta *meta,
                                                       guint texture_id[4]);


/**
 * GstVideoRegionOfInterestMeta:
 * @meta: parent #GstMeta
 * @roi_type: GQuark describing the semantic of the Roi (f.i. a face, a pedestrian)
 * @id: identifier of this particular ROI
 * @parent_id: identifier of its parent ROI, used f.i. for ROI hierarchisation.
 * @x: x component of upper-left corner
 * @y: y component of upper-left corner
 * @w: bounding box width
 * @h: bounding box height
 * @params: list of #GstStructure containing element-specific params for downstream,
 *          see gst_video_region_of_interest_meta_add_param(). (Since: 1.14)
 *
 * Extra buffer metadata describing an image region of interest
 */
typedef struct {
  GstMeta meta;

  GQuark roi_type;
  gint id;
  gint parent_id;

  guint x;
  guint y;
  guint w;
  guint h;

  GList *params;
} GstVideoRegionOfInterestMeta;

GST_VIDEO_API
GType              gst_video_region_of_interest_meta_api_get_type (void);
#define GST_VIDEO_REGION_OF_INTEREST_META_API_TYPE (gst_video_region_of_interest_meta_api_get_type())
GST_VIDEO_API
const GstMetaInfo *gst_video_region_of_interest_meta_get_info (void);
#define GST_VIDEO_REGION_OF_INTEREST_META_INFO (gst_video_region_of_interest_meta_get_info())

#define gst_buffer_get_video_region_of_interest_meta(b) \
        ((GstVideoRegionOfInterestMeta*)gst_buffer_get_meta((b),GST_VIDEO_REGION_OF_INTEREST_META_API_TYPE))
GST_VIDEO_API
GstVideoRegionOfInterestMeta *gst_buffer_get_video_region_of_interest_meta_id (GstBuffer   * buffer,
                                                                               gint          id);

GST_VIDEO_API
GstVideoRegionOfInterestMeta *gst_buffer_add_video_region_of_interest_meta    (GstBuffer   * buffer,
									       const gchar * roi_type,
									       guint         x,
                                                                               guint         y,
                                                                               guint         w,
                                                                               guint         h);

GST_VIDEO_API
GstVideoRegionOfInterestMeta *gst_buffer_add_video_region_of_interest_meta_id (GstBuffer   * buffer, 
									       GQuark        roi_type,
                                                                               guint         x,
                                                                               guint         y,
                                                                               guint         w,
                                                                               guint         h);
GST_VIDEO_API
void gst_video_region_of_interest_meta_add_param (GstVideoRegionOfInterestMeta * meta,
                                                  GstStructure * s);

GST_VIDEO_API
GstStructure *gst_video_region_of_interest_meta_get_param (GstVideoRegionOfInterestMeta * meta,
                                                           const gchar * name);

/**
 * GstVideoTimeCodeMeta:
 * @meta: parent #GstMeta
 * @tc: the GstVideoTimeCode to attach
 *
 * Extra buffer metadata describing the GstVideoTimeCode of the frame.
 *
 * Each frame is assumed to have its own timecode, i.e. they are not
 * automatically incremented/interpolated.
 *
 * Since: 1.10
 */
typedef struct {
  GstMeta meta;

  GstVideoTimeCode tc;
} GstVideoTimeCodeMeta;

GST_VIDEO_API
GType              gst_video_time_code_meta_api_get_type (void);
#define GST_VIDEO_TIME_CODE_META_API_TYPE (gst_video_time_code_meta_api_get_type())

GST_VIDEO_API
const GstMetaInfo *gst_video_time_code_meta_get_info (void);
#define GST_VIDEO_TIME_CODE_META_INFO (gst_video_time_code_meta_get_info())

#define gst_buffer_get_video_time_code_meta(b) \
        ((GstVideoTimeCodeMeta*)gst_buffer_get_meta((b),GST_VIDEO_TIME_CODE_META_API_TYPE))

GST_VIDEO_API
GstVideoTimeCodeMeta *gst_buffer_add_video_time_code_meta    (GstBuffer             * buffer,
                                                              const GstVideoTimeCode* tc);

GST_VIDEO_API
GstVideoTimeCodeMeta *
gst_buffer_add_video_time_code_meta_full                     (GstBuffer             * buffer,
                                                              guint fps_n,
                                                              guint fps_d,
                                                              GDateTime             * latest_daily_jam,
                                                              GstVideoTimeCodeFlags   flags,
                                                              guint                   hours,
                                                              guint                   minutes,
                                                              guint                   seconds,
                                                              guint                   frames,
                                                              guint                   field_count);

G_END_DECLS

#endif /* __GST_VIDEO_META_H__ */
