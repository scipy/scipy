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

#ifndef __GST_VIDEO_INFO_H__
#define __GST_VIDEO_INFO_H__

#include <gst/gst.h>
#include <gst/video/video-format.h>
#include <gst/video/video-color.h>

G_BEGIN_DECLS

#include <gst/video/video-enumtypes.h>

typedef struct _GstVideoInfo GstVideoInfo;

/**
 * GST_CAPS_FEATURE_FORMAT_INTERLACED:
 *
 * Name of the caps feature indicating that the stream is interlaced.
 *
 * Currently it is only used for video with 'interlace-mode=alternate'
 * to ensure backwards compatibility for this new mode.
 * In this mode each buffer carries a single field of interlaced video.
 * @GST_VIDEO_BUFFER_FLAG_TOP_FIELD and @GST_VIDEO_BUFFER_FLAG_BOTTOM_FIELD
 * indicate whether the buffer carries a top or bottom field. The order of
 * buffers/fields in the stream and the timestamps on the buffers indicate the
 * temporal order of the fields.
 * Top and bottom fields are expected to alternate in this mode.
 * The frame rate in the caps still signals the frame rate, so the notional field
 * rate will be twice the frame rate from the caps
 * (see @GST_VIDEO_INFO_FIELD_RATE_N).
 *
 * Since: 1.16.
 */
#define GST_CAPS_FEATURE_FORMAT_INTERLACED "format:Interlaced"

/**
 * GstVideoInterlaceMode:
 * @GST_VIDEO_INTERLACE_MODE_PROGRESSIVE: all frames are progressive
 * @GST_VIDEO_INTERLACE_MODE_INTERLEAVED: 2 fields are interleaved in one video
 *     frame. Extra buffer flags describe the field order.
 * @GST_VIDEO_INTERLACE_MODE_MIXED: frames contains both interlaced and
 *     progressive video, the buffer flags describe the frame and fields.
 * @GST_VIDEO_INTERLACE_MODE_FIELDS: 2 fields are stored in one buffer, use the
 *     frame ID to get access to the required field. For multiview (the
 *     'views' property > 1) the fields of view N can be found at frame ID
 *     (N * 2) and (N * 2) + 1.
 *     Each field has only half the amount of lines as noted in the
 *     height property. This mode requires multiple GstVideoMeta metadata
 *     to describe the fields.
 * @GST_VIDEO_INTERLACE_MODE_ALTERNATE: 1 field is stored in one buffer,
 *     @GST_VIDEO_BUFFER_FLAG_TF or @GST_VIDEO_BUFFER_FLAG_BF indicates if
 *     the buffer is carrying the top or bottom field, respectively. The top and
 *     bottom buffers must alternate in the pipeline, with this mode
 *     (Since: 1.16).
 *
 * The possible values of the #GstVideoInterlaceMode describing the interlace
 * mode of the stream.
 */
typedef enum {
  GST_VIDEO_INTERLACE_MODE_PROGRESSIVE = 0,
  GST_VIDEO_INTERLACE_MODE_INTERLEAVED,
  GST_VIDEO_INTERLACE_MODE_MIXED,
  GST_VIDEO_INTERLACE_MODE_FIELDS,
  GST_VIDEO_INTERLACE_MODE_ALTERNATE,
} GstVideoInterlaceMode;

GST_VIDEO_API
const gchar *          gst_video_interlace_mode_to_string    (GstVideoInterlaceMode mode);

GST_VIDEO_API
GstVideoInterlaceMode  gst_video_interlace_mode_from_string  (const gchar * mode);

/**
 * GstVideoMultiviewMode:
 * @GST_VIDEO_MULTIVIEW_MODE_NONE: A special value indicating
 * no multiview information. Used in GstVideoInfo and other places to
 * indicate that no specific multiview handling has been requested or
 * provided. This value is never carried on caps.
 * @GST_VIDEO_MULTIVIEW_MODE_MONO: All frames are monoscopic.
 * @GST_VIDEO_MULTIVIEW_MODE_LEFT: All frames represent a left-eye view.
 * @GST_VIDEO_MULTIVIEW_MODE_RIGHT: All frames represent a right-eye view.
 * @GST_VIDEO_MULTIVIEW_MODE_SIDE_BY_SIDE: Left and right eye views are
 * provided in the left and right half of the frame respectively.
 * @GST_VIDEO_MULTIVIEW_MODE_SIDE_BY_SIDE_QUINCUNX: Left and right eye
 * views are provided in the left and right half of the frame, but
 * have been sampled using quincunx method, with half-pixel offset
 * between the 2 views.
 * @GST_VIDEO_MULTIVIEW_MODE_COLUMN_INTERLEAVED: Alternating vertical
 * columns of pixels represent the left and right eye view respectively.
 * @GST_VIDEO_MULTIVIEW_MODE_ROW_INTERLEAVED: Alternating horizontal
 * rows of pixels represent the left and right eye view respectively.
 * @GST_VIDEO_MULTIVIEW_MODE_TOP_BOTTOM: The top half of the frame
 * contains the left eye, and the bottom half the right eye.
 * @GST_VIDEO_MULTIVIEW_MODE_CHECKERBOARD: Pixels are arranged with
 * alternating pixels representing left and right eye views in a
 * checkerboard fashion.
 * @GST_VIDEO_MULTIVIEW_MODE_FRAME_BY_FRAME: Left and right eye views
 * are provided in separate frames alternately.
 * @GST_VIDEO_MULTIVIEW_MODE_MULTIVIEW_FRAME_BY_FRAME: Multiple
 * independent views are provided in separate frames in sequence.
 * This method only applies to raw video buffers at the moment.
 * Specific view identification is via the `GstVideoMultiviewMeta`
 * and #GstVideoMeta(s) on raw video buffers.
 * @GST_VIDEO_MULTIVIEW_MODE_SEPARATED: Multiple views are
 * provided as separate #GstMemory framebuffers attached to each
 * #GstBuffer, described by the `GstVideoMultiviewMeta`
 * and #GstVideoMeta(s)
 *
 * All possible stereoscopic 3D and multiview representations.
 * In conjunction with #GstVideoMultiviewFlags, describes how
 * multiview content is being transported in the stream.
 */
typedef enum {
  GST_VIDEO_MULTIVIEW_MODE_NONE = -1,
  GST_VIDEO_MULTIVIEW_MODE_MONO = 0,
  /* Single view modes */
  GST_VIDEO_MULTIVIEW_MODE_LEFT,
  GST_VIDEO_MULTIVIEW_MODE_RIGHT,
  /* Stereo view modes */
  GST_VIDEO_MULTIVIEW_MODE_SIDE_BY_SIDE,
  GST_VIDEO_MULTIVIEW_MODE_SIDE_BY_SIDE_QUINCUNX,
  GST_VIDEO_MULTIVIEW_MODE_COLUMN_INTERLEAVED,
  GST_VIDEO_MULTIVIEW_MODE_ROW_INTERLEAVED,
  GST_VIDEO_MULTIVIEW_MODE_TOP_BOTTOM,
  GST_VIDEO_MULTIVIEW_MODE_CHECKERBOARD,
  /* Padding for new frame packing modes */

  GST_VIDEO_MULTIVIEW_MODE_FRAME_BY_FRAME = 32,
  /* Multivew mode(s) */
  GST_VIDEO_MULTIVIEW_MODE_MULTIVIEW_FRAME_BY_FRAME,
  GST_VIDEO_MULTIVIEW_MODE_SEPARATED
  /* future expansion for annotated modes */
} GstVideoMultiviewMode;

/**
 * GstVideoMultiviewFramePacking:
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_NONE: A special value indicating
 * no frame packing info.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_MONO: All frames are monoscopic.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_LEFT: All frames represent a left-eye view.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_RIGHT: All frames represent a right-eye view.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_SIDE_BY_SIDE: Left and right eye views are
 * provided in the left and right half of the frame respectively.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_SIDE_BY_SIDE_QUINCUNX: Left and right eye
 * views are provided in the left and right half of the frame, but
 * have been sampled using quincunx method, with half-pixel offset
 * between the 2 views.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_COLUMN_INTERLEAVED: Alternating vertical
 * columns of pixels represent the left and right eye view respectively.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_ROW_INTERLEAVED: Alternating horizontal
 * rows of pixels represent the left and right eye view respectively.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_TOP_BOTTOM: The top half of the frame
 * contains the left eye, and the bottom half the right eye.
 * @GST_VIDEO_MULTIVIEW_FRAME_PACKING_CHECKERBOARD: Pixels are arranged with
 * alternating pixels representing left and right eye views in a
 * checkerboard fashion.
 *
 * #GstVideoMultiviewFramePacking represents the subset of #GstVideoMultiviewMode
 * values that can be applied to any video frame without needing extra metadata.
 * It can be used by elements that provide a property to override the
 * multiview interpretation of a video stream when the video doesn't contain
 * any markers.
 *
 * This enum is used (for example) on playbin, to re-interpret a played
 * video stream as a stereoscopic video. The individual enum values are
 * equivalent to and have the same value as the matching #GstVideoMultiviewMode.
 *
 */
typedef enum {
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_NONE = GST_VIDEO_MULTIVIEW_MODE_NONE,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_MONO = GST_VIDEO_MULTIVIEW_MODE_MONO,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_LEFT = GST_VIDEO_MULTIVIEW_MODE_LEFT,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_RIGHT = GST_VIDEO_MULTIVIEW_MODE_RIGHT,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_SIDE_BY_SIDE = GST_VIDEO_MULTIVIEW_MODE_SIDE_BY_SIDE,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_SIDE_BY_SIDE_QUINCUNX = GST_VIDEO_MULTIVIEW_MODE_SIDE_BY_SIDE_QUINCUNX,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_COLUMN_INTERLEAVED = GST_VIDEO_MULTIVIEW_MODE_COLUMN_INTERLEAVED,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_ROW_INTERLEAVED = GST_VIDEO_MULTIVIEW_MODE_ROW_INTERLEAVED,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_TOP_BOTTOM = GST_VIDEO_MULTIVIEW_MODE_TOP_BOTTOM,
  GST_VIDEO_MULTIVIEW_FRAME_PACKING_CHECKERBOARD = GST_VIDEO_MULTIVIEW_MODE_CHECKERBOARD
} GstVideoMultiviewFramePacking;

#define GST_VIDEO_MULTIVIEW_MAX_FRAME_PACKING GST_VIDEO_MULTIVIEW_FRAME_PACKING_CHECKERBOARD

/**
 * GstVideoMultiviewFlags:
 * @GST_VIDEO_MULTIVIEW_FLAGS_NONE: No flags
 * @GST_VIDEO_MULTIVIEW_FLAGS_RIGHT_VIEW_FIRST: For stereo streams, the
 *     normal arrangement of left and right views is reversed.
 * @GST_VIDEO_MULTIVIEW_FLAGS_LEFT_FLIPPED: The left view is vertically
 *     mirrored.
 * @GST_VIDEO_MULTIVIEW_FLAGS_LEFT_FLOPPED: The left view is horizontally
 *     mirrored.
 * @GST_VIDEO_MULTIVIEW_FLAGS_RIGHT_FLIPPED: The right view is
 *     vertically mirrored.
 * @GST_VIDEO_MULTIVIEW_FLAGS_RIGHT_FLOPPED: The right view is
 *     horizontally mirrored.
 * @GST_VIDEO_MULTIVIEW_FLAGS_HALF_ASPECT: For frame-packed
 *     multiview modes, indicates that the individual
 *     views have been encoded with half the true width or height
 *     and should be scaled back up for display. This flag
 *     is used for overriding input layout interpretation
 *     by adjusting pixel-aspect-ratio.
 *     For side-by-side, column interleaved or checkerboard packings, the
 *     pixel width will be doubled. For row interleaved and top-bottom
 *     encodings, pixel height will be doubled.
 * @GST_VIDEO_MULTIVIEW_FLAGS_MIXED_MONO: The video stream contains both
 *     mono and multiview portions, signalled on each buffer by the
 *     absence or presence of the @GST_VIDEO_BUFFER_FLAG_MULTIPLE_VIEW
 *     buffer flag.
 *
 * GstVideoMultiviewFlags are used to indicate extra properties of a
 * stereo/multiview stream beyond the frame layout and buffer mapping
 * that is conveyed in the #GstVideoMultiviewMode.
 */
typedef enum {
  GST_VIDEO_MULTIVIEW_FLAGS_NONE             = 0,
  GST_VIDEO_MULTIVIEW_FLAGS_RIGHT_VIEW_FIRST = (1 << 0),
  GST_VIDEO_MULTIVIEW_FLAGS_LEFT_FLIPPED     = (1 << 1),
  GST_VIDEO_MULTIVIEW_FLAGS_LEFT_FLOPPED     = (1 << 2),
  GST_VIDEO_MULTIVIEW_FLAGS_RIGHT_FLIPPED    = (1 << 3),
  GST_VIDEO_MULTIVIEW_FLAGS_RIGHT_FLOPPED    = (1 << 4),
  GST_VIDEO_MULTIVIEW_FLAGS_HALF_ASPECT      = (1 << 14),
  GST_VIDEO_MULTIVIEW_FLAGS_MIXED_MONO       = (1 << 15)
} GstVideoMultiviewFlags;

/**
 * GstVideoFlags:
 * @GST_VIDEO_FLAG_NONE: no flags
 * @GST_VIDEO_FLAG_VARIABLE_FPS: a variable fps is selected, fps_n and fps_d
 *     denote the maximum fps of the video
 * @GST_VIDEO_FLAG_PREMULTIPLIED_ALPHA: Each color has been scaled by the alpha
 *     value.
 *
 * Extra video flags
 */
typedef enum {
  GST_VIDEO_FLAG_NONE                = 0,
  GST_VIDEO_FLAG_VARIABLE_FPS        = (1 << 0),
  GST_VIDEO_FLAG_PREMULTIPLIED_ALPHA = (1 << 1)
} GstVideoFlags;

/**
 * GstVideoFieldOrder:
 * @GST_VIDEO_FIELD_ORDER_UNKNOWN: unknown field order for interlaced content.
 *     The actual field order is signalled via buffer flags.
 * @GST_VIDEO_FIELD_ORDER_TOP_FIELD_FIRST: top field is first
 * @GST_VIDEO_FIELD_ORDER_BOTTOM_FIELD_FIRST: bottom field is first
 *
 * Field order of interlaced content. This is only valid for
 * interlace-mode=interleaved and not interlace-mode=mixed. In the case of
 * mixed or GST_VIDEO_FIELD_ORDER_UNKOWN, the field order is signalled via
 * buffer flags.
 *
 * Since: 1.12
 */
typedef enum {
  GST_VIDEO_FIELD_ORDER_UNKNOWN            = 0,
  GST_VIDEO_FIELD_ORDER_TOP_FIELD_FIRST    = 1,
  GST_VIDEO_FIELD_ORDER_BOTTOM_FIELD_FIRST = 2,
} GstVideoFieldOrder;

GST_VIDEO_API
const gchar *      gst_video_field_order_to_string    (GstVideoFieldOrder order);

GST_VIDEO_API
GstVideoFieldOrder gst_video_field_order_from_string  (const gchar * order);

/**
 * GstVideoInfo:
 * @finfo: the format info of the video
 * @interlace_mode: the interlace mode
 * @flags: additional video flags
 * @width: the width of the video
 * @height: the height of the video
 * @views: the number of views for multiview video
 * @size: the default size of one frame
 * @chroma_site: a #GstVideoChromaSite.
 * @colorimetry: the colorimetry info
 * @par_n: the pixel-aspect-ratio numerator
 * @par_d: the pixel-aspect-ratio denominator
 * @fps_n: the framerate numerator
 * @fps_d: the framerate denominator
 * @offset: offsets of the planes
 * @stride: strides of the planes
 * @multiview_mode: delivery mode for multiple views. (Since: 1.6)
 * @multiview_flags: flags for multiple views configuration (Since: 1.6)
 *
 * Information describing image properties. This information can be filled
 * in from GstCaps with gst_video_info_from_caps(). The information is also used
 * to store the specific video info when mapping a video frame with
 * gst_video_frame_map().
 *
 * Use the provided macros to access the info in this structure.
 */
struct _GstVideoInfo {
  const GstVideoFormatInfo *finfo;

  GstVideoInterlaceMode     interlace_mode;
  GstVideoFlags             flags;
  gint                      width;
  gint                      height;
  gsize                     size;
  gint                      views;

  GstVideoChromaSite        chroma_site;
  GstVideoColorimetry       colorimetry;

  gint                      par_n;
  gint                      par_d;
  gint                      fps_n;
  gint                      fps_d;

  gsize                     offset[GST_VIDEO_MAX_PLANES];
  gint                      stride[GST_VIDEO_MAX_PLANES];

  /* Union preserves padded struct size for backwards compat
   * Consumer code should use the accessor macros for fields */
  union {
    struct { /* < skip > */
      GstVideoMultiviewMode     multiview_mode;
      GstVideoMultiviewFlags    multiview_flags;
      GstVideoFieldOrder        field_order;
    } abi;
    /*< private >*/
    gpointer _gst_reserved[GST_PADDING];
  } ABI;
};

#define GST_TYPE_VIDEO_INFO              (gst_video_info_get_type ())
GST_VIDEO_API
GType gst_video_info_get_type            (void);

/* general info */
#define GST_VIDEO_INFO_FORMAT(i)         (GST_VIDEO_FORMAT_INFO_FORMAT((i)->finfo))
#define GST_VIDEO_INFO_NAME(i)           (GST_VIDEO_FORMAT_INFO_NAME((i)->finfo))
#define GST_VIDEO_INFO_IS_YUV(i)         (GST_VIDEO_FORMAT_INFO_IS_YUV((i)->finfo))
#define GST_VIDEO_INFO_IS_RGB(i)         (GST_VIDEO_FORMAT_INFO_IS_RGB((i)->finfo))
#define GST_VIDEO_INFO_IS_GRAY(i)        (GST_VIDEO_FORMAT_INFO_IS_GRAY((i)->finfo))
#define GST_VIDEO_INFO_HAS_ALPHA(i)      (GST_VIDEO_FORMAT_INFO_HAS_ALPHA((i)->finfo))

#define GST_VIDEO_INFO_INTERLACE_MODE(i) ((i)->interlace_mode)
#define GST_VIDEO_INFO_IS_INTERLACED(i)  ((i)->interlace_mode != GST_VIDEO_INTERLACE_MODE_PROGRESSIVE)
#define GST_VIDEO_INFO_FIELD_ORDER(i)    ((i)->ABI.abi.field_order)
#define GST_VIDEO_INFO_FLAGS(i)          ((i)->flags)
#define GST_VIDEO_INFO_WIDTH(i)          ((i)->width)
#define GST_VIDEO_INFO_HEIGHT(i)         ((i)->height)
/**
 * GST_VIDEO_INFO_FIELD_HEIGHT:
 *
 * The height of a field. It's the height of the full frame unless split-field
 * (alternate) interlacing is in use.
 *
 * Since: 1.16.
 */
#define GST_VIDEO_INFO_FIELD_HEIGHT(i)   ((i)->interlace_mode == GST_VIDEO_INTERLACE_MODE_ALTERNATE? GST_ROUND_UP_2 ((i)->height) / 2 : (i)->height)
#define GST_VIDEO_INFO_SIZE(i)           ((i)->size)
#define GST_VIDEO_INFO_VIEWS(i)          ((i)->views)
#define GST_VIDEO_INFO_PAR_N(i)          ((i)->par_n)
#define GST_VIDEO_INFO_PAR_D(i)          ((i)->par_d)
#define GST_VIDEO_INFO_FPS_N(i)          ((i)->fps_n)
#define GST_VIDEO_INFO_FIELD_RATE_N(i)   ((GST_VIDEO_INFO_INTERLACE_MODE ((i)) == \
                                           GST_VIDEO_INTERLACE_MODE_ALTERNATE) ? \
                                           (i)->fps_n * 2 : (i)->fps_n)
#define GST_VIDEO_INFO_FPS_D(i)          ((i)->fps_d)

#define GST_VIDEO_INFO_COLORIMETRY(i) ((i)->colorimetry)
#define GST_VIDEO_INFO_CHROMA_SITE(i) ((i)->chroma_site)

#define GST_VIDEO_INFO_MULTIVIEW_MODE(i)          ((i)->ABI.abi.multiview_mode)
#define GST_VIDEO_INFO_MULTIVIEW_FLAGS(i)          ((i)->ABI.abi.multiview_flags)

/* dealing with GstVideoInfo flags */
#define GST_VIDEO_INFO_FLAG_IS_SET(i,flag) ((GST_VIDEO_INFO_FLAGS(i) & (flag)) == (flag))
#define GST_VIDEO_INFO_FLAG_SET(i,flag)    (GST_VIDEO_INFO_FLAGS(i) |= (flag))
#define GST_VIDEO_INFO_FLAG_UNSET(i,flag)  (GST_VIDEO_INFO_FLAGS(i) &= ~(flag))

/* dealing with planes */
#define GST_VIDEO_INFO_N_PLANES(i)       (GST_VIDEO_FORMAT_INFO_N_PLANES((i)->finfo))
#define GST_VIDEO_INFO_PLANE_OFFSET(i,p) ((i)->offset[p])
#define GST_VIDEO_INFO_PLANE_STRIDE(i,p) ((i)->stride[p])
/**
 * GST_VIDEO_INFO_PLANE_HEIGHT:
 *
 * The padded height in pixels of a plane (padded size divided by the plane stride).
 * In case of GST_VIDEO_INTERLACE_MODE_ALTERNATE info, this macro returns the
 * plane heights used to hold a single field, not the full frame.
 *
 * The size passed as third argument is the size of the pixel data and should
 * not contain any extra metadata padding.
 *
 * It is not valid to use this macro with a TILED format.
 *
 * Since: 1.18
 */
#define GST_VIDEO_INFO_PLANE_HEIGHT(i,p,sizes) ((i)->stride[p] == 0 ? 0 : sizes[p] / (i)->stride[p])

/* dealing with components */
#define GST_VIDEO_INFO_N_COMPONENTS(i)   GST_VIDEO_FORMAT_INFO_N_COMPONENTS((i)->finfo)
#define GST_VIDEO_INFO_COMP_DEPTH(i,c)   GST_VIDEO_FORMAT_INFO_DEPTH((i)->finfo,(c))
#define GST_VIDEO_INFO_COMP_DATA(i,d,c)  GST_VIDEO_FORMAT_INFO_DATA((i)->finfo,d,(c))
#define GST_VIDEO_INFO_COMP_OFFSET(i,c)  GST_VIDEO_FORMAT_INFO_OFFSET((i)->finfo,(i)->offset,(c))
#define GST_VIDEO_INFO_COMP_STRIDE(i,c)  GST_VIDEO_FORMAT_INFO_STRIDE((i)->finfo,(i)->stride,(c))
#define GST_VIDEO_INFO_COMP_WIDTH(i,c)   GST_VIDEO_FORMAT_INFO_SCALE_WIDTH((i)->finfo,(c),(i)->width)
#define GST_VIDEO_INFO_COMP_HEIGHT(i,c)  GST_VIDEO_FORMAT_INFO_SCALE_HEIGHT((i)->finfo,(c),GST_VIDEO_INFO_FIELD_HEIGHT(i))
#define GST_VIDEO_INFO_COMP_PLANE(i,c)   GST_VIDEO_FORMAT_INFO_PLANE((i)->finfo,(c))
#define GST_VIDEO_INFO_COMP_PSTRIDE(i,c) GST_VIDEO_FORMAT_INFO_PSTRIDE((i)->finfo,(c))
#define GST_VIDEO_INFO_COMP_POFFSET(i,c) GST_VIDEO_FORMAT_INFO_POFFSET((i)->finfo,(c))

GST_VIDEO_API
GstVideoInfo * gst_video_info_new         (void);

GST_VIDEO_API
void           gst_video_info_init        (GstVideoInfo *info);

GST_VIDEO_API
GstVideoInfo * gst_video_info_copy        (const GstVideoInfo *info);

GST_VIDEO_API
void           gst_video_info_free        (GstVideoInfo *info);

GST_VIDEO_API
GstVideoInfo * gst_video_info_new_from_caps (const GstCaps * caps);

GST_VIDEO_API
gboolean       gst_video_info_set_format  (GstVideoInfo *info, GstVideoFormat format,
                                           guint width, guint height);

GST_VIDEO_API
gboolean       gst_video_info_set_interlaced_format
                                          (GstVideoInfo         *info,
                                           GstVideoFormat        format,
                                           GstVideoInterlaceMode mode,
                                           guint                 width,
                                           guint                 height);

GST_VIDEO_API
gboolean       gst_video_info_from_caps   (GstVideoInfo *info, const GstCaps  * caps);

GST_VIDEO_API
GstCaps *      gst_video_info_to_caps     (const GstVideoInfo *info);

GST_VIDEO_API
gboolean       gst_video_info_convert     (const GstVideoInfo *info,
                                           GstFormat     src_format,
                                           gint64        src_value,
                                           GstFormat     dest_format,
                                           gint64       *dest_value);

GST_VIDEO_API
gboolean       gst_video_info_is_equal    (const GstVideoInfo *info,
                                           const GstVideoInfo *other);

#include <gst/video/video.h>

GST_VIDEO_API
gboolean       gst_video_info_align       (GstVideoInfo * info, GstVideoAlignment * align);

GST_VIDEO_API
gboolean       gst_video_info_align_full  (GstVideoInfo * info, GstVideoAlignment * align, gsize plane_size[GST_VIDEO_MAX_PLANES]);


G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstVideoInfo, gst_video_info_free)

G_END_DECLS

#endif /* __GST_VIDEO_INFO_H__ */
