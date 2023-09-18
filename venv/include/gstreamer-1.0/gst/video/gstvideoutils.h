/* GStreamer
 * Copyright (C) 2008 David Schleef <ds@schleef.org>
 * Copyright (C) 2012 Collabora Ltd.
 *	Author : Edward Hervey <edward@collabora.com>
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
#include <gst/video/video.h>
#endif

#ifndef _GST_VIDEO_UTILS_H_
#define _GST_VIDEO_UTILS_H_

#include <gst/gst.h>
#include <gst/video/video-prelude.h>
#include <gst/video/video-hdr.h>

G_BEGIN_DECLS
#define GST_TYPE_VIDEO_CODEC_STATE \
  (gst_video_codec_state_get_type())

#define GST_TYPE_VIDEO_CODEC_FRAME \
  (gst_video_codec_frame_get_type())

typedef struct _GstVideoCodecState GstVideoCodecState;
typedef struct _GstVideoCodecFrame GstVideoCodecFrame;

/**
 * GstVideoCodecState:
 * @info: The #GstVideoInfo describing the stream
 * @caps: The #GstCaps used in the caps negotiation of the pad.
 * @codec_data: a #GstBuffer corresponding to the
 *     'codec_data' field of a stream, or NULL.
 * @allocation_caps: The #GstCaps for allocation query and pool
 *     negotiation. Since: 1.10
 * @mastering_display_info: Mastering display color volume information
 *     (HDR metadata) for the stream. Since: 1.20
 * @content_light_level: Content light level information for the stream.
 *     Since: 1.20
 *
 * Structure representing the state of an incoming or outgoing video
 * stream for encoders and decoders.
 *
 * Decoders and encoders will receive such a state through their
 * respective @set_format vmethods.
 *
 * Decoders and encoders can set the downstream state, by using the
 * gst_video_decoder_set_output_state() or
 * gst_video_encoder_set_output_state() methods.
 */
/**
 * GstVideoCodecState.mastering_display_info:
 *
 * Mastering display color volume information (HDR metadata) for the stream.
 *
 * Since: 1.20
 */
/**
 * GstVideoCodecState.content_light_level:
 *
 * Content light level information for the stream.
 *
 * Since: 1.20
 */
struct _GstVideoCodecState
{
  /*< private >*/
  gint ref_count;

  /*< public >*/
  GstVideoInfo info;

  GstCaps *caps;

  GstBuffer *codec_data;

  GstCaps *allocation_caps;

  GstVideoMasteringDisplayInfo *mastering_display_info;
  GstVideoContentLightLevel *content_light_level;

  /*< private >*/
  gpointer padding[GST_PADDING_LARGE - 3];
};

/**
 * GstVideoCodecFrameFlags:
 * @GST_VIDEO_CODEC_FRAME_FLAG_DECODE_ONLY: is the frame only meant to be decoded
 * @GST_VIDEO_CODEC_FRAME_FLAG_SYNC_POINT: is the frame a synchronization point (keyframe)
 * @GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME: should the output frame be made a keyframe
 * @GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME_HEADERS: should the encoder output stream headers
 * @GST_VIDEO_CODEC_FRAME_FLAG_CORRUPTED: the buffer data is corrupted (Since: 1.20)
 *
 * Flags for #GstVideoCodecFrame
 */
typedef enum
{
  GST_VIDEO_CODEC_FRAME_FLAG_DECODE_ONLY            = (1<<0),
  GST_VIDEO_CODEC_FRAME_FLAG_SYNC_POINT             = (1<<1),
  GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME         = (1<<2),
  GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME_HEADERS = (1<<3),
  /**
   * GST_VIDEO_CODEC_FRAME_FLAG_CORRUPTED:
   *
   * The buffer data is corrupted.
   *
   * Since: 1.20
   */
  GST_VIDEO_CODEC_FRAME_FLAG_CORRUPTED = (1<<4),
} GstVideoCodecFrameFlags;

/**
 * GST_VIDEO_CODEC_FRAME_FLAGS:
 * @frame: a #GstVideoCodecFrame
 *
 * The entire set of flags for the @frame
 */
#define GST_VIDEO_CODEC_FRAME_FLAGS(frame) ((frame)->flags)

/**
 * GST_VIDEO_CODEC_FRAME_FLAG_IS_SET:
 * @frame: a #GstVideoCodecFrame
 * @flag: a flag to check for
 *
 * Checks whether the given @flag is set
 */
#define GST_VIDEO_CODEC_FRAME_FLAG_IS_SET(frame,flag)   !!(GST_VIDEO_CODEC_FRAME_FLAGS(frame) & (flag))

/**
 * GST_VIDEO_CODEC_FRAME_FLAG_SET:
 * @frame: a #GstVideoCodecFrame
 * @flag: Flag to set, can be any number of bits in guint32.
 *
 * This macro sets the given bits
 */
#define GST_VIDEO_CODEC_FRAME_FLAG_SET(frame,flag)     (GST_VIDEO_CODEC_FRAME_FLAGS(frame) |= (flag))

/**
 * GST_VIDEO_CODEC_FRAME_FLAG_UNSET:
 * @frame: a #GstVideoCodecFrame
 * @flag: Flag to unset
 *
 * This macro usets the given bits.
 */
#define GST_VIDEO_CODEC_FRAME_FLAG_UNSET(frame,flag)   (GST_VIDEO_CODEC_FRAME_FLAGS(frame) &= ~(flag))

/**
 * GST_VIDEO_CODEC_FRAME_IS_DECODE_ONLY:
 * @frame: a #GstVideoCodecFrame
 *
 * Tests if the buffer should only be decoded but not sent downstream.
 */
#define GST_VIDEO_CODEC_FRAME_IS_DECODE_ONLY(frame)     (GST_VIDEO_CODEC_FRAME_FLAG_IS_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_DECODE_ONLY))

/**
 * GST_VIDEO_CODEC_FRAME_SET_DECODE_ONLY:
 * @frame: a #GstVideoCodecFrame
 *
 * Sets the buffer to not be sent downstream.
 *
 * Decoder implementation can use this if they have frames that
 * are not meant to be displayed.
 *
 * Encoder implementation can safely ignore this field.
 */
#define GST_VIDEO_CODEC_FRAME_SET_DECODE_ONLY(frame)    (GST_VIDEO_CODEC_FRAME_FLAG_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_DECODE_ONLY))

/**
 * GST_VIDEO_CODEC_FRAME_IS_SYNC_POINT:
 * @frame: a #GstVideoCodecFrame
 *
 * Tests if the frame is a synchronization point (like a keyframe).
 *
 * Decoder implementations can use this to detect keyframes.
 */
#define GST_VIDEO_CODEC_FRAME_IS_SYNC_POINT(frame)      (GST_VIDEO_CODEC_FRAME_FLAG_IS_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_SYNC_POINT))

/**
 * GST_VIDEO_CODEC_FRAME_SET_SYNC_POINT:
 * @frame: a #GstVideoCodecFrame
 *
 * Sets the frame to be a synchronization point (like a keyframe).
 *
 * Encoder implementations should set this accordingly.
 *
 * Decoder implementing parsing features should set this when they
 * detect such a synchronization point.
 */
#define GST_VIDEO_CODEC_FRAME_SET_SYNC_POINT(frame)     (GST_VIDEO_CODEC_FRAME_FLAG_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_SYNC_POINT))
#define GST_VIDEO_CODEC_FRAME_UNSET_SYNC_POINT(frame)   (GST_VIDEO_CODEC_FRAME_FLAG_UNSET(frame, GST_VIDEO_CODEC_FRAME_FLAG_SYNC_POINT))


/**
 * GST_VIDEO_CODEC_FRAME_IS_FORCE_KEYFRAME:
 * @frame: a #GstVideoCodecFrame
 *
 * Tests if the frame must be encoded as a keyframe. Applies only to
 * frames provided to encoders. Decoders can safely ignore this field.
 */
#define GST_VIDEO_CODEC_FRAME_IS_FORCE_KEYFRAME(frame)      (GST_VIDEO_CODEC_FRAME_FLAG_IS_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME))
#define GST_VIDEO_CODEC_FRAME_SET_FORCE_KEYFRAME(frame)     (GST_VIDEO_CODEC_FRAME_FLAG_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME))
#define GST_VIDEO_CODEC_FRAME_UNSET_FORCE_KEYFRAME(frame)   (GST_VIDEO_CODEC_FRAME_FLAG_UNSET(frame, GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME))

/**
 * GST_VIDEO_CODEC_FRAME_IS_FORCE_KEYFRAME_HEADERS:
 * @frame: a #GstVideoCodecFrame
 *
 * Tests if encoder should output stream headers before outputting the
 * resulting encoded buffer for the given frame.
 *
 * Applies only to frames provided to encoders. Decoders can safely
 * ignore this field.
 */
#define GST_VIDEO_CODEC_FRAME_IS_FORCE_KEYFRAME_HEADERS(frame)      (GST_VIDEO_CODEC_FRAME_FLAG_IS_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME_HEADERS))
#define GST_VIDEO_CODEC_FRAME_SET_FORCE_KEYFRAME_HEADERS(frame)     (GST_VIDEO_CODEC_FRAME_FLAG_SET(frame, GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME_HEADERS))
#define GST_VIDEO_CODEC_FRAME_UNSET_FORCE_KEYFRAME_HEADERS(frame)   (GST_VIDEO_CODEC_FRAME_FLAG_UNSET(frame, GST_VIDEO_CODEC_FRAME_FLAG_FORCE_KEYFRAME_HEADERS))

/**
 * GstVideoCodecFrame:
 * @pts: Presentation timestamp
 * @dts: Decoding timestamp
 * @duration: Duration of the frame
 * @system_frame_number: Unique identifier for the frame. Use this if you need
 *       to get hold of the frame later (like when data is being decoded).
 *       Typical usage in decoders is to set this on the opaque value provided
 *       to the library and get back the frame using gst_video_decoder_get_frame()
 * @distance_from_sync: Distance in frames from the last synchronization point.
 * @input_buffer: the input #GstBuffer that created this frame. The buffer is owned
 *           by the frame and references to the frame instead of the buffer should
 *           be kept.
 * @output_buffer: the output #GstBuffer. Implementations should set this either
 *           directly, or by using the
 *           gst_video_decoder_allocate_output_frame() or
 *           gst_video_decoder_allocate_output_buffer() methods. The buffer is
 *           owned by the frame and references to the frame instead of the
 *           buffer should be kept.
 * @deadline: Running time when the frame will be used.
 *
 * A #GstVideoCodecFrame represents a video frame both in raw and
 * encoded form.
 */
struct _GstVideoCodecFrame
{
  /*< private >*/
  gint ref_count;
  guint32 flags;

  /*< public >*/
  guint32 system_frame_number;	/* ED */

  /*< private >*/
  guint32 decode_frame_number;	/* ED */
  guint32 presentation_frame_number; /* ED */

  /*< public >*/
  GstClockTime dts;       /* ED */
  GstClockTime pts;       /* ED */
  GstClockTime duration;  /* ED */

  int distance_from_sync;	/* ED */

  GstBuffer *input_buffer;	/* ED */
  GstBuffer *output_buffer;	/* ED */

  GstClockTime deadline;	/* D */

  /*< private >*/

  /* Events that should be pushed downstream *before*
   * the next output_buffer */
  /* FIXME 2.0: Use a GQueue or similar */
  GList *events;		/* ED */

  gpointer       user_data;
  GDestroyNotify user_data_destroy_notify;

  union {
    struct {
      /*< private >*/
      GstClockTime ts;
      GstClockTime ts2;
      guint num_subframes;
      guint subframes_processed;
    } ABI;
    gpointer padding[GST_PADDING_LARGE];
  } abidata;
};

/* GstVideoCodecState */

GST_VIDEO_API
GType           gst_video_codec_state_get_type (void);

GST_VIDEO_API
GstVideoCodecState *gst_video_codec_state_ref (GstVideoCodecState * state);

GST_VIDEO_API
void                gst_video_codec_state_unref (GstVideoCodecState * state);


/* GstVideoCodecFrame */

GST_VIDEO_API
GType                gst_video_codec_frame_get_type (void);

GST_VIDEO_API
GstVideoCodecFrame  *gst_video_codec_frame_ref (GstVideoCodecFrame * frame);

GST_VIDEO_API
void                 gst_video_codec_frame_unref (GstVideoCodecFrame * frame);

GST_VIDEO_API
void                 gst_video_codec_frame_set_user_data (GstVideoCodecFrame *frame,
						          gpointer user_data,
				                          GDestroyNotify notify);

GST_VIDEO_API
gpointer             gst_video_codec_frame_get_user_data (GstVideoCodecFrame *frame);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstVideoCodecFrame, gst_video_codec_frame_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstVideoCodecState, gst_video_codec_state_unref)

G_END_DECLS

#endif
