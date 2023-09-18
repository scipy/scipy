/* GStreamer
 * Copyright (C) 2008 David Schleef <ds@schleef.org>
 * Copyright (C) 2011 Mark Nauwelaerts <mark.nauwelaerts@collabora.co.uk>.
 * Copyright (C) 2011 Nokia Corporation. All rights reserved.
 *   Contact: Stefan Kost <stefan.kost@nokia.com>
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

#ifndef _GST_VIDEO_DECODER_H_
#define _GST_VIDEO_DECODER_H_

#include <gst/base/gstadapter.h>
#include <gst/video/gstvideoutils.h>

G_BEGIN_DECLS

#define GST_TYPE_VIDEO_DECODER \
  (gst_video_decoder_get_type())
#define GST_VIDEO_DECODER(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_VIDEO_DECODER,GstVideoDecoder))
#define GST_VIDEO_DECODER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_VIDEO_DECODER,GstVideoDecoderClass))
#define GST_VIDEO_DECODER_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_VIDEO_DECODER,GstVideoDecoderClass))
#define GST_IS_VIDEO_DECODER(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_VIDEO_DECODER))
#define GST_IS_VIDEO_DECODER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_VIDEO_DECODER))
#define GST_VIDEO_DECODER_CAST(obj) ((GstVideoDecoder *)(obj))

/**
 * GST_VIDEO_DECODER_SINK_NAME:
 *
 * The name of the templates for the sink pad.
 */
#define GST_VIDEO_DECODER_SINK_NAME    "sink"
/**
 * GST_VIDEO_DECODER_SRC_NAME:
 *
 * The name of the templates for the source pad.
 */
#define GST_VIDEO_DECODER_SRC_NAME     "src"

/**
 * GST_VIDEO_DECODER_SRC_PAD:
 * @obj: a #GstVideoDecoder
 *
 * Gives the pointer to the source #GstPad object of the element.
 */
#define GST_VIDEO_DECODER_SRC_PAD(obj)         (((GstVideoDecoder *) (obj))->srcpad)

/**
 * GST_VIDEO_DECODER_SINK_PAD:
 * @obj: a #GstVideoDecoder
 *
 * Gives the pointer to the sink #GstPad object of the element.
 */
#define GST_VIDEO_DECODER_SINK_PAD(obj)        (((GstVideoDecoder *) (obj))->sinkpad)
/**
 * GST_VIDEO_DECODER_FLOW_NEED_DATA:
 *
 * Returned while parsing to indicate more data is needed.
 **/
#define GST_VIDEO_DECODER_FLOW_NEED_DATA GST_FLOW_CUSTOM_SUCCESS

/**
 * GST_VIDEO_DECODER_INPUT_SEGMENT:
 * @obj: base decoder instance
 *
 * Gives the segment of the element.
 */
#define GST_VIDEO_DECODER_INPUT_SEGMENT(obj)     (GST_VIDEO_DECODER_CAST (obj)->input_segment)

/**
 * GST_VIDEO_DECODER_OUTPUT_SEGMENT:
 * @obj: base decoder instance
 *
 * Gives the segment of the element.
 */
#define GST_VIDEO_DECODER_OUTPUT_SEGMENT(obj)     (GST_VIDEO_DECODER_CAST (obj)->output_segment)

/**
 * GST_VIDEO_DECODER_STREAM_LOCK:
 * @decoder: video decoder instance
 *
 * Obtain a lock to protect the decoder function from concurrent access.
 */
#define GST_VIDEO_DECODER_STREAM_LOCK(decoder) g_rec_mutex_lock (&GST_VIDEO_DECODER (decoder)->stream_lock)

/**
 * GST_VIDEO_DECODER_STREAM_UNLOCK:
 * @decoder: video decoder instance
 *
 * Release the lock that protects the decoder function from concurrent access.
 */
#define GST_VIDEO_DECODER_STREAM_UNLOCK(decoder) g_rec_mutex_unlock (&GST_VIDEO_DECODER (decoder)->stream_lock)

typedef struct _GstVideoDecoder GstVideoDecoder;
typedef struct _GstVideoDecoderClass GstVideoDecoderClass;
typedef struct _GstVideoDecoderPrivate GstVideoDecoderPrivate;


/* do not use this one, use macro below */

GST_VIDEO_API
GstFlowReturn _gst_video_decoder_error (GstVideoDecoder *dec, gint weight,
                                             GQuark domain, gint code,
                                             gchar *txt, gchar *debug,
                                             const gchar *file, const gchar *function,
                                             gint line);

/**
 * GST_VIDEO_DECODER_ERROR:
 * @el:     the base video decoder element that generates the error
 * @w:      element defined weight of the error, added to error count
 * @domain: like CORE, LIBRARY, RESOURCE or STREAM (see #gstreamer-GstGError)
 * @code:   error code defined for that domain (see #gstreamer-GstGError)
 * @text:   the message to display (format string and args enclosed in
 *          parentheses)
 * @debug:  debugging information for the message (format string and args
 *          enclosed in parentheses)
 * @ret:    variable to receive return value
 *
 * Utility function that video decoder elements can use in case they encountered
 * a data processing error that may be fatal for the current "data unit" but
 * need not prevent subsequent decoding.  Such errors are counted and if there
 * are too many, as configured in the context's max_errors, the pipeline will
 * post an error message and the application will be requested to stop further
 * media processing.  Otherwise, it is considered a "glitch" and only a warning
 * is logged. In either case, @ret is set to the proper value to
 * return to upstream/caller (indicating either GST_FLOW_ERROR or GST_FLOW_OK).
 */
#define GST_VIDEO_DECODER_ERROR(el, w, domain, code, text, debug, ret) \
G_STMT_START {                                                              \
  gchar *__txt = _gst_element_error_printf text;                            \
  gchar *__dbg = _gst_element_error_printf debug;                           \
  GstVideoDecoder *__dec = GST_VIDEO_DECODER (el);                   \
  ret = _gst_video_decoder_error (__dec, w, GST_ ## domain ## _ERROR,    \
      GST_ ## domain ## _ERROR_ ## code, __txt, __dbg, __FILE__,            \
      GST_FUNCTION, __LINE__);                                              \
} G_STMT_END

/**
 * GST_VIDEO_DECODER_MAX_ERRORS:
 *
 * Default maximum number of errors tolerated before signaling error.
 */
#define GST_VIDEO_DECODER_MAX_ERRORS     -1


/**
 * GstVideoDecoder:
 *
 * The opaque #GstVideoDecoder data structure.
 */
struct _GstVideoDecoder
{
  /*< private >*/
  GstElement     element;

  /*< protected >*/
  GstPad         *sinkpad;
  GstPad         *srcpad;

  /* protects all data processing, i.e. is locked
   * in the chain function, finish_frame and when
   * processing serialized events */
  GRecMutex stream_lock;

  /* MT-protected (with STREAM_LOCK) */
  GstSegment      input_segment;
  GstSegment      output_segment;

  GstVideoDecoderPrivate *priv;

  /*< private >*/
  gpointer padding[GST_PADDING_LARGE];
};

/**
 * GstVideoDecoderClass:
 * @open:           Optional.
 *                  Called when the element changes to GST_STATE_READY.
 *                  Allows opening external resources.
 * @close:          Optional.
 *                  Called when the element changes to GST_STATE_NULL.
 *                  Allows closing external resources.
 * @start:          Optional.
 *                  Called when the element starts processing.
 *                  Allows opening external resources.
 * @stop:           Optional.
 *                  Called when the element stops processing.
 *                  Allows closing external resources.
 * @set_format:     Notifies subclass of incoming data format (caps).
 * @parse:          Required for non-packetized input.
 *                  Allows chopping incoming data into manageable units (frames)
 *                  for subsequent decoding.
 * @reset:          Optional.
 *                  Allows subclass (decoder) to perform post-seek semantics reset.
 *                  Deprecated.
 * @handle_frame:   Provides input data frame to subclass. In subframe mode, the subclass needs
 *                  to take ownership of @GstVideoCodecFrame.input_buffer as it will be modified
 *                  by the base class on the next subframe buffer receiving.
 * @finish:         Optional.
 *                  Called to request subclass to dispatch any pending remaining
 *                  data at EOS. Sub-classes can refuse to decode new data after.
 * @drain:	    Optional.
 *                  Called to request subclass to decode any data it can at this
 *                  point, but that more data may arrive after. (e.g. at segment end).
 *                  Sub-classes should be prepared to handle new data afterward,
 *                  or seamless segment processing will break. Since: 1.6
 * @sink_event:     Optional.
 *                  Event handler on the sink pad. This function should return
 *                  TRUE if the event was handled and should be discarded
 *                  (i.e. not unref'ed).
 *                  Subclasses should chain up to the parent implementation to
 *                  invoke the default handler.
 * @src_event:      Optional.
 *                  Event handler on the source pad. This function should return
 *                  TRUE if the event was handled and should be discarded
 *                  (i.e. not unref'ed).
 *                  Subclasses should chain up to the parent implementation to
 *                  invoke the default handler.
 * @negotiate:      Optional.
 *                  Negotiate with downstream and configure buffer pools, etc.
 *                  Subclasses should chain up to the parent implementation to
 *                  invoke the default handler.
 * @decide_allocation: Optional.
 *                     Setup the allocation parameters for allocating output
 *                     buffers. The passed in query contains the result of the
 *                     downstream allocation query.
 *                     Subclasses should chain up to the parent implementation to
 *                     invoke the default handler.
 * @propose_allocation: Optional.
 *                      Propose buffer allocation parameters for upstream elements.
 *                      Subclasses should chain up to the parent implementation to
 *                      invoke the default handler.
 * @flush:              Optional.
 *                      Flush all remaining data from the decoder without
 *                      pushing it downstream. Since: 1.2
 * @sink_query:     Optional.
 *                  Query handler on the sink pad. This function should
 *                  return TRUE if the query could be performed. Subclasses
 *                  should chain up to the parent implementation to invoke the
 *                  default handler. Since: 1.4
 * @src_query:      Optional.
 *                  Query handler on the source pad. This function should
 *                  return TRUE if the query could be performed. Subclasses
 *                  should chain up to the parent implementation to invoke the
 *                  default handler. Since: 1.4
 * @getcaps:        Optional.
 *                  Allows for a custom sink getcaps implementation.
 *                  If not implemented, default returns
 *                  gst_video_decoder_proxy_getcaps
 *                  applied to sink template caps.
 * @transform_meta: Optional. Transform the metadata on the input buffer to the
 *                  output buffer. By default this method is copies all meta without
 *                  tags and meta with only the "video" tag. subclasses can
 *                  implement this method and return %TRUE if the metadata is to be
 *                  copied. Since: 1.6
 *
 * Subclasses can override any of the available virtual methods or not, as
 * needed. At minimum @handle_frame needs to be overridden, and @set_format
 * and likely as well.  If non-packetized input is supported or expected,
 * @parse needs to be overridden as well.
 */
struct _GstVideoDecoderClass
{
  /*< private >*/
  GstElementClass element_class;

  /*< public >*/
  gboolean      (*open)           (GstVideoDecoder *decoder);

  gboolean      (*close)          (GstVideoDecoder *decoder);

  gboolean      (*start)          (GstVideoDecoder *decoder);

  gboolean      (*stop)           (GstVideoDecoder *decoder);

  GstFlowReturn (*parse)          (GstVideoDecoder *decoder,
				   GstVideoCodecFrame *frame,
				   GstAdapter *adapter,
				   gboolean at_eos);

  gboolean      (*set_format)     (GstVideoDecoder *decoder,
				   GstVideoCodecState * state);

  gboolean      (*reset)          (GstVideoDecoder *decoder,
				   gboolean hard);

  GstFlowReturn (*finish)         (GstVideoDecoder *decoder);

  /**
   * GstVideoDecoderClass::handle_frame:
   * @decoder: The #GstVideoDecoder
   * @frame: (transfer full): The frame to handle
   */
  GstFlowReturn (*handle_frame)   (GstVideoDecoder *decoder,
				   GstVideoCodecFrame *frame);

  gboolean      (*sink_event)     (GstVideoDecoder *decoder,
				   GstEvent *event);

  gboolean      (*src_event)      (GstVideoDecoder *decoder,
				   GstEvent *event);

  gboolean      (*negotiate)      (GstVideoDecoder *decoder);

  gboolean      (*decide_allocation)  (GstVideoDecoder *decoder, GstQuery *query);

  gboolean      (*propose_allocation) (GstVideoDecoder *decoder, GstQuery * query);

  gboolean      (*flush)              (GstVideoDecoder *decoder);

  gboolean      (*sink_query)     (GstVideoDecoder *decoder,
				   GstQuery *query);

  gboolean      (*src_query)      (GstVideoDecoder *decoder,
				   GstQuery *query);

  GstCaps*      (*getcaps)        (GstVideoDecoder *decoder,
                                   GstCaps *filter);

  GstFlowReturn (*drain)          (GstVideoDecoder *decoder);

  gboolean      (*transform_meta) (GstVideoDecoder *decoder,
                                   GstVideoCodecFrame *frame,
                                   GstMeta * meta);

  /**
   * GstVideoDecoderClass::handle_missing_data:
   * @decoder: The #GstVideoDecoder
   * @timestamp: Timestamp of the missing data
   * @duration: Duration of the missing data
   *
   * Returns: %TRUE if the decoder should be drained afterwards.
   *
   * Since: 1.20
   */
  gboolean      (*handle_missing_data) (GstVideoDecoder *decoder,
                                        GstClockTime timestamp,
                                        GstClockTime duration);

  /*< private >*/
  gpointer padding[GST_PADDING_LARGE-7];
};

/**
 * GstVideoDecoderRequestSyncPointFlags:
 * @GST_VIDEO_DECODER_REQUEST_SYNC_POINT_DISCARD_INPUT: discard all following
 *     input until the next sync point.
 * @GST_VIDEO_DECODER_REQUEST_SYNC_POINT_CORRUPT_OUTPUT: discard all following
 *     output until the next sync point.
 *
 * Flags to be used in combination with gst_video_decoder_request_sync_point().
 * See the function documentation for more details.
 *
 * Since: 1.20
 */
typedef enum {
  GST_VIDEO_DECODER_REQUEST_SYNC_POINT_DISCARD_INPUT  = (1<<0),
  GST_VIDEO_DECODER_REQUEST_SYNC_POINT_CORRUPT_OUTPUT = (1<<1),
} GstVideoDecoderRequestSyncPointFlags;

GST_VIDEO_API
GType    gst_video_decoder_get_type (void);

/* Context parameters */

GST_VIDEO_API
void     gst_video_decoder_set_packetized (GstVideoDecoder * decoder,
					   gboolean packetized);

GST_VIDEO_API
gboolean gst_video_decoder_get_packetized (GstVideoDecoder * decoder);

GST_VIDEO_API
void     gst_video_decoder_set_subframe_mode (GstVideoDecoder * decoder,
                                              gboolean subframe_mode);

GST_VIDEO_API
gboolean gst_video_decoder_get_subframe_mode (GstVideoDecoder * decoder);

GST_VIDEO_API
guint gst_video_decoder_get_input_subframe_index (GstVideoDecoder * decoder, GstVideoCodecFrame * frame);

GST_VIDEO_API
guint gst_video_decoder_get_processed_subframe_index (GstVideoDecoder * decoder, GstVideoCodecFrame * frame);

GST_VIDEO_API
void     gst_video_decoder_set_estimate_rate (GstVideoDecoder * dec,
					      gboolean          enabled);

GST_VIDEO_API
gint     gst_video_decoder_get_estimate_rate (GstVideoDecoder * dec);

GST_VIDEO_API
void     gst_video_decoder_set_max_errors (GstVideoDecoder * dec,
					   gint              num);

GST_VIDEO_API
gint     gst_video_decoder_get_max_errors (GstVideoDecoder * dec);

GST_VIDEO_API
void     gst_video_decoder_set_needs_format (GstVideoDecoder * dec,
                                             gboolean enabled);

GST_VIDEO_API
gboolean gst_video_decoder_get_needs_format (GstVideoDecoder * dec);

GST_VIDEO_API
void     gst_video_decoder_set_needs_sync_point (GstVideoDecoder * dec,
                                                 gboolean enabled);

GST_VIDEO_API
gboolean gst_video_decoder_get_needs_sync_point (GstVideoDecoder * dec);

GST_VIDEO_API
void     gst_video_decoder_set_latency (GstVideoDecoder *decoder,
					GstClockTime min_latency,
					GstClockTime max_latency);

GST_VIDEO_API
void     gst_video_decoder_get_latency (GstVideoDecoder *decoder,
					GstClockTime *min_latency,
					GstClockTime *max_latency);

GST_VIDEO_API
void     gst_video_decoder_get_allocator (GstVideoDecoder *decoder,
                                          GstAllocator **allocator,
                                          GstAllocationParams *params);

GST_VIDEO_API
GstBufferPool *gst_video_decoder_get_buffer_pool (GstVideoDecoder *decoder);

/* Object methods */

GST_VIDEO_API
GstVideoCodecFrame *gst_video_decoder_get_frame        (GstVideoDecoder *decoder,
						        int frame_number);

GST_VIDEO_API
GstVideoCodecFrame *gst_video_decoder_get_oldest_frame (GstVideoDecoder *decoder);

GST_VIDEO_API
GList *             gst_video_decoder_get_frames       (GstVideoDecoder *decoder);

/* Parsing related methods */

GST_VIDEO_API
void           gst_video_decoder_add_to_frame     (GstVideoDecoder *decoder,
						   int n_bytes);

GST_VIDEO_API
GstFlowReturn  gst_video_decoder_have_frame       (GstVideoDecoder *decoder);

GST_VIDEO_API
GstFlowReturn  gst_video_decoder_have_last_subframe (GstVideoDecoder *decoder,
                                                     GstVideoCodecFrame * frame);

GST_VIDEO_API
gsize          gst_video_decoder_get_pending_frame_size (GstVideoDecoder *decoder);

GST_VIDEO_API
GstBuffer     *gst_video_decoder_allocate_output_buffer (GstVideoDecoder * decoder);

GST_VIDEO_API
GstFlowReturn  gst_video_decoder_allocate_output_frame_with_params (GstVideoDecoder *decoder,
                                                                    GstVideoCodecFrame * frame,
                                                                    GstBufferPoolAcquireParams *params);

GST_VIDEO_API
GstFlowReturn  gst_video_decoder_allocate_output_frame  (GstVideoDecoder *decoder,
						         GstVideoCodecFrame *frame);

GST_VIDEO_API
GstVideoCodecState *gst_video_decoder_set_output_state (GstVideoDecoder *decoder,
							GstVideoFormat fmt, guint width, guint height,
							GstVideoCodecState *reference);

GST_VIDEO_API
GstVideoCodecState *gst_video_decoder_set_interlaced_output_state (GstVideoDecoder *decoder,
                                                                   GstVideoFormat fmt, GstVideoInterlaceMode interlace_mode,
                                                                   guint width, guint height, GstVideoCodecState *reference);

GST_VIDEO_API
GstVideoCodecState *gst_video_decoder_get_output_state (GstVideoDecoder *decoder);

GST_VIDEO_API
gboolean         gst_video_decoder_negotiate           (GstVideoDecoder * decoder);

GST_VIDEO_API
GstClockTimeDiff gst_video_decoder_get_max_decode_time (GstVideoDecoder *decoder,
							GstVideoCodecFrame *frame);

GST_VIDEO_API
gdouble          gst_video_decoder_get_qos_proportion (GstVideoDecoder * decoder);

GST_VIDEO_API
GstFlowReturn    gst_video_decoder_finish_frame (GstVideoDecoder *decoder,
						 GstVideoCodecFrame *frame);
GST_VIDEO_API
GstFlowReturn    gst_video_decoder_finish_subframe (GstVideoDecoder *decoder,
                                                 GstVideoCodecFrame *frame);

GST_VIDEO_API
GstFlowReturn    gst_video_decoder_drop_frame (GstVideoDecoder *dec,
					       GstVideoCodecFrame *frame);
GST_VIDEO_API
GstFlowReturn    gst_video_decoder_drop_subframe (GstVideoDecoder *dec,
                                               GstVideoCodecFrame *frame);

GST_VIDEO_API
void             gst_video_decoder_request_sync_point (GstVideoDecoder *dec,
                                                       GstVideoCodecFrame *frame,
                                                       GstVideoDecoderRequestSyncPointFlags flags);

GST_VIDEO_API
void             gst_video_decoder_release_frame (GstVideoDecoder * dec,
						  GstVideoCodecFrame * frame);

GST_VIDEO_API
void             gst_video_decoder_merge_tags (GstVideoDecoder *decoder,
                                               const GstTagList *tags,
                                               GstTagMergeMode mode);

GST_VIDEO_API
GstCaps *        gst_video_decoder_proxy_getcaps (GstVideoDecoder * decoder,
						  GstCaps         * caps,
                                                  GstCaps         * filter);

GST_VIDEO_API
void             gst_video_decoder_set_use_default_pad_acceptcaps (GstVideoDecoder * decoder,
                                                                   gboolean use);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstVideoDecoder, gst_object_unref)

G_END_DECLS

#endif

