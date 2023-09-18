/* GStreamer
 * Copyright (C) 2009 Igalia S.L.
 * Author: Iago Toral Quiroga <itoral@igalia.com>
 * Copyright (C) 2011 Mark Nauwelaerts <mark.nauwelaerts@collabora.co.uk>.
 * Copyright (C) 2011 Nokia Corporation. All rights reserved.
 *   Contact: Stefan Kost <stefan.kost@nokia.com>
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

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef _GST_AUDIO_DECODER_H_
#define _GST_AUDIO_DECODER_H_

#include <gst/gst.h>
#include <gst/base/gstadapter.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_DECODER \
  (gst_audio_decoder_get_type())
#define GST_AUDIO_DECODER(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_DECODER,GstAudioDecoder))
#define GST_AUDIO_DECODER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_DECODER,GstAudioDecoderClass))
#define GST_AUDIO_DECODER_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_AUDIO_DECODER,GstAudioDecoderClass))
#define GST_IS_AUDIO_DECODER(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_DECODER))
#define GST_IS_AUDIO_DECODER_CLASS(obj) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_DECODER))
#define GST_AUDIO_DECODER_CAST(obj) \
  ((GstAudioDecoder *)(obj))

/**
 * GST_AUDIO_DECODER_SINK_NAME:
 *
 * The name of the templates for the sink pad.
 */
#define GST_AUDIO_DECODER_SINK_NAME    "sink"
/**
 * GST_AUDIO_DECODER_SRC_NAME:
 *
 * The name of the templates for the source pad.
 */
#define GST_AUDIO_DECODER_SRC_NAME     "src"

/**
 * GST_AUDIO_DECODER_SRC_PAD:
 * @obj: base audio codec instance
 *
 * Gives the pointer to the source #GstPad object of the element.
 */
#define GST_AUDIO_DECODER_SRC_PAD(obj)         (((GstAudioDecoder *) (obj))->srcpad)

/**
 * GST_AUDIO_DECODER_SINK_PAD:
 * @obj: base audio codec instance
 *
 * Gives the pointer to the sink #GstPad object of the element.
 */
#define GST_AUDIO_DECODER_SINK_PAD(obj)        (((GstAudioDecoder *) (obj))->sinkpad)

#define GST_AUDIO_DECODER_STREAM_LOCK(dec)   g_rec_mutex_lock (&GST_AUDIO_DECODER (dec)->stream_lock)
#define GST_AUDIO_DECODER_STREAM_UNLOCK(dec) g_rec_mutex_unlock (&GST_AUDIO_DECODER (dec)->stream_lock)

/**
 * GST_AUDIO_DECODER_INPUT_SEGMENT:
 * @obj: audio decoder instance
 *
 * Gives the input segment of the element.
 */
#define GST_AUDIO_DECODER_INPUT_SEGMENT(obj)   (GST_AUDIO_DECODER_CAST (obj)->input_segment)

/**
 * GST_AUDIO_DECODER_OUTPUT_SEGMENT:
 * @obj: audio decoder instance
 *
 * Gives the output segment of the element.
 */
#define GST_AUDIO_DECODER_OUTPUT_SEGMENT(obj)   (GST_AUDIO_DECODER_CAST (obj)->output_segment)

typedef struct _GstAudioDecoder GstAudioDecoder;
typedef struct _GstAudioDecoderClass GstAudioDecoderClass;

typedef struct _GstAudioDecoderPrivate GstAudioDecoderPrivate;

/* do not use this one, use macro below */

GST_AUDIO_API
GstFlowReturn _gst_audio_decoder_error (GstAudioDecoder *dec, gint weight,
                                        GQuark domain, gint code,
                                        gchar *txt, gchar *debug,
                                        const gchar *file, const gchar *function,
                                        gint line);

/**
 * GST_AUDIO_DECODER_ERROR:
 * @el:     the base audio decoder element that generates the error
 * @weight: element defined weight of the error, added to error count
 * @domain: like CORE, LIBRARY, RESOURCE or STREAM (see #gstreamer-GstGError)
 * @code:   error code defined for that domain (see #gstreamer-GstGError)
 * @text:   the message to display (format string and args enclosed in
 *          parentheses)
 * @debug:  debugging information for the message (format string and args
 *          enclosed in parentheses)
 * @ret:    variable to receive return value
 *
 * Utility function that audio decoder elements can use in case they encountered
 * a data processing error that may be fatal for the current "data unit" but
 * need not prevent subsequent decoding.  Such errors are counted and if there
 * are too many, as configured in the context's max_errors, the pipeline will
 * post an error message and the application will be requested to stop further
 * media processing.  Otherwise, it is considered a "glitch" and only a warning
 * is logged. In either case, @ret is set to the proper value to
 * return to upstream/caller (indicating either GST_FLOW_ERROR or GST_FLOW_OK).
 */
#define GST_AUDIO_DECODER_ERROR(el, weight, domain, code, text, debug, ret) \
G_STMT_START {                                                              \
  gchar *__txt = _gst_element_error_printf text;                            \
  gchar *__dbg = _gst_element_error_printf debug;                           \
  GstAudioDecoder *__dec = GST_AUDIO_DECODER (el);                   \
  ret = _gst_audio_decoder_error (__dec, weight, GST_ ## domain ## _ERROR,    \
      GST_ ## domain ## _ERROR_ ## code, __txt, __dbg, __FILE__,            \
      GST_FUNCTION, __LINE__);                                              \
} G_STMT_END


/**
 * GST_AUDIO_DECODER_MAX_ERRORS:
 *
 * Default maximum number of errors tolerated before signaling error.
 */
#define GST_AUDIO_DECODER_MAX_ERRORS     -1

/**
 * GstAudioDecoder:
 *
 * The opaque #GstAudioDecoder data structure.
 */
struct _GstAudioDecoder
{
  GstElement element;

  /*< protected >*/
  /* source and sink pads */
  GstPad         *sinkpad;
  GstPad         *srcpad;

  /* protects all data processing, i.e. is locked
   * in the chain function, finish_frame and when
   * processing serialized events */
  GRecMutex       stream_lock;

  /* MT-protected (with STREAM_LOCK) */
  GstSegment      input_segment;
  GstSegment      output_segment;

  /*< private >*/
  GstAudioDecoderPrivate *priv;

  gpointer       _gst_reserved[GST_PADDING_LARGE];
};

/**
 * GstAudioDecoderClass:
 * @element_class:  The parent class structure
 * @start:          Optional.
 *                  Called when the element starts processing.
 *                  Allows opening external resources.
 * @stop:           Optional.
 *                  Called when the element stops processing.
 *                  Allows closing external resources.
 * @set_format:     Notifies subclass of incoming data format (caps).
 * @parse:          Optional.
 *                  Allows chopping incoming data into manageable units (frames)
 *                  for subsequent decoding.  This division is at subclass
 *                  discretion and may or may not correspond to 1 (or more)
 *                  frames as defined by audio format.
 * @handle_frame:   Provides input data (or NULL to clear any remaining data)
 *                  to subclass.  Input data ref management is performed by
 *                  base class, subclass should not care or intervene,
 *                  and input data is only valid until next call to base class,
 *                  most notably a call to gst_audio_decoder_finish_frame().
 * @flush:          Optional.
 *                  Instructs subclass to clear any codec caches and discard
 *                  any pending samples and not yet returned decoded data.
 *                  @hard indicates whether a FLUSH is being processed,
 *                  or otherwise a DISCONT (or conceptually similar).
 * @sink_event:     Optional.
 *                  Event handler on the sink pad. Subclasses should chain up to
 *                  the parent implementation to invoke the default handler.
 * @src_event:      Optional.
 *                  Event handler on the src pad. Subclasses should chain up to
 *                  the parent implementation to invoke the default handler.
 * @pre_push:       Optional.
 *                  Called just prior to pushing (encoded data) buffer downstream.
 *                  Subclass has full discretionary access to buffer,
 *                  and a not OK flow return will abort downstream pushing.
 * @open:           Optional.
 *                  Called when the element changes to GST_STATE_READY.
 *                  Allows opening external resources.
 * @close:          Optional.
 *                  Called when the element changes to GST_STATE_NULL.
 *                  Allows closing external resources.
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
 * @sink_query:     Optional.
 *                  Query handler on the sink pad. This function should
 *                  return TRUE if the query could be performed. Subclasses
 *                  should chain up to the parent implementation to invoke the
 *                  default handler. Since: 1.6
 * @src_query:      Optional.
 *                  Query handler on the source pad. This function should
 *                  return TRUE if the query could be performed. Subclasses
 *                  should chain up to the parent implementation to invoke the
 *                  default handler. Since: 1.6
 * @getcaps:        Optional.
 *                  Allows for a custom sink getcaps implementation.
 *                  If not implemented,
 *                  default returns gst_audio_decoder_proxy_getcaps
 *                  applied to sink template caps.
 * @transform_meta: Optional. Transform the metadata on the input buffer to the
 *                  output buffer. By default this method copies all meta without
 *                  tags and meta with only the "audio" tag. subclasses can
 *                  implement this method and return %TRUE if the metadata is to be
 *                  copied. Since: 1.6
 *
 * Subclasses can override any of the available virtual methods or not, as
 * needed. At minimum @handle_frame (and likely @set_format) needs to be
 * overridden.
 */
struct _GstAudioDecoderClass
{
  GstElementClass element_class;

  /*< public >*/
  /* virtual methods for subclasses */

  gboolean      (*start)              (GstAudioDecoder *dec);

  gboolean      (*stop)               (GstAudioDecoder *dec);

  gboolean      (*set_format)         (GstAudioDecoder *dec,
                                       GstCaps *caps);

  /**
   * GstAudioDecoderClass::parse:
   * @offset: (out):
   * @length: (out):
   */
  GstFlowReturn (*parse)              (GstAudioDecoder *dec,
                                       GstAdapter *adapter,
                                       gint *offset, gint *length);

  GstFlowReturn (*handle_frame)       (GstAudioDecoder *dec,
                                       GstBuffer *buffer);

  void          (*flush)              (GstAudioDecoder *dec, gboolean hard);

  GstFlowReturn (*pre_push)           (GstAudioDecoder *dec,
                                       GstBuffer **buffer);

  gboolean      (*sink_event)         (GstAudioDecoder *dec,
                                       GstEvent *event);
  gboolean      (*src_event)          (GstAudioDecoder *dec,
                                       GstEvent *event);

  gboolean      (*open)               (GstAudioDecoder *dec);

  gboolean      (*close)              (GstAudioDecoder *dec);

  gboolean      (*negotiate)          (GstAudioDecoder *dec);

  gboolean      (*decide_allocation)  (GstAudioDecoder *dec, GstQuery *query);

  gboolean      (*propose_allocation) (GstAudioDecoder *dec,
                                       GstQuery * query);

  gboolean      (*sink_query)         (GstAudioDecoder *dec, GstQuery *query);

  gboolean      (*src_query)          (GstAudioDecoder *dec, GstQuery *query);

  GstCaps *     (*getcaps)            (GstAudioDecoder * dec,
                                       GstCaps * filter);

  gboolean      (*transform_meta)     (GstAudioDecoder *enc, GstBuffer *outbuf,
                                       GstMeta *meta, GstBuffer *inbuf);

  /*< private >*/
  gpointer       _gst_reserved[GST_PADDING_LARGE - 4];
};

GST_AUDIO_API
GType             gst_audio_decoder_get_type (void);

GST_AUDIO_API
gboolean          gst_audio_decoder_set_output_format  (GstAudioDecoder    * dec,
                                                        const GstAudioInfo * info);

GST_AUDIO_API
gboolean          gst_audio_decoder_set_output_caps  (GstAudioDecoder * dec,
                                                      GstCaps         * caps);
GST_AUDIO_API
GstCaps *         gst_audio_decoder_proxy_getcaps (GstAudioDecoder * decoder,
                                                   GstCaps         * caps,
                                                   GstCaps         * filter);

GST_AUDIO_API
gboolean          gst_audio_decoder_negotiate (GstAudioDecoder * dec);

GST_AUDIO_API
GstFlowReturn     gst_audio_decoder_finish_subframe (GstAudioDecoder * dec,
                                                     GstBuffer       * buf);

GST_AUDIO_API
GstFlowReturn     gst_audio_decoder_finish_frame (GstAudioDecoder * dec,
                                                  GstBuffer * buf, gint frames);

GST_AUDIO_API
GstBuffer *       gst_audio_decoder_allocate_output_buffer (GstAudioDecoder * dec,
                                                            gsize              size);

/* context parameters */

GST_AUDIO_API
GstAudioInfo    * gst_audio_decoder_get_audio_info (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_plc_aware  (GstAudioDecoder * dec,
                                                    gboolean          plc);

GST_AUDIO_API
gint              gst_audio_decoder_get_plc_aware  (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_estimate_rate  (GstAudioDecoder * dec,
                                                    gboolean          enabled);

GST_AUDIO_API
gint              gst_audio_decoder_get_estimate_rate  (GstAudioDecoder * dec);

GST_AUDIO_API
gint              gst_audio_decoder_get_delay      (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_max_errors (GstAudioDecoder * dec,
                                                   gint               num);

GST_AUDIO_API
gint              gst_audio_decoder_get_max_errors (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_latency (GstAudioDecoder * dec,
                                                 GstClockTime      min,
                                                 GstClockTime      max);

GST_AUDIO_API
void              gst_audio_decoder_get_latency (GstAudioDecoder * dec,
                                                 GstClockTime    * min,
                                                 GstClockTime    * max);

GST_AUDIO_API
void              gst_audio_decoder_get_parse_state (GstAudioDecoder * dec,
                                                     gboolean        * sync,
                                                     gboolean        * eos);

GST_AUDIO_API
void              gst_audio_decoder_set_allocation_caps (GstAudioDecoder * dec,
                                                         GstCaps         * allocation_caps);

/* object properties */

GST_AUDIO_API
void              gst_audio_decoder_set_plc (GstAudioDecoder * dec,
                                             gboolean          enabled);

GST_AUDIO_API
gboolean          gst_audio_decoder_get_plc (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_min_latency (GstAudioDecoder * dec,
                                                     GstClockTime      num);

GST_AUDIO_API
GstClockTime      gst_audio_decoder_get_min_latency (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_tolerance   (GstAudioDecoder * dec,
                                                     GstClockTime      tolerance);

GST_AUDIO_API
GstClockTime      gst_audio_decoder_get_tolerance   (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_drainable (GstAudioDecoder * dec,
                                                   gboolean enabled);

GST_AUDIO_API
gboolean          gst_audio_decoder_get_drainable (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_set_needs_format (GstAudioDecoder * dec,
                                                      gboolean enabled);

GST_AUDIO_API
gboolean          gst_audio_decoder_get_needs_format (GstAudioDecoder * dec);

GST_AUDIO_API
void              gst_audio_decoder_get_allocator (GstAudioDecoder * dec,
                                                   GstAllocator ** allocator,
                                                   GstAllocationParams * params);

GST_AUDIO_API
void              gst_audio_decoder_merge_tags (GstAudioDecoder * dec,
                                                const GstTagList * tags, GstTagMergeMode mode);

GST_AUDIO_API
void              gst_audio_decoder_set_use_default_pad_acceptcaps (GstAudioDecoder * decoder,
                                                                   gboolean use);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioDecoder, gst_object_unref)

G_END_DECLS

#endif /* _GST_AUDIO_DECODER_H_ */
