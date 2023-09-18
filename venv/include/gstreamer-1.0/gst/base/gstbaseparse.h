/* GStreamer
 * Copyright (C) 2008 Nokia Corporation. All rights reserved.
 *
 * Contact: Stefan Kost <stefan.kost@nokia.com>
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

#ifndef __GST_BASE_PARSE_H__
#define __GST_BASE_PARSE_H__

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_BASE_PARSE            (gst_base_parse_get_type())
#define GST_BASE_PARSE(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_BASE_PARSE,GstBaseParse))
#define GST_BASE_PARSE_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_BASE_PARSE,GstBaseParseClass))
#define GST_BASE_PARSE_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_BASE_PARSE,GstBaseParseClass))
#define GST_IS_BASE_PARSE(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_BASE_PARSE))
#define GST_IS_BASE_PARSE_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_BASE_PARSE))
#define GST_BASE_PARSE_CAST(obj)       ((GstBaseParse *)(obj))

/**
 * GST_BASE_PARSE_SRC_PAD:
 * @obj: base parse instance
 *
 * Gives the pointer to the source #GstPad object of the element.
 */
#define GST_BASE_PARSE_SRC_PAD(obj)    (GST_BASE_PARSE_CAST (obj)->srcpad)

/**
 * GST_BASE_PARSE_SINK_PAD:
 * @obj: base parse instance
 *
 * Gives the pointer to the sink #GstPad object of the element.
 */
#define GST_BASE_PARSE_SINK_PAD(obj)    (GST_BASE_PARSE_CAST (obj)->sinkpad)

/**
 * GST_BASE_PARSE_FLOW_DROPPED:
 *
 * A #GstFlowReturn that can be returned from
 * #GstBaseParseClass::handle_frame to indicate that no output buffer was
 * generated, or from #GstBaseParseClass::pre_push_frame to to forego
 * pushing buffer.
 */
#define GST_BASE_PARSE_FLOW_DROPPED     GST_FLOW_CUSTOM_SUCCESS

/* not public API, use accessor macros below */
#define GST_BASE_PARSE_FLAG_LOST_SYNC (1 << 0)
#define GST_BASE_PARSE_FLAG_DRAINING  (1 << 1)

/**
 * GST_BASE_PARSE_LOST_SYNC:
 * @parse: base parse instance
 *
 * Obtains current sync status.
 */
#define GST_BASE_PARSE_LOST_SYNC(parse) (!!(GST_BASE_PARSE_CAST(parse)->flags & GST_BASE_PARSE_FLAG_LOST_SYNC))

/**
 * GST_BASE_PARSE_DRAINING:
 * @parse: base parse instance
 *
 * Obtains current drain status (ie. whether EOS has been received and
 * the parser is now processing the frames at the end of the stream)
 */
#define GST_BASE_PARSE_DRAINING(parse)  (!!(GST_BASE_PARSE_CAST(parse)->flags & GST_BASE_PARSE_FLAG_DRAINING))

/**
 * GstBaseParseFrameFlags:
 * @GST_BASE_PARSE_FRAME_FLAG_NONE: no flag
 * @GST_BASE_PARSE_FRAME_FLAG_NEW_FRAME: set by baseclass if current frame
 *   is passed for processing to the subclass for the first time
 *   (and not set on subsequent calls with same data).
 * @GST_BASE_PARSE_FRAME_FLAG_NO_FRAME: set to indicate this buffer should not be
 *   counted as frame, e.g. if this frame is dependent on a previous one.
 *   As it is not counted as a frame, bitrate increases but frame to time
 *   conversions are maintained.
 * @GST_BASE_PARSE_FRAME_FLAG_CLIP: @pre_push_frame can set this to indicate
 *    that regular segment clipping can still be performed (as opposed to
 *    any custom one having been done).
 * @GST_BASE_PARSE_FRAME_FLAG_DROP: indicates to @finish_frame that the
 *    the frame should be dropped (and might be handled internally by subclass)
 * @GST_BASE_PARSE_FRAME_FLAG_QUEUE: indicates to @finish_frame that the
 *    the frame should be queued for now and processed fully later
 *    when the first non-queued frame is finished
 *
 * Flags to be used in a #GstBaseParseFrame.
 */
typedef enum {
  GST_BASE_PARSE_FRAME_FLAG_NONE         = 0,
  GST_BASE_PARSE_FRAME_FLAG_NEW_FRAME    = (1 << 0),
  GST_BASE_PARSE_FRAME_FLAG_NO_FRAME     = (1 << 1),
  GST_BASE_PARSE_FRAME_FLAG_CLIP         = (1 << 2),
  GST_BASE_PARSE_FRAME_FLAG_DROP         = (1 << 3),
  GST_BASE_PARSE_FRAME_FLAG_QUEUE        = (1 << 4)
} GstBaseParseFrameFlags;

/**
 * GstBaseParseFrame:
 * @buffer: input data to be parsed for frames.
 * @out_buffer: output data.
 * @offset: media specific offset of input frame
 *   Note that a converter may have a different one on the frame's buffer.
 * @overhead: subclass can set this to indicates the metadata overhead
 *   for the given frame, which is then used to enable more accurate bitrate
 *   computations. If this is -1, it is assumed that this frame should be
 *   skipped in bitrate calculation.
 * @flags: a combination of input and output #GstBaseParseFrameFlags that
 *  convey additional context to subclass or allow subclass to tune
 *  subsequent #GstBaseParse actions.
 *
 * Frame (context) data passed to each frame parsing virtual methods.  In
 * addition to providing the data to be checked for a valid frame or an already
 * identified frame, it conveys additional metadata or control information
 * from and to the subclass w.r.t. the particular frame in question (rather
 * than global parameters).  Some of these may apply to each parsing stage, others
 * only to some a particular one.  These parameters are effectively zeroed at start
 * of each frame's processing, i.e. parsing virtual method invocation sequence.
 */
typedef struct {
  GstBuffer * buffer;
  GstBuffer * out_buffer;
  guint       flags;
  guint64     offset;
  gint        overhead;
  /*< private >*/
  gint        size;
  guint       _gst_reserved_i[2];
  gpointer    _gst_reserved_p[2];
  guint       _private_flags;
} GstBaseParseFrame;

typedef struct _GstBaseParse GstBaseParse;
typedef struct _GstBaseParseClass GstBaseParseClass;
typedef struct _GstBaseParsePrivate GstBaseParsePrivate;

/**
 * GstBaseParse:
 * @element: the parent element.
 *
 * The opaque #GstBaseParse data structure.
 */
struct _GstBaseParse {
  /*< public >*/
  GstElement     element;

  /*< protected >*/
  /* source and sink pads */
  GstPad         *sinkpad;
  GstPad         *srcpad;

  guint           flags;

  /* MT-protected (with STREAM_LOCK) */
  GstSegment      segment;

  /*< private >*/
  gpointer       _gst_reserved[GST_PADDING_LARGE];
  GstBaseParsePrivate *priv;
};

/**
 * GstBaseParseClass:
 * @parent_class: the parent class
 * @start:          Optional.
 *                  Called when the element starts processing.
 *                  Allows opening external resources.
 * @stop:           Optional.
 *                  Called when the element stops processing.
 *                  Allows closing external resources.
 * @set_sink_caps:  Optional.
 *                  Allows the subclass to be notified of the actual caps set.
 * @get_sink_caps:  Optional.
 *                  Allows the subclass to do its own sink get caps if needed.
 * @handle_frame:   Parses the input data into valid frames as defined by subclass
 *                  which should be passed to gst_base_parse_finish_frame().
 *                  The frame's input buffer is guaranteed writable,
 *                  whereas the input frame ownership is held by caller
 *                  (so subclass should make a copy if it needs to hang on).
 *                  Input buffer (data) is provided by baseclass with as much
 *                  metadata set as possible by baseclass according to upstream
 *                  information and/or subclass settings,
 *                  though subclass may still set buffer timestamp and duration
 *                  if desired.
 * @convert:        Optional.
 *                  Convert between formats.
 * @sink_event:     Optional.
 *                  Event handler on the sink pad. This function should chain
 *                  up to the parent implementation to let the default handler
 *                  run.
 * @src_event:      Optional.
 *                  Event handler on the source pad. Should chain up to the
 *                  parent to let the default handler run.
 * @pre_push_frame: Optional.
 *                   Called just prior to pushing a frame (after any pending
 *                   events have been sent) to give subclass a chance to perform
 *                   additional actions at this time (e.g. tag sending) or to
 *                   decide whether this buffer should be dropped or not
 *                   (e.g. custom segment clipping).
 * @detect:         Optional.
 *                   Called until it doesn't return GST_FLOW_OK anymore for
 *                   the first buffers. Can be used by the subclass to detect
 *                   the stream format.
 * @sink_query:     Optional.
 *                   Query handler on the sink pad. This function should chain
 *                   up to the parent implementation to let the default handler
 *                   run (Since: 1.2)
 * @src_query:      Optional.
 *                   Query handler on the source pad. Should chain up to the
 *                   parent to let the default handler run (Since: 1.2)
 *
 * Subclasses can override any of the available virtual methods or not, as
 * needed. At minimum @handle_frame needs to be overridden.
 */
struct _GstBaseParseClass {
  GstElementClass parent_class;

  /*< public >*/
  /* virtual methods for subclasses */

  gboolean      (*start)              (GstBaseParse * parse);

  gboolean      (*stop)               (GstBaseParse * parse);

  gboolean      (*set_sink_caps)      (GstBaseParse * parse,
                                       GstCaps      * caps);

  /**
   * GstBaseParseClass::handle_frame:
   * @skipsize: (out):
   *
   * Parses the input data into valid frames as defined by subclass
   * which should be passed to gst_base_parse_finish_frame().
   * The frame's input buffer is guaranteed writable,
   * whereas the input frame ownership is held by caller
   * (so subclass should make a copy if it needs to hang on).
   * Input buffer (data) is provided by baseclass with as much
   * metadata set as possible by baseclass according to upstream
   * information and/or subclass settings,
   * though subclass may still set buffer timestamp and duration
   * if desired.
   */
  GstFlowReturn (*handle_frame)       (GstBaseParse      * parse,
                                       GstBaseParseFrame * frame,
                                       gint              * skipsize);

  GstFlowReturn (*pre_push_frame)     (GstBaseParse      * parse,
                                       GstBaseParseFrame * frame);

  gboolean      (*convert)            (GstBaseParse * parse,
                                       GstFormat      src_format,
                                       gint64         src_value,
                                       GstFormat      dest_format,
                                       gint64       * dest_value);

  gboolean      (*sink_event)         (GstBaseParse * parse,
                                       GstEvent     * event);

  gboolean      (*src_event)          (GstBaseParse * parse,
                                       GstEvent     * event);

  GstCaps *     (*get_sink_caps)      (GstBaseParse * parse,
                                       GstCaps      * filter);

  GstFlowReturn (*detect)             (GstBaseParse * parse,
                                       GstBuffer    * buffer);

  gboolean      (*sink_query)         (GstBaseParse * parse,
                                       GstQuery     * query);

  gboolean      (*src_query)          (GstBaseParse * parse,
                                       GstQuery     * query);

  /*< private >*/
  gpointer       _gst_reserved[GST_PADDING_LARGE - 2];
};

GST_BASE_API
GType           gst_base_parse_get_type (void);

GST_BASE_API
GType           gst_base_parse_frame_get_type (void);

GST_BASE_API
GstBaseParseFrame * gst_base_parse_frame_new  (GstBuffer              * buffer,
                                               GstBaseParseFrameFlags   flags,
                                               gint                     overhead);
GST_BASE_API
void            gst_base_parse_frame_init      (GstBaseParseFrame * frame);

GST_BASE_API
GstBaseParseFrame * gst_base_parse_frame_copy  (GstBaseParseFrame * frame);
GST_BASE_API
void            gst_base_parse_frame_free      (GstBaseParseFrame * frame);

GST_BASE_API
GstFlowReturn   gst_base_parse_push_frame      (GstBaseParse      * parse,
                                                GstBaseParseFrame * frame);
GST_BASE_API
GstFlowReturn   gst_base_parse_finish_frame    (GstBaseParse * parse,
                                                GstBaseParseFrame * frame,
                                                gint size);
GST_BASE_API
void            gst_base_parse_set_duration    (GstBaseParse      * parse,
                                                GstFormat           fmt,
                                                gint64              duration,
                                                gint                interval);
GST_BASE_API
void            gst_base_parse_set_average_bitrate (GstBaseParse   * parse,
                                                    guint            bitrate);
GST_BASE_API
void            gst_base_parse_set_min_frame_size (GstBaseParse    * parse,
                                                   guint             min_size);
GST_BASE_API
void            gst_base_parse_set_has_timing_info (GstBaseParse   * parse,
                                                    gboolean         has_timing);
GST_BASE_API
void            gst_base_parse_drain           (GstBaseParse * parse);

GST_BASE_API
void            gst_base_parse_set_syncable    (GstBaseParse * parse,
                                                gboolean       syncable);
GST_BASE_API
void            gst_base_parse_set_passthrough (GstBaseParse * parse,
                                                gboolean       passthrough);
GST_BASE_API
void            gst_base_parse_set_pts_interpolation (GstBaseParse * parse,
                                                      gboolean pts_interpolate);
GST_BASE_API
void            gst_base_parse_set_infer_ts (GstBaseParse * parse,
                                             gboolean infer_ts);
GST_BASE_API
void            gst_base_parse_set_frame_rate  (GstBaseParse * parse,
                                                guint          fps_num,
                                                guint          fps_den,
                                                guint          lead_in,
                                                guint          lead_out);
GST_BASE_API
void            gst_base_parse_set_latency     (GstBaseParse * parse,
                                                GstClockTime min_latency,
                                                GstClockTime max_latency);
GST_BASE_API
gboolean        gst_base_parse_convert_default (GstBaseParse * parse,
                                                GstFormat      src_format,
                                                gint64         src_value,
                                                GstFormat      dest_format,
                                                gint64       * dest_value);
GST_BASE_API
gboolean        gst_base_parse_add_index_entry (GstBaseParse * parse,
                                                guint64        offset,
                                                GstClockTime   ts,
                                                gboolean       key,
                                                gboolean       force);
GST_BASE_API
void            gst_base_parse_set_ts_at_offset (GstBaseParse *parse,
                                                 gsize offset);
GST_BASE_API
void            gst_base_parse_merge_tags       (GstBaseParse  * parse,
                                                 GstTagList    * tags,
                                                 GstTagMergeMode mode);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBaseParseFrame, gst_base_parse_frame_free)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBaseParse, gst_object_unref)

G_END_DECLS

#endif /* __GST_BASE_PARSE_H__ */
