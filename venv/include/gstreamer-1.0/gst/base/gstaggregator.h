/* GStreamer aggregator base class
 * Copyright (C) 2014 Mathieu Duponchelle <mathieu.duponchelle@oencreed.com>
 * Copyright (C) 2014 Thibault Saunier <tsaunier@gnome.org>
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

#ifndef __GST_AGGREGATOR_H__
#define __GST_AGGREGATOR_H__

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

G_BEGIN_DECLS

/**************************
 * GstAggregator Structs  *
 *************************/

typedef struct _GstAggregator GstAggregator;
typedef struct _GstAggregatorPrivate GstAggregatorPrivate;
typedef struct _GstAggregatorClass GstAggregatorClass;

/************************
 * GstAggregatorPad API *
 ***********************/

#define GST_TYPE_AGGREGATOR_PAD            (gst_aggregator_pad_get_type())
#define GST_AGGREGATOR_PAD(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AGGREGATOR_PAD, GstAggregatorPad))
#define GST_AGGREGATOR_PAD_CAST(obj)       ((GstAggregatorPad *)(obj))
#define GST_AGGREGATOR_PAD_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AGGREGATOR_PAD, GstAggregatorPadClass))
#define GST_AGGREGATOR_PAD_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AGGREGATOR_PAD, GstAggregatorPadClass))
#define GST_IS_AGGREGATOR_PAD(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AGGREGATOR_PAD))
#define GST_IS_AGGREGATOR_PAD_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AGGREGATOR_PAD))

/****************************
 * GstAggregatorPad Structs *
 ***************************/

typedef struct _GstAggregatorPad GstAggregatorPad;
typedef struct _GstAggregatorPadClass GstAggregatorPadClass;
typedef struct _GstAggregatorPadPrivate GstAggregatorPadPrivate;

/**
 * GstAggregatorPad:
 * @segment: last segment received.
 *
 * The implementation the GstPad to use with #GstAggregator
 *
 * Since: 1.14
 */
struct _GstAggregatorPad
{
  GstPad                       parent;

  /*< public >*/
  /* Protected by the OBJECT_LOCK */
  GstSegment segment;

  /* < private > */
  GstAggregatorPadPrivate   *  priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAggregatorPadClass:
 * @flush:       Optional
 *               Called when the pad has received a flush stop, this is the place
 *               to flush any information specific to the pad, it allows for individual
 *               pads to be flushed while others might not be.
 * @skip_buffer: Optional
 *               Called before input buffers are queued in the pad, return %TRUE
 *               if the buffer should be skipped.
 *
 * Since: 1.14
 */
struct _GstAggregatorPadClass
{
  GstPadClass   parent_class;

  GstFlowReturn (*flush)       (GstAggregatorPad * aggpad, GstAggregator * aggregator);
  gboolean      (*skip_buffer) (GstAggregatorPad * aggpad, GstAggregator * aggregator, GstBuffer * buffer);

  /*< private >*/
  gpointer      _gst_reserved[GST_PADDING_LARGE];
};

GST_BASE_API
GType gst_aggregator_pad_get_type           (void);

/****************************
 * GstAggregatorPad methods *
 ***************************/

GST_BASE_API
GstBuffer * gst_aggregator_pad_pop_buffer   (GstAggregatorPad *  pad);

GST_BASE_API
GstBuffer * gst_aggregator_pad_peek_buffer  (GstAggregatorPad *  pad);

GST_BASE_API
gboolean    gst_aggregator_pad_drop_buffer  (GstAggregatorPad *  pad);

GST_BASE_API
gboolean    gst_aggregator_pad_has_buffer   (GstAggregatorPad * pad);

GST_BASE_API
gboolean    gst_aggregator_pad_is_eos       (GstAggregatorPad *  pad);

GST_BASE_API
gboolean    gst_aggregator_pad_is_inactive  (GstAggregatorPad * pad);

/*********************
 * GstAggregator API *
 ********************/

#define GST_TYPE_AGGREGATOR            (gst_aggregator_get_type())
#define GST_AGGREGATOR(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AGGREGATOR,GstAggregator))
#define GST_AGGREGATOR_CAST(obj)       ((GstAggregator *)(obj))
#define GST_AGGREGATOR_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AGGREGATOR,GstAggregatorClass))
#define GST_AGGREGATOR_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AGGREGATOR,GstAggregatorClass))
#define GST_IS_AGGREGATOR(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AGGREGATOR))
#define GST_IS_AGGREGATOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AGGREGATOR))

#define GST_AGGREGATOR_FLOW_NEED_DATA             GST_FLOW_CUSTOM_ERROR

/**
 * GstAggregator:
 * @srcpad: the aggregator's source pad
 *
 * Aggregator base class object structure.
 *
 * Since: 1.14
 */
struct _GstAggregator
{
  GstElement               parent;

  /*< public >*/
  GstPad                *  srcpad;

  /*< private >*/
  GstAggregatorPrivate  *  priv;

  gpointer                 _gst_reserved[GST_PADDING_LARGE];
};

/**
 * GstAggregatorClass:
 * @flush:          Optional.
 *                  Called after a successful flushing seek, once all the flush
 *                  stops have been received. Flush pad-specific data in
 *                  #GstAggregatorPad->flush.
 * @clip:           Optional.
 *                  Called when a buffer is received on a sink pad, the task of
 *                  clipping it and translating it to the current segment falls
 *                  on the subclass. The function should use the segment of data
 *                  and the negotiated media type on the pad to perform
 *                  clipping of input buffer. This function takes ownership of
 *                  buf and should output a buffer or return NULL in
 *                  if the buffer should be dropped.
 * @finish_buffer:  Optional.
 *                  Called when a subclass calls gst_aggregator_finish_buffer()
 *                  from their aggregate function to push out a buffer.
 *                  Subclasses can override this to modify or decorate buffers
 *                  before they get pushed out. This function takes ownership
 *                  of the buffer passed. Subclasses that override this method
 *                  should always chain up to the parent class virtual method.
 * @sink_event:     Optional.
 *                  Called when an event is received on a sink pad, the subclass
 *                  should always chain up.
 * @sink_query:     Optional.
 *                  Called when a query is received on a sink pad, the subclass
 *                  should always chain up.
 * @src_event:      Optional.
 *                  Called when an event is received on the src pad, the subclass
 *                  should always chain up.
 * @src_query:      Optional.
 *                  Called when a query is received on the src pad, the subclass
 *                  should always chain up.
 * @src_activate:   Optional.
 *                  Called when the src pad is activated, it will start/stop its
 *                  pad task right after that call.
 * @aggregate:      Mandatory.
 *                  Called when buffers are queued on all sinkpads. Classes
 *                  should iterate the GstElement->sinkpads and peek or steal
 *                  buffers from the #GstAggregatorPads. If the subclass returns
 *                  GST_FLOW_EOS, sending of the eos event will be taken care
 *                  of. Once / if a buffer has been constructed from the
 *                  aggregated buffers, the subclass should call _finish_buffer.
 * @stop:           Optional.
 *                  Called when the element goes from PAUSED to READY.
 *                  The subclass should free all resources and reset its state.
 * @start:          Optional.
 *                  Called when the element goes from READY to PAUSED.
 *                  The subclass should get ready to process
 *                  aggregated buffers.
 * @get_next_time:  Optional.
 *                  Called when the element needs to know the running time of the next
 *                  rendered buffer for live pipelines. This causes deadline
 *                  based aggregation to occur. Defaults to returning
 *                  GST_CLOCK_TIME_NONE causing the element to wait for buffers
 *                  on all sink pads before aggregating.
 * @create_new_pad: Optional.
 *                  Called when a new pad needs to be created. Allows subclass that
 *                  don't have a single sink pad template to provide a pad based
 *                  on the provided information.
 * @update_src_caps: Lets subclasses update the #GstCaps representing
 *                   the src pad caps before usage.  The result should end up
 *                   in @ret. Return %GST_AGGREGATOR_FLOW_NEED_DATA to indicate that the
 *                   element needs more information (caps, a buffer, etc) to
 *                   choose the correct caps. Should return ANY caps if the
 *                   stream has not caps at all.
 * @fixate_src_caps: Optional.
 *                   Fixate and return the src pad caps provided.  The function takes
 *                   ownership of @caps and returns a fixated version of
 *                   @caps. @caps is not guaranteed to be writable.
 * @negotiated_src_caps: Optional.
 *                       Notifies subclasses what caps format has been negotiated
 * @decide_allocation: Optional.
 *                     Allows the subclass to influence the allocation choices.
 *                     Setup the allocation parameters for allocating output
 *                     buffers. The passed in query contains the result of the
 *                     downstream allocation query.
 * @propose_allocation: Optional.
 *                     Allows the subclass to handle the allocation query from upstream.
 * @negotiate: Optional.
 *             Negotiate the caps with the peer (Since: 1.18).
 * @sink_event_pre_queue: Optional.
 *                        Called when an event is received on a sink pad before queueing up
 *                        serialized events. The subclass should always chain up (Since: 1.18).
 * @sink_query_pre_queue: Optional.
 *                        Called when a query is received on a sink pad before queueing up
 *                        serialized queries. The subclass should always chain up (Since: 1.18).
 *
 * The aggregator base class will handle in a thread-safe way all manners of
 * concurrent flushes, seeks, pad additions and removals, leaving to the
 * subclass the responsibility of clipping buffers, and aggregating buffers in
 * the way the implementor sees fit.
 *
 * It will also take care of event ordering (stream-start, segment, eos).
 *
 * Basically, a simple implementation will override @aggregate, and call
 * _finish_buffer from inside that function.
 *
 * Since: 1.14
 */
struct _GstAggregatorClass {
  GstElementClass   parent_class;

  GstFlowReturn     (*flush)          (GstAggregator    *  aggregator);

  GstBuffer *       (*clip)           (GstAggregator    *  aggregator,
                                       GstAggregatorPad *  aggregator_pad,
                                       GstBuffer        *  buf);

  GstFlowReturn     (*finish_buffer)  (GstAggregator    * aggregator,
                                       GstBuffer        * buffer);

  /* sinkpads virtual methods */
  gboolean          (*sink_event)     (GstAggregator    *  aggregator,
                                       GstAggregatorPad *  aggregator_pad,
                                       GstEvent         *  event);

  gboolean          (*sink_query)     (GstAggregator    *  aggregator,
                                       GstAggregatorPad *  aggregator_pad,
                                       GstQuery         *  query);

  /* srcpad virtual methods */
  gboolean          (*src_event)      (GstAggregator    *  aggregator,
                                       GstEvent         *  event);

  gboolean          (*src_query)      (GstAggregator    *  aggregator,
                                       GstQuery         *  query);

  gboolean          (*src_activate)   (GstAggregator    *  aggregator,
                                       GstPadMode          mode,
                                       gboolean            active);

  GstFlowReturn     (*aggregate)      (GstAggregator    *  aggregator,
                                       gboolean            timeout);

  gboolean          (*stop)           (GstAggregator    *  aggregator);

  gboolean          (*start)          (GstAggregator    *  aggregator);

  GstClockTime      (*get_next_time)  (GstAggregator    *  aggregator);

  GstAggregatorPad * (*create_new_pad) (GstAggregator  * self,
                                        GstPadTemplate * templ,
                                        const gchar    * req_name,
                                        const GstCaps  * caps);

  /**
   * GstAggregatorClass::update_src_caps:
   * @ret: (out) (allow-none):
   */
  GstFlowReturn     (*update_src_caps) (GstAggregator *  self,
                                        GstCaps       *  caps,
                                        GstCaps       ** ret);
  GstCaps *         (*fixate_src_caps) (GstAggregator *  self,
                                        GstCaps       *  caps);
  gboolean          (*negotiated_src_caps) (GstAggregator *  self,
                                            GstCaps      *  caps);
  gboolean          (*decide_allocation) (GstAggregator * self,
                                          GstQuery * query);
  gboolean          (*propose_allocation) (GstAggregator * self,
                                           GstAggregatorPad * pad,
                                           GstQuery * decide_query,
                                           GstQuery * query);

  gboolean          (*negotiate) (GstAggregator * self);

  GstFlowReturn     (*sink_event_pre_queue)     (GstAggregator    *  aggregator,
                                                 GstAggregatorPad *  aggregator_pad,
                                                 GstEvent         *  event);

  gboolean          (*sink_query_pre_queue)     (GstAggregator    *  aggregator,
                                                 GstAggregatorPad *  aggregator_pad,
                                                 GstQuery         *  query);

  /**
   * GstAggregatorClass::finish_buffer_list:
   *
   * Optional. Equivalent of #GstAggregatorClass::finish_buffer for
   * buffer lists.
   *
   * Since: 1.18
   */
  GstFlowReturn     (*finish_buffer_list) (GstAggregator    * aggregator,
                                           GstBufferList    * bufferlist);
  /**
   * GstAggregatorClass::peek_next_sample:
   *
   * See gst_aggregator_peek_next_sample().
   *
   * Since: 1.18
   */
  GstSample *       (*peek_next_sample)         (GstAggregator *aggregator,
                                                 GstAggregatorPad * aggregator_pad);

  /*< private >*/
  gpointer          _gst_reserved[GST_PADDING_LARGE-5];
};

/************************************
 * GstAggregator convenience macros *
 ***********************************/

/**
 * GST_AGGREGATOR_SRC_PAD:
 * @agg: a #GstAggregator
 *
 * Convenience macro to access the source pad of #GstAggregator
 *
 * Since: 1.6
 */
#define GST_AGGREGATOR_SRC_PAD(agg) (((GstAggregator *)(agg))->srcpad)

/*************************
 * GstAggregator methods *
 ************************/

GST_BASE_API
GstFlowReturn  gst_aggregator_finish_buffer         (GstAggregator                *  aggregator,
                                                     GstBuffer                    *  buffer);

GST_BASE_API
GstFlowReturn  gst_aggregator_finish_buffer_list    (GstAggregator                *  aggregator,
                                                     GstBufferList                *  bufferlist);

GST_BASE_API
void           gst_aggregator_set_src_caps          (GstAggregator                *  self,
                                                     GstCaps                      *  caps);

GST_BASE_API
gboolean        gst_aggregator_negotiate            (GstAggregator                * self);

GST_BASE_API
void           gst_aggregator_set_latency           (GstAggregator                *  self,
                                                     GstClockTime                    min_latency,
                                                     GstClockTime                    max_latency);

GST_BASE_API
GType gst_aggregator_get_type(void);

GST_BASE_API
GstClockTime  gst_aggregator_get_latency           (GstAggregator                 *  self);

GST_BASE_API
GstBufferPool * gst_aggregator_get_buffer_pool     (GstAggregator                 * self);

GST_BASE_API
void            gst_aggregator_get_allocator       (GstAggregator                 * self,
                                                    GstAllocator                 ** allocator,
                                                    GstAllocationParams           * params);

GST_BASE_API
GstClockTime    gst_aggregator_simple_get_next_time (GstAggregator                * self);

GST_BASE_API
void            gst_aggregator_update_segment       (GstAggregator                * self,
                                                     const GstSegment             * segment);

GST_BASE_API
GstSample     * gst_aggregator_peek_next_sample     (GstAggregator *self,
                                                     GstAggregatorPad * pad);

GST_BASE_API
void            gst_aggregator_selected_samples     (GstAggregator                * self,
                                                     GstClockTime                   pts,
                                                     GstClockTime                   dts,
                                                     GstClockTime                   duration,
                                                     GstStructure                 * info);

GST_BASE_API
void            gst_aggregator_set_ignore_inactive_pads (GstAggregator * self,
                                                         gboolean ignore);

GST_BASE_API
gboolean        gst_aggregator_get_ignore_inactive_pads (GstAggregator * self);

GST_BASE_API
gboolean        gst_aggregator_get_force_live       (GstAggregator *self);

GST_BASE_API
void            gst_aggregator_set_force_live       (GstAggregator *self,
                                                     gboolean force_live);

/**
 * GstAggregatorStartTimeSelection:
 * @GST_AGGREGATOR_START_TIME_SELECTION_ZERO: Start at running time 0.
 * @GST_AGGREGATOR_START_TIME_SELECTION_FIRST: Start at the running time of
 * the first buffer that is received.
 * @GST_AGGREGATOR_START_TIME_SELECTION_SET: Start at the running time
 * selected by the `start-time` property.
 *
 * Since: 1.18
 */
typedef enum
{
  GST_AGGREGATOR_START_TIME_SELECTION_ZERO,
  GST_AGGREGATOR_START_TIME_SELECTION_FIRST,
  GST_AGGREGATOR_START_TIME_SELECTION_SET
} GstAggregatorStartTimeSelection;

GST_BASE_API
GType           gst_aggregator_start_time_selection_get_type (void);

G_END_DECLS

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAggregator, gst_object_unref)
G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAggregatorPad, gst_object_unref)

#endif /* __GST_AGGREGATOR_H__ */
