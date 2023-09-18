/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
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

#ifndef __GST_BASE_TRANSFORM_H__
#define __GST_BASE_TRANSFORM_H__

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_BASE_TRANSFORM		   (gst_base_transform_get_type())
#define GST_BASE_TRANSFORM(obj)		   (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_BASE_TRANSFORM,GstBaseTransform))
#define GST_BASE_TRANSFORM_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_BASE_TRANSFORM,GstBaseTransformClass))
#define GST_BASE_TRANSFORM_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_BASE_TRANSFORM,GstBaseTransformClass))
#define GST_IS_BASE_TRANSFORM(obj)	   (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_BASE_TRANSFORM))
#define GST_IS_BASE_TRANSFORM_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_BASE_TRANSFORM))
#define GST_BASE_TRANSFORM_CAST(obj)	((GstBaseTransform *)(obj))

/**
 * GST_BASE_TRANSFORM_SINK_NAME:
 *
 * The name of the templates for the sink pad.
 */
#define GST_BASE_TRANSFORM_SINK_NAME	"sink"
/**
 * GST_BASE_TRANSFORM_SRC_NAME:
 *
 * The name of the templates for the source pad.
 */
#define GST_BASE_TRANSFORM_SRC_NAME	"src"

/**
 * GST_BASE_TRANSFORM_SRC_PAD:
 * @obj: base transform instance
 *
 * Gives the pointer to the source #GstPad object of the element.
 */
#define GST_BASE_TRANSFORM_SRC_PAD(obj)		(GST_BASE_TRANSFORM_CAST (obj)->srcpad)

/**
 * GST_BASE_TRANSFORM_SINK_PAD:
 * @obj: base transform instance
 *
 * Gives the pointer to the sink #GstPad object of the element.
 */
#define GST_BASE_TRANSFORM_SINK_PAD(obj)	(GST_BASE_TRANSFORM_CAST (obj)->sinkpad)

/**
 * GST_BASE_TRANSFORM_FLOW_DROPPED:
 *
 * A #GstFlowReturn that can be returned from transform and transform_ip to
 * indicate that no output buffer was generated.
 */
#define GST_BASE_TRANSFORM_FLOW_DROPPED   GST_FLOW_CUSTOM_SUCCESS

typedef struct _GstBaseTransform GstBaseTransform;
typedef struct _GstBaseTransformClass GstBaseTransformClass;
typedef struct _GstBaseTransformPrivate GstBaseTransformPrivate;

/**
 * GstBaseTransform:
 *
 * The opaque #GstBaseTransform data structure.
 */
struct _GstBaseTransform {
  GstElement	 element;

  /*< protected >*/
  /* source and sink pads */
  GstPad	*sinkpad;
  GstPad	*srcpad;

  /* MT-protected (with STREAM_LOCK) */
  gboolean       have_segment;
  GstSegment     segment;
  /* Default submit_input_buffer places the buffer here,
   * for consumption by the generate_output method: */
  GstBuffer      *queued_buf;

  /*< private >*/
  GstBaseTransformPrivate *priv;

  gpointer       _gst_reserved[GST_PADDING_LARGE-1];
};

/**
 * GstBaseTransformClass:
 * @parent_class:   Element parent class
 * @passthrough_on_same_caps: If set to %TRUE, passthrough mode will be
 *                            automatically enabled if the caps are the same.
 *                            Set to %FALSE by default.
 * @transform_ip_on_passthrough: If set to %TRUE, @transform_ip will be called in
 *                           passthrough mode. The passed buffer might not be
 *                           writable. When %FALSE, neither @transform nor
 *                           @transform_ip will be called in passthrough mode.
 *                           Set to %TRUE by default.
 * @transform_caps: Optional.  Given the pad in this direction and the given
 *                  caps, what caps are allowed on the other pad in this
 *                  element ?
 * @fixate_caps:    Optional. Given the pad in this direction and the given
 *                  caps, fixate the caps on the other pad. The function takes
 *                  ownership of @othercaps and returns a fixated version of
 *                  @othercaps. @othercaps is not guaranteed to be writable.
 * @accept_caps:    Optional.
 *                  Subclasses can override this method to check if @caps can be
 *                  handled by the element. The default implementation might not be
 *                  the most optimal way to check this in all cases.
 * @set_caps:       Allows the subclass to be notified of the actual caps set.
 * @query:          Optional.
 *                  Handle a requested query. Subclasses that implement this
 *                  must chain up to the parent if they didn't handle the
 *                  query
 * @decide_allocation: Setup the allocation parameters for allocating output
 *                    buffers. The passed in query contains the result of the
 *                    downstream allocation query. This function is only called
 *                    when not operating in passthrough mode. The default
 *                    implementation will remove all memory dependent metadata.
 *                    If there is a @filter_meta method implementation, it will
 *                    be called for all metadata API in the downstream query,
 *                    otherwise the metadata API is removed.
 * @filter_meta: Return %TRUE if the metadata API should be proposed in the
 *               upstream allocation query. The default implementation is %NULL
 *               and will cause all metadata to be removed.
 * @propose_allocation: Propose buffer allocation parameters for upstream elements.
 *                      This function must be implemented if the element reads or
 *                      writes the buffer content. The query that was passed to
 *                      the decide_allocation is passed in this method (or %NULL
 *                      when the element is in passthrough mode). The default
 *                      implementation will pass the query downstream when in
 *                      passthrough mode and will copy all the filtered metadata
 *                      API in non-passthrough mode.
 * @transform_size: Optional. Given the size of a buffer in the given direction
 *                  with the given caps, calculate the size in bytes of a buffer
 *                  on the other pad with the given other caps.
 *                  The default implementation uses get_unit_size and keeps
 *                  the number of units the same.
 * @get_unit_size:  Required if the transform is not in-place.
 *                  Get the size in bytes of one unit for the given caps.
 * @start:          Optional.
 *                  Called when the element starts processing.
 *                  Allows opening external resources.
 * @stop:           Optional.
 *                  Called when the element stops processing.
 *                  Allows closing external resources.
 * @sink_event:     Optional.
 *                  Event handler on the sink pad. The default implementation
 *                  handles the event and forwards it downstream.
 * @src_event:      Optional.
 *                  Event handler on the source pad. The default implementation
 *                  handles the event and forwards it upstream.
 * @prepare_output_buffer: Optional.
 *                         Subclasses can override this to do their own
 *                         allocation of output buffers.  Elements that only do
 *                         analysis can return a subbuffer or even just
 *                         return a reference to the input buffer (if in
 *                         passthrough mode). The default implementation will
 *                         use the negotiated allocator or bufferpool and
 *                         transform_size to allocate an output buffer or it
 *                         will return the input buffer in passthrough mode.
 * @copy_metadata: Optional.
 *                 Copy the metadata from the input buffer to the output buffer.
 *                 The default implementation will copy the flags, timestamps and
 *                 offsets of the buffer.
 * @transform_meta: Optional. Transform the metadata on the input buffer to the
 *                  output buffer. By default this method copies all meta without
 *                  tags. Subclasses can implement this method and return %TRUE if
 *                  the metadata is to be copied.
 * @before_transform: Optional.
 *                    This method is called right before the base class will
 *                    start processing. Dynamic properties or other delayed
 *                    configuration could be performed in this method.
 * @transform:      Required if the element does not operate in-place.
 *                  Transforms one incoming buffer to one outgoing buffer.
 *                  The function is allowed to change size/timestamp/duration
 *                  of the outgoing buffer.
 * @transform_ip:   Required if the element operates in-place.
 *                  Transform the incoming buffer in-place.
 * @submit_input_buffer: Function which accepts a new input buffer and pre-processes it.
 *                  The default implementation performs caps (re)negotiation, then
 *                  QoS if needed, and places the input buffer into the @queued_buf
 *                  member variable. If the buffer is dropped due to QoS, it returns
 *                  GST_BASE_TRANSFORM_FLOW_DROPPED. If this input buffer is not
 *                  contiguous with any previous input buffer, then @is_discont
 *                  is set to %TRUE. (Since: 1.6)
 * @generate_output: Called after each new input buffer is submitted repeatedly
 *                   until it either generates an error or fails to generate an output
 *                   buffer. The default implementation takes the contents of the
 *                   @queued_buf variable, generates an output buffer if needed
 *                   by calling the class @prepare_output_buffer, and then
 *                   calls either @transform or @transform_ip. Elements that don't
 *                   do 1-to-1 transformations of input to output buffers can either
 *                   return GST_BASE_TRANSFORM_FLOW_DROPPED or simply not generate
 *                   an output buffer until they are ready to do so. (Since: 1.6)
 *
 * Subclasses can override any of the available virtual methods or not, as
 * needed. At minimum either @transform or @transform_ip need to be overridden.
 * If the element can overwrite the input data with the results (data is of the
 * same type and quantity) it should provide @transform_ip.
 */
struct _GstBaseTransformClass {
  GstElementClass parent_class;

  /*< public >*/
  gboolean       passthrough_on_same_caps;
  gboolean       transform_ip_on_passthrough;

  /* virtual methods for subclasses */
  GstCaps*	(*transform_caps) (GstBaseTransform *trans,
                                   GstPadDirection direction,
                                   GstCaps *caps, GstCaps *filter);
  /**
   * GstBaseTransformClass::fixate_caps:
   * @othercaps: (transfer full):
   */
  GstCaps*	(*fixate_caps)	  (GstBaseTransform *trans,
                                   GstPadDirection direction, GstCaps *caps,
                                   GstCaps *othercaps);
  gboolean      (*accept_caps)    (GstBaseTransform *trans, GstPadDirection direction,
                                   GstCaps *caps);
  gboolean      (*set_caps)       (GstBaseTransform *trans, GstCaps *incaps,
                                   GstCaps *outcaps);
  gboolean      (*query)          (GstBaseTransform *trans, GstPadDirection direction,
                                   GstQuery *query);

  /* decide allocation query for output buffers */
  gboolean      (*decide_allocation)  (GstBaseTransform *trans, GstQuery *query);
  gboolean      (*filter_meta)        (GstBaseTransform *trans, GstQuery *query,
                                       GType api, const GstStructure *params);

  /* propose allocation query parameters for input buffers */
  gboolean      (*propose_allocation) (GstBaseTransform *trans, GstQuery *decide_query,
                                       GstQuery *query);

  /**
   * GstBaseTransformClass::transform_size:
   * @othersize: (out):
   */
  gboolean      (*transform_size) (GstBaseTransform *trans,
                                   GstPadDirection direction,
                                   GstCaps *caps, gsize size,
                                   GstCaps *othercaps, gsize *othersize);

  /**
   * GstBaseTransformClass::get_unit_size:
   * @size: (out):
   */
  gboolean      (*get_unit_size)  (GstBaseTransform *trans, GstCaps *caps,
                                   gsize *size);

  /* states */
  gboolean      (*start)        (GstBaseTransform *trans);
  gboolean      (*stop)         (GstBaseTransform *trans);

  /* sink and src pad event handlers */
  /**
   * GstBaseTransformClass::sink_event:
   * @event: (transfer full):
   */
  gboolean      (*sink_event)   (GstBaseTransform *trans, GstEvent *event);
  /**
   * GstBaseTransformClass::src_event:
   * @event: (transfer full):
   */
  gboolean      (*src_event)    (GstBaseTransform *trans, GstEvent *event);

  /**
   * GstBaseTransformClass::prepare_output_buffer:
   * @outbuf: (out):
   */
  GstFlowReturn (*prepare_output_buffer) (GstBaseTransform * trans,
                                          GstBuffer *input, GstBuffer **outbuf);

  /* metadata */
  gboolean      (*copy_metadata)     (GstBaseTransform *trans, GstBuffer *input,
                                      GstBuffer *outbuf);
  gboolean      (*transform_meta)    (GstBaseTransform *trans, GstBuffer *outbuf,
                                      GstMeta *meta, GstBuffer *inbuf);

  void          (*before_transform)  (GstBaseTransform *trans, GstBuffer *buffer);

  /* transform */
  GstFlowReturn (*transform)    (GstBaseTransform *trans, GstBuffer *inbuf,
                                 GstBuffer *outbuf);
  GstFlowReturn (*transform_ip) (GstBaseTransform *trans, GstBuffer *buf);

  GstFlowReturn (*submit_input_buffer) (GstBaseTransform *trans, gboolean is_discont, GstBuffer *input);

  /**
   * GstBaseTransformClass::generate_output:
   * @outbuf: (out):
   */
  GstFlowReturn (*generate_output) (GstBaseTransform *trans, GstBuffer **outbuf);

  /*< private >*/
  gpointer       _gst_reserved[GST_PADDING_LARGE - 2];
};

GST_BASE_API
GType           gst_base_transform_get_type         (void);

GST_BASE_API
void		gst_base_transform_set_passthrough  (GstBaseTransform *trans,
	                                             gboolean passthrough);
GST_BASE_API
gboolean	gst_base_transform_is_passthrough   (GstBaseTransform *trans);

GST_BASE_API
void		gst_base_transform_set_in_place     (GstBaseTransform *trans,
	                                             gboolean in_place);
GST_BASE_API
gboolean	gst_base_transform_is_in_place      (GstBaseTransform *trans);

GST_BASE_API
void		gst_base_transform_update_qos       (GstBaseTransform *trans,
						     gdouble proportion,
						     GstClockTimeDiff diff,
						     GstClockTime timestamp);
GST_BASE_API
void		gst_base_transform_set_qos_enabled  (GstBaseTransform *trans,
		                                     gboolean enabled);
GST_BASE_API
gboolean	gst_base_transform_is_qos_enabled   (GstBaseTransform *trans);

GST_BASE_API
void            gst_base_transform_set_gap_aware    (GstBaseTransform *trans,
                                                     gboolean gap_aware);
GST_BASE_API
void            gst_base_transform_set_prefer_passthrough (GstBaseTransform *trans,
                                                           gboolean prefer_passthrough);
GST_BASE_API
GstBufferPool * gst_base_transform_get_buffer_pool  (GstBaseTransform *trans);

GST_BASE_API
void            gst_base_transform_get_allocator    (GstBaseTransform *trans,
                                                     GstAllocator **allocator,
                                                     GstAllocationParams *params);
GST_BASE_API
void		gst_base_transform_reconfigure_sink (GstBaseTransform *trans);

GST_BASE_API
void		gst_base_transform_reconfigure_src  (GstBaseTransform *trans);

GST_BASE_API
gboolean gst_base_transform_update_src_caps (GstBaseTransform *trans,
                                             GstCaps *updated_caps);

GST_BASE_API
gboolean gst_base_transform_reconfigure (GstBaseTransform * trans);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBaseTransform, gst_object_unref)

G_END_DECLS

#endif /* __GST_BASE_TRANSFORM_H__ */
