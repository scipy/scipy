/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wim.taymans@chello.be>
 *
 * gstpad.h: Header for GstPad object
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


#ifndef __GST_PAD_H__
#define __GST_PAD_H__

#include <gst/gstconfig.h>

typedef struct _GstPad GstPad;
typedef struct _GstPadPrivate GstPadPrivate;
typedef struct _GstPadClass GstPadClass;
typedef struct _GstPadProbeInfo GstPadProbeInfo;

/**
 * GstPadDirection:
 * @GST_PAD_UNKNOWN: direction is unknown.
 * @GST_PAD_SRC: the pad is a source pad.
 * @GST_PAD_SINK: the pad is a sink pad.
 *
 * The direction of a pad.
 */
typedef enum {
  GST_PAD_UNKNOWN,
  GST_PAD_SRC,
  GST_PAD_SINK
} GstPadDirection;

/**
 * GstPadMode:
 * @GST_PAD_MODE_NONE: Pad will not handle dataflow
 * @GST_PAD_MODE_PUSH: Pad handles dataflow in downstream push mode
 * @GST_PAD_MODE_PULL: Pad handles dataflow in upstream pull mode
 *
 * The status of a GstPad. After activating a pad, which usually happens when the
 * parent element goes from READY to PAUSED, the GstPadMode defines if the
 * pad operates in push or pull mode.
 */
typedef enum {
  GST_PAD_MODE_NONE,
  GST_PAD_MODE_PUSH,
  GST_PAD_MODE_PULL
} GstPadMode;

#include <glib.h>

GST_API
const gchar   * gst_pad_mode_get_name (GstPadMode mode);

#include <gst/gstobject.h>
#include <gst/gstbuffer.h>
#include <gst/gstbufferlist.h>
#include <gst/gstcaps.h>
#include <gst/gstpadtemplate.h>
#include <gst/gstevent.h>
#include <gst/gstquery.h>
#include <gst/gsttask.h>

G_BEGIN_DECLS

/*
 * Pad base class
 */
#define GST_TYPE_PAD			(gst_pad_get_type ())
#define GST_IS_PAD(obj)			(G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_PAD))
#define GST_IS_PAD_CLASS(klass)		(G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_PAD))
#define GST_PAD(obj)			(G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_PAD, GstPad))
#define GST_PAD_CLASS(klass)		(G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_PAD, GstPadClass))
#define GST_PAD_CAST(obj)		((GstPad*)(obj))



/**
 * GstPadLinkReturn:
 * @GST_PAD_LINK_OK		: link succeeded
 * @GST_PAD_LINK_WRONG_HIERARCHY: pads have no common grandparent
 * @GST_PAD_LINK_WAS_LINKED	: pad was already linked
 * @GST_PAD_LINK_WRONG_DIRECTION: pads have wrong direction
 * @GST_PAD_LINK_NOFORMAT	: pads do not have common format
 * @GST_PAD_LINK_NOSCHED	: pads cannot cooperate in scheduling
 * @GST_PAD_LINK_REFUSED	: refused for some reason
 *
 * Result values from gst_pad_link and friends.
 */
typedef enum {
  GST_PAD_LINK_OK               =  0,
  GST_PAD_LINK_WRONG_HIERARCHY  = -1,
  GST_PAD_LINK_WAS_LINKED       = -2,
  GST_PAD_LINK_WRONG_DIRECTION  = -3,
  GST_PAD_LINK_NOFORMAT         = -4,
  GST_PAD_LINK_NOSCHED          = -5,
  GST_PAD_LINK_REFUSED          = -6
} GstPadLinkReturn;

/**
 * GST_PAD_LINK_FAILED:
 * @ret: the #GstPadLinkReturn value
 *
 * Macro to test if the given #GstPadLinkReturn value indicates a failed
 * link step.
 */
#define GST_PAD_LINK_FAILED(ret) ((ret) < GST_PAD_LINK_OK)

/**
 * GST_PAD_LINK_SUCCESSFUL:
 * @ret: the #GstPadLinkReturn value
 *
 * Macro to test if the given #GstPadLinkReturn value indicates a successful
 * link step.
 */
#define GST_PAD_LINK_SUCCESSFUL(ret) ((ret) >= GST_PAD_LINK_OK)

/**
 * GstFlowReturn:
 * @GST_FLOW_OK:		 Data passing was ok.
 * @GST_FLOW_NOT_LINKED:	 Pad is not linked.
 * @GST_FLOW_FLUSHING:	         Pad is flushing.
 * @GST_FLOW_EOS:                Pad is EOS.
 * @GST_FLOW_NOT_NEGOTIATED:	 Pad is not negotiated.
 * @GST_FLOW_ERROR:		 Some (fatal) error occurred. Element generating
 *                               this error should post an error message using
 *                               GST_ELEMENT_ERROR() with more details.
 * @GST_FLOW_NOT_SUPPORTED:	 This operation is not supported.
 * @GST_FLOW_CUSTOM_SUCCESS:	 Elements can use values starting from
 *                               this (and higher) to define custom success
 *                               codes.
 * @GST_FLOW_CUSTOM_SUCCESS_1:	 Pre-defined custom success code (define your
 *                               custom success code to this to avoid compiler
 *                               warnings).
 * @GST_FLOW_CUSTOM_SUCCESS_2:	 Pre-defined custom success code.
 * @GST_FLOW_CUSTOM_ERROR:	 Elements can use values starting from
 *                               this (and lower) to define custom error codes.
 * @GST_FLOW_CUSTOM_ERROR_1:	 Pre-defined custom error code (define your
 *                               custom error code to this to avoid compiler
 *                               warnings).
 * @GST_FLOW_CUSTOM_ERROR_2:	 Pre-defined custom error code.
 *
 * The result of passing data to a pad.
 *
 * Note that the custom return values should not be exposed outside of the
 * element scope.
 */
typedef enum {
  /* custom success starts here */
  GST_FLOW_CUSTOM_SUCCESS_2 = 102,
  GST_FLOW_CUSTOM_SUCCESS_1 = 101,
  GST_FLOW_CUSTOM_SUCCESS = 100,

  /* core predefined */
  GST_FLOW_OK		  =  0,
  /* expected failures */
  GST_FLOW_NOT_LINKED     = -1,
  GST_FLOW_FLUSHING       = -2,
  /* error cases */
  GST_FLOW_EOS            = -3,
  GST_FLOW_NOT_NEGOTIATED = -4,
  GST_FLOW_ERROR	  = -5,
  GST_FLOW_NOT_SUPPORTED  = -6,

  /* custom error starts here */
  GST_FLOW_CUSTOM_ERROR   = -100,
  GST_FLOW_CUSTOM_ERROR_1 = -101,
  GST_FLOW_CUSTOM_ERROR_2 = -102
} GstFlowReturn;

GST_API
const gchar *         gst_flow_get_name (GstFlowReturn ret);

GST_API
GQuark                gst_flow_to_quark (GstFlowReturn ret);

GST_API
const gchar *         gst_pad_link_get_name (GstPadLinkReturn ret);

/**
 * GstPadLinkCheck:
 * @GST_PAD_LINK_CHECK_NOTHING: Don't check hierarchy or caps compatibility.
 * @GST_PAD_LINK_CHECK_HIERARCHY: Check the pads have same parents/grandparents.
 *   Could be omitted if it is already known that the two elements that own the
 *   pads are in the same bin.
 * @GST_PAD_LINK_CHECK_TEMPLATE_CAPS: Check if the pads are compatible by using
 *   their template caps. This is much faster than @GST_PAD_LINK_CHECK_CAPS, but
 *   would be unsafe e.g. if one pad has %GST_CAPS_ANY.
 * @GST_PAD_LINK_CHECK_CAPS: Check if the pads are compatible by comparing the
 *   caps returned by gst_pad_query_caps().
 * @GST_PAD_LINK_CHECK_NO_RECONFIGURE: Disables pushing a reconfigure event when pads are
 *   linked.
 * @GST_PAD_LINK_CHECK_DEFAULT: The default checks done when linking
 *   pads (i.e. the ones used by gst_pad_link()).
 *
 * The amount of checking to be done when linking pads. @GST_PAD_LINK_CHECK_CAPS
 * and @GST_PAD_LINK_CHECK_TEMPLATE_CAPS are mutually exclusive. If both are
 * specified, expensive but safe @GST_PAD_LINK_CHECK_CAPS are performed.
 *
 * > Only disable some of the checks if you are 100% certain you know the link
 * > will not fail because of hierarchy/caps compatibility failures. If uncertain,
 * > use the default checks (%GST_PAD_LINK_CHECK_DEFAULT) or the regular methods
 * > for linking the pads.
 */

typedef enum {
  GST_PAD_LINK_CHECK_NOTHING       = 0,
  GST_PAD_LINK_CHECK_HIERARCHY     = 1 << 0,
  GST_PAD_LINK_CHECK_TEMPLATE_CAPS = 1 << 1,
  GST_PAD_LINK_CHECK_CAPS          = 1 << 2,


  /* Not really checks, more like flags
   * Added here to avoid creating a new gst_pad_link_variant */
  GST_PAD_LINK_CHECK_NO_RECONFIGURE = 1 << 3,

  GST_PAD_LINK_CHECK_DEFAULT       = GST_PAD_LINK_CHECK_HIERARCHY | GST_PAD_LINK_CHECK_CAPS
} GstPadLinkCheck;

/* pad states */
/**
 * GstPadActivateFunction:
 * @pad: a #GstPad
 * @parent: the parent of @pad
 *
 * This function is called when the pad is activated during the element
 * READY to PAUSED state change. By default this function will call the
 * activate function that puts the pad in push mode but elements can
 * override this function to activate the pad in pull mode if they wish.
 *
 * Returns: %TRUE if the pad could be activated.
 */
typedef gboolean		(*GstPadActivateFunction)	(GstPad *pad, GstObject *parent);
/**
 * GstPadActivateModeFunction:
 * @pad: a #GstPad
 * @parent: the parent of @pad
 * @mode: the requested activation mode of @pad
 * @active: activate or deactivate the pad.
 *
 * The prototype of the push and pull activate functions.
 *
 * Returns: %TRUE if the pad could be activated or deactivated.
 */
typedef gboolean		(*GstPadActivateModeFunction)	(GstPad *pad, GstObject *parent,
                                                                 GstPadMode mode, gboolean active);


/* data passing */
/**
 * GstPadChainFunction:
 * @pad: the sink #GstPad that performed the chain.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @buffer: (transfer full): the #GstBuffer that is chained, not %NULL.
 *
 * A function that will be called on sinkpads when chaining buffers.
 * The function typically processes the data contained in the buffer and
 * either consumes the data or passes it on to the internally linked pad(s).
 *
 * The implementer of this function receives a refcount to @buffer and should
 * gst_buffer_unref() when the buffer is no longer needed.
 *
 * When a chain function detects an error in the data stream, it must post an
 * error on the bus and return an appropriate #GstFlowReturn value.
 *
 * Returns: #GST_FLOW_OK for success
 */
typedef GstFlowReturn		(*GstPadChainFunction)		(GstPad *pad, GstObject *parent,
                                                                 GstBuffer *buffer);

/**
 * GstPadChainListFunction:
 * @pad: the sink #GstPad that performed the chain.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @list: (transfer full): the #GstBufferList that is chained, not %NULL.
 *
 * A function that will be called on sinkpads when chaining buffer lists.
 * The function typically processes the data contained in the buffer list and
 * either consumes the data or passes it on to the internally linked pad(s).
 *
 * The implementer of this function receives a refcount to @list and
 * should gst_buffer_list_unref() when the list is no longer needed.
 *
 * When a chainlist function detects an error in the data stream, it must
 * post an error on the bus and return an appropriate #GstFlowReturn value.
 *
 * Returns: #GST_FLOW_OK for success
 */
typedef GstFlowReturn		(*GstPadChainListFunction)	(GstPad *pad, GstObject *parent,
                                                                 GstBufferList *list);

/**
 * GstPadGetRangeFunction:
 * @pad: the src #GstPad to perform the getrange on.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @offset: the offset of the range
 * @length: the length of the range
 * @buffer: a memory location to hold the result buffer, cannot be %NULL.
 *
 * This function will be called on source pads when a peer element
 * request a buffer at the specified @offset and @length. If this function
 * returns #GST_FLOW_OK, the result buffer will be stored in @buffer. The
 * contents of @buffer is invalid for any other return value.
 *
 * This function is installed on a source pad with
 * gst_pad_set_getrange_function() and can only be called on source pads after
 * they are successfully activated with gst_pad_activate_mode() with the
 * #GST_PAD_MODE_PULL.
 *
 * @offset and @length are always given in byte units. @offset must normally be a value
 * between 0 and the length in bytes of the data available on @pad. The
 * length (duration in bytes) can be retrieved with a #GST_QUERY_DURATION or with a
 * #GST_QUERY_SEEKING.
 *
 * Any @offset larger or equal than the length will make the function return
 * #GST_FLOW_EOS, which corresponds to EOS. In this case @buffer does not
 * contain a valid buffer.
 *
 * The buffer size of @buffer will only be smaller than @length when @offset is
 * near the end of the stream. In all other cases, the size of @buffer must be
 * exactly the requested size.
 *
 * It is allowed to call this function with a 0 @length and valid @offset, in
 * which case @buffer will contain a 0-sized buffer and the function returns
 * #GST_FLOW_OK.
 *
 * When this function is called with a -1 @offset, the sequentially next buffer
 * of length @length in the stream is returned.
 *
 * When this function is called with a -1 @length, a buffer with a default
 * optimal length is returned in @buffer. The length might depend on the value
 * of @offset.
 *
 * Returns: #GST_FLOW_OK for success and a valid buffer in @buffer. Any other
 * return value leaves @buffer undefined.
 */
typedef GstFlowReturn		(*GstPadGetRangeFunction)	(GstPad *pad, GstObject *parent,
                                                                 guint64 offset, guint length,
                                                                 GstBuffer **buffer);

/**
 * GstPadEventFunction:
 * @pad: the #GstPad to handle the event.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @event: (transfer full): the #GstEvent to handle.
 *
 * Function signature to handle an event for the pad.
 *
 * Returns: %TRUE if the pad could handle the event.
 */
typedef gboolean		(*GstPadEventFunction)		(GstPad *pad, GstObject *parent,
                                                                 GstEvent *event);

/**
 * GstPadEventFullFunction:
 * @pad: the #GstPad to handle the event.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @event: (transfer full): the #GstEvent to handle.
 *
 * Function signature to handle an event for the pad.
 *
 * This variant is for specific elements that will take into account the
 * last downstream flow return (from a pad push), in which case they can
 * return it.
 *
 * Returns: %GST_FLOW_OK if the event was handled properly, or any other
 * #GstFlowReturn dependent on downstream state.
 *
 * Since: 1.8
 */
typedef GstFlowReturn		(*GstPadEventFullFunction)	(GstPad *pad, GstObject *parent,
                                                                 GstEvent *event);


/* internal links */
/**
 * GstPadIterIntLinkFunction:
 * @pad: The #GstPad to query.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 *
 * The signature of the internal pad link iterator function.
 *
 * Returns: a new #GstIterator that will iterate over all pads that are
 * linked to the given pad on the inside of the parent element.
 *
 * the caller must call gst_iterator_free() after usage.
 */
typedef GstIterator*           (*GstPadIterIntLinkFunction)    (GstPad *pad, GstObject *parent);

/* generic query function */
/**
 * GstPadQueryFunction:
 * @pad: the #GstPad to query.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @query: the #GstQuery object to execute
 *
 * The signature of the query function.
 *
 * Returns: %TRUE if the query could be performed.
 */
typedef gboolean		(*GstPadQueryFunction)		(GstPad *pad, GstObject *parent,
                                                                 GstQuery *query);


/* linking */
/**
 * GstPadLinkFunction:
 * @pad: the #GstPad that is linked.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 * @peer: the peer #GstPad of the link
 *
 * Function signature to handle a new link on the pad.
 *
 * Returns: the result of the link with the specified peer.
 */
typedef GstPadLinkReturn	(*GstPadLinkFunction)		(GstPad *pad, GstObject *parent, GstPad *peer);
/**
 * GstPadUnlinkFunction:
 * @pad: the #GstPad that is linked.
 * @parent: (allow-none): the parent of @pad. If the #GST_PAD_FLAG_NEED_PARENT
 *          flag is set, @parent is guaranteed to be not-%NULL and remain valid
 *          during the execution of this function.
 *
 * Function signature to handle a unlinking the pad prom its peer.
 *
 * The pad's lock is already held when the unlink function is called, so most
 * pad functions cannot be called from within the callback.
 */
typedef void			(*GstPadUnlinkFunction)		(GstPad *pad, GstObject *parent);


/* misc */
/**
 * GstPadForwardFunction:
 * @pad: the #GstPad that is forwarded.
 * @user_data: the gpointer to optional user data.
 *
 * A forward function is called for all internally linked pads, see
 * gst_pad_forward().
 *
 * Returns: %TRUE if the dispatching procedure has to be stopped.
 */
typedef gboolean		(*GstPadForwardFunction)	(GstPad *pad, gpointer user_data);

/**
 * GstPadProbeType:
 * @GST_PAD_PROBE_TYPE_INVALID: invalid probe type
 * @GST_PAD_PROBE_TYPE_IDLE: probe idle pads and block while the callback is called
 * @GST_PAD_PROBE_TYPE_BLOCK: probe and block pads
 * @GST_PAD_PROBE_TYPE_BUFFER: probe buffers
 * @GST_PAD_PROBE_TYPE_BUFFER_LIST: probe buffer lists
 * @GST_PAD_PROBE_TYPE_EVENT_DOWNSTREAM: probe downstream events
 * @GST_PAD_PROBE_TYPE_EVENT_UPSTREAM: probe upstream events
 * @GST_PAD_PROBE_TYPE_EVENT_FLUSH: probe flush events. This probe has to be
 *     explicitly enabled and is not included in the
 *     @@GST_PAD_PROBE_TYPE_EVENT_DOWNSTREAM or
 *     @@GST_PAD_PROBE_TYPE_EVENT_UPSTREAM probe types.
 * @GST_PAD_PROBE_TYPE_QUERY_DOWNSTREAM: probe downstream queries
 * @GST_PAD_PROBE_TYPE_QUERY_UPSTREAM: probe upstream queries
 * @GST_PAD_PROBE_TYPE_PUSH: probe push
 * @GST_PAD_PROBE_TYPE_PULL: probe pull
 * @GST_PAD_PROBE_TYPE_BLOCKING: probe and block at the next opportunity, at data flow or when idle
 * @GST_PAD_PROBE_TYPE_DATA_DOWNSTREAM: probe downstream data (buffers, buffer lists, and events)
 * @GST_PAD_PROBE_TYPE_DATA_UPSTREAM: probe upstream data (events)
 * @GST_PAD_PROBE_TYPE_DATA_BOTH: probe upstream and downstream data (buffers, buffer lists, and events)
 * @GST_PAD_PROBE_TYPE_BLOCK_DOWNSTREAM: probe and block downstream data (buffers, buffer lists, and events)
 * @GST_PAD_PROBE_TYPE_BLOCK_UPSTREAM: probe and block upstream data (events)
 * @GST_PAD_PROBE_TYPE_EVENT_BOTH: probe upstream and downstream events
 * @GST_PAD_PROBE_TYPE_QUERY_BOTH: probe upstream and downstream queries
 * @GST_PAD_PROBE_TYPE_ALL_BOTH: probe upstream events and queries and downstream buffers, buffer lists, events and queries
 * @GST_PAD_PROBE_TYPE_SCHEDULING: probe push and pull
 *
 * The different probing types that can occur. When either one of
 * @GST_PAD_PROBE_TYPE_IDLE or @GST_PAD_PROBE_TYPE_BLOCK is used, the probe will be a
 * blocking probe.
 */
typedef enum
{
  GST_PAD_PROBE_TYPE_INVALID          = 0,
  /* flags to control blocking */
  GST_PAD_PROBE_TYPE_IDLE             = (1 << 0),
  GST_PAD_PROBE_TYPE_BLOCK            = (1 << 1),
  /* flags to select datatypes */
  GST_PAD_PROBE_TYPE_BUFFER           = (1 << 4),
  GST_PAD_PROBE_TYPE_BUFFER_LIST      = (1 << 5),
  GST_PAD_PROBE_TYPE_EVENT_DOWNSTREAM = (1 << 6),
  GST_PAD_PROBE_TYPE_EVENT_UPSTREAM   = (1 << 7),
  GST_PAD_PROBE_TYPE_EVENT_FLUSH      = (1 << 8),
  GST_PAD_PROBE_TYPE_QUERY_DOWNSTREAM = (1 << 9),
  GST_PAD_PROBE_TYPE_QUERY_UPSTREAM   = (1 << 10),
  /* flags to select scheduling mode */
  GST_PAD_PROBE_TYPE_PUSH             = (1 << 12),
  GST_PAD_PROBE_TYPE_PULL             = (1 << 13),

  /* flag combinations */
  GST_PAD_PROBE_TYPE_BLOCKING         = GST_PAD_PROBE_TYPE_IDLE | GST_PAD_PROBE_TYPE_BLOCK,
  GST_PAD_PROBE_TYPE_DATA_DOWNSTREAM  = GST_PAD_PROBE_TYPE_BUFFER | GST_PAD_PROBE_TYPE_BUFFER_LIST | GST_PAD_PROBE_TYPE_EVENT_DOWNSTREAM,
  GST_PAD_PROBE_TYPE_DATA_UPSTREAM    = GST_PAD_PROBE_TYPE_EVENT_UPSTREAM,
  GST_PAD_PROBE_TYPE_DATA_BOTH        = GST_PAD_PROBE_TYPE_DATA_DOWNSTREAM | GST_PAD_PROBE_TYPE_DATA_UPSTREAM,
  GST_PAD_PROBE_TYPE_BLOCK_DOWNSTREAM = GST_PAD_PROBE_TYPE_BLOCK | GST_PAD_PROBE_TYPE_DATA_DOWNSTREAM,
  GST_PAD_PROBE_TYPE_BLOCK_UPSTREAM   = GST_PAD_PROBE_TYPE_BLOCK | GST_PAD_PROBE_TYPE_DATA_UPSTREAM,
  GST_PAD_PROBE_TYPE_EVENT_BOTH       = GST_PAD_PROBE_TYPE_EVENT_DOWNSTREAM | GST_PAD_PROBE_TYPE_EVENT_UPSTREAM,
  GST_PAD_PROBE_TYPE_QUERY_BOTH       = GST_PAD_PROBE_TYPE_QUERY_DOWNSTREAM | GST_PAD_PROBE_TYPE_QUERY_UPSTREAM,
  GST_PAD_PROBE_TYPE_ALL_BOTH         = GST_PAD_PROBE_TYPE_DATA_BOTH | GST_PAD_PROBE_TYPE_QUERY_BOTH,
  GST_PAD_PROBE_TYPE_SCHEDULING       = GST_PAD_PROBE_TYPE_PUSH | GST_PAD_PROBE_TYPE_PULL
} GstPadProbeType;


/**
 * GstPadProbeReturn:
 * @GST_PAD_PROBE_OK: normal probe return value. This leaves the probe in
 *        place, and defers decisions about dropping or passing data to other
 *        probes, if any. If there are no other probes, the default behaviour
 *        for the probe type applies ('block' for blocking probes,
 *        and 'pass' for non-blocking probes).
 * @GST_PAD_PROBE_DROP: drop data in data probes. For push mode this means that
 *        the data item is not sent downstream. For pull mode, it means that
 *        the data item is not passed upstream. In both cases, no other probes
 *        are called for this item and %GST_FLOW_OK or %TRUE is returned to the
 *        caller.
 * @GST_PAD_PROBE_REMOVE: remove this probe, passing the data. For blocking probes
 *        this will cause data flow to unblock, unless there are also other
 *        blocking probes installed.
 * @GST_PAD_PROBE_PASS: pass the data item in the block probe and block on the
 *        next item. Note, that if there are multiple pad probes installed and
 *        any probe returns PASS, the data will be passed.
 * @GST_PAD_PROBE_HANDLED: Data has been handled in the probe and will not be
 *        forwarded further. For events and buffers this is the same behaviour as
 *        %GST_PAD_PROBE_DROP (except that in this case you need to unref the buffer
 *        or event yourself). For queries it will also return %TRUE to the caller.
 *        The probe can also modify the #GstFlowReturn value by using the
 *        #GST_PAD_PROBE_INFO_FLOW_RETURN() accessor.
 *        Note that the resulting query must contain valid entries.
 *        Since: 1.6
 *
 * Different return values for the #GstPadProbeCallback.
 */
typedef enum
{
  GST_PAD_PROBE_DROP,
  GST_PAD_PROBE_OK,
  GST_PAD_PROBE_REMOVE,
  GST_PAD_PROBE_PASS,
  GST_PAD_PROBE_HANDLED
} GstPadProbeReturn;


/**
 * GstPadProbeInfo:
 * @type: the current probe type
 * @id: the id of the probe
 * @data: (allow-none): type specific data, check the @type field to know the
 *    datatype.  This field can be %NULL.
 * @offset: offset of pull probe, this field is valid when @type contains
 *    #GST_PAD_PROBE_TYPE_PULL
 * @size: size of pull probe, this field is valid when @type contains
 *    #GST_PAD_PROBE_TYPE_PULL
 *
 * Info passed in the #GstPadProbeCallback.
 */
struct _GstPadProbeInfo
{
  GstPadProbeType type;
  gulong id;
  gpointer data;
  guint64 offset;
  guint size;

  /*< private >*/
  union {
    gpointer _gst_reserved[GST_PADDING];
    struct {
      GstFlowReturn flow_ret;
    } abi;
  } ABI;
};

#define GST_PAD_PROBE_INFO_TYPE(d)         ((d)->type)
#define GST_PAD_PROBE_INFO_ID(d)           ((d)->id)
#define GST_PAD_PROBE_INFO_DATA(d)         ((d)->data)
#define GST_PAD_PROBE_INFO_FLOW_RETURN(d)  ((d)->ABI.abi.flow_ret)

#define GST_PAD_PROBE_INFO_BUFFER(d)       GST_BUFFER_CAST(GST_PAD_PROBE_INFO_DATA(d))
#define GST_PAD_PROBE_INFO_BUFFER_LIST(d)  GST_BUFFER_LIST_CAST(GST_PAD_PROBE_INFO_DATA(d))
#define GST_PAD_PROBE_INFO_EVENT(d)        GST_EVENT_CAST(GST_PAD_PROBE_INFO_DATA(d))
#define GST_PAD_PROBE_INFO_QUERY(d)        GST_QUERY_CAST(GST_PAD_PROBE_INFO_DATA(d))

#define GST_PAD_PROBE_INFO_OFFSET(d)       ((d)->offset)
#define GST_PAD_PROBE_INFO_SIZE(d)         ((d)->size)

GST_API
GstEvent*      gst_pad_probe_info_get_event       (GstPadProbeInfo * info);

GST_API
GstQuery*      gst_pad_probe_info_get_query       (GstPadProbeInfo * info);

GST_API
GstBuffer*     gst_pad_probe_info_get_buffer      (GstPadProbeInfo * info);

GST_API
GstBufferList* gst_pad_probe_info_get_buffer_list (GstPadProbeInfo * info);

/**
 * GstPadProbeCallback:
 * @pad: the #GstPad that is blocked
 * @info: #GstPadProbeInfo
 * @user_data: the gpointer to optional user data.
 *
 * Callback used by gst_pad_add_probe(). Gets called to notify about the current
 * blocking type.
 *
 * The callback is allowed to modify the data pointer in @info.
 *
 * Returns: a #GstPadProbeReturn
 */
typedef GstPadProbeReturn   (*GstPadProbeCallback)   (GstPad *pad, GstPadProbeInfo *info,
                                                      gpointer user_data);

/**
 * GstPadStickyEventsForeachFunction:
 * @pad: the #GstPad.
 * @event: (allow-none): a sticky #GstEvent.
 * @user_data: the #gpointer to optional user data.
 *
 * Callback used by gst_pad_sticky_events_foreach().
 *
 * When this function returns %TRUE, the next event will be
 * returned. When %FALSE is returned, gst_pad_sticky_events_foreach() will return.
 *
 * When @event is set to %NULL, the item will be removed from the list of sticky events.
 * @event can be replaced by assigning a new reference to it.
 * This function is responsible for unreffing the old event when
 * removing or modifying.
 *
 * Returns: %TRUE if the iteration should continue
 */
typedef gboolean  (*GstPadStickyEventsForeachFunction) (GstPad *pad, GstEvent **event,
                                                        gpointer user_data);

/**
 * GstPadFlags:
 * @GST_PAD_FLAG_BLOCKED: is dataflow on a pad blocked
 * @GST_PAD_FLAG_FLUSHING: is pad flushing
 * @GST_PAD_FLAG_EOS: is pad in EOS state
 * @GST_PAD_FLAG_BLOCKING: is pad currently blocking on a buffer or event
 * @GST_PAD_FLAG_NEED_PARENT: ensure that there is a parent object before calling
 *                       into the pad callbacks.
 * @GST_PAD_FLAG_NEED_RECONFIGURE: the pad should be reconfigured/renegotiated.
 *                            The flag has to be unset manually after
 *                            reconfiguration happened.
 * @GST_PAD_FLAG_PENDING_EVENTS: the pad has pending events
 * @GST_PAD_FLAG_FIXED_CAPS: the pad is using fixed caps. This means that
 *     once the caps are set on the pad, the default caps query function
 *     will only return those caps.
 * @GST_PAD_FLAG_PROXY_CAPS: the default event and query handler will forward
 *                      all events and queries to the internally linked pads
 *                      instead of discarding them.
 * @GST_PAD_FLAG_PROXY_ALLOCATION: the default query handler will forward
 *                      allocation queries to the internally linked pads
 *                      instead of discarding them.
 * @GST_PAD_FLAG_PROXY_SCHEDULING: the default query handler will forward
 *                      scheduling queries to the internally linked pads
 *                      instead of discarding them.
 * @GST_PAD_FLAG_ACCEPT_INTERSECT: the default accept-caps handler will check
 *                      it the caps intersect the query-caps result instead
 *                      of checking for a subset. This is interesting for
 *                      parsers that can accept incompletely specified caps.
 * @GST_PAD_FLAG_ACCEPT_TEMPLATE: the default accept-caps handler will use
 *                      the template pad caps instead of query caps to
 *                      compare with the accept caps. Use this in combination
 *                      with %GST_PAD_FLAG_ACCEPT_INTERSECT. (Since: 1.6)
 * @GST_PAD_FLAG_LAST: offset to define more flags
 *
 * Pad state flags
 */
typedef enum {
  GST_PAD_FLAG_BLOCKED          = (GST_OBJECT_FLAG_LAST << 0),
  GST_PAD_FLAG_FLUSHING         = (GST_OBJECT_FLAG_LAST << 1),
  GST_PAD_FLAG_EOS              = (GST_OBJECT_FLAG_LAST << 2),
  GST_PAD_FLAG_BLOCKING         = (GST_OBJECT_FLAG_LAST << 3),
  GST_PAD_FLAG_NEED_PARENT      = (GST_OBJECT_FLAG_LAST << 4),
  GST_PAD_FLAG_NEED_RECONFIGURE = (GST_OBJECT_FLAG_LAST << 5),
  GST_PAD_FLAG_PENDING_EVENTS   = (GST_OBJECT_FLAG_LAST << 6),
  GST_PAD_FLAG_FIXED_CAPS       = (GST_OBJECT_FLAG_LAST << 7),
  GST_PAD_FLAG_PROXY_CAPS       = (GST_OBJECT_FLAG_LAST << 8),
  GST_PAD_FLAG_PROXY_ALLOCATION = (GST_OBJECT_FLAG_LAST << 9),
  GST_PAD_FLAG_PROXY_SCHEDULING = (GST_OBJECT_FLAG_LAST << 10),
  GST_PAD_FLAG_ACCEPT_INTERSECT = (GST_OBJECT_FLAG_LAST << 11),
  GST_PAD_FLAG_ACCEPT_TEMPLATE  = (GST_OBJECT_FLAG_LAST << 12),
  /* padding */
  GST_PAD_FLAG_LAST        = (GST_OBJECT_FLAG_LAST << 16)
} GstPadFlags;

/**
 * GstPad:
 * @element_private: private data owned by the parent element
 * @padtemplate: padtemplate for this pad
 * @direction: the direction of the pad, cannot change after creating
 *             the pad.
 *
 * The #GstPad structure. Use the functions to update the variables.
 */
struct _GstPad {
  GstObject                      object;

  /*< public >*/
  gpointer                       element_private;

  GstPadTemplate                *padtemplate;

  GstPadDirection                direction;

  /*< private >*/
  /* streaming rec_lock */
  GRecMutex		         stream_rec_lock;
  GstTask			*task;

  /* block cond, mutex is from the object */
  GCond				 block_cond;
  GHookList                      probes;

  GstPadMode		         mode;
  GstPadActivateFunction	 activatefunc;
  gpointer                       activatedata;
  GDestroyNotify                 activatenotify;
  GstPadActivateModeFunction	 activatemodefunc;
  gpointer                       activatemodedata;
  GDestroyNotify                 activatemodenotify;

  /* pad link */
  GstPad			*peer;
  GstPadLinkFunction		 linkfunc;
  gpointer                       linkdata;
  GDestroyNotify                 linknotify;
  GstPadUnlinkFunction		 unlinkfunc;
  gpointer                       unlinkdata;
  GDestroyNotify                 unlinknotify;

  /* data transport functions */
  GstPadChainFunction		 chainfunc;
  gpointer                       chaindata;
  GDestroyNotify                 chainnotify;
  GstPadChainListFunction        chainlistfunc;
  gpointer                       chainlistdata;
  GDestroyNotify                 chainlistnotify;
  GstPadGetRangeFunction	 getrangefunc;
  gpointer                       getrangedata;
  GDestroyNotify                 getrangenotify;
  GstPadEventFunction		 eventfunc;
  gpointer                       eventdata;
  GDestroyNotify                 eventnotify;

  /* pad offset */
  gint64                         offset;

  /* generic query method */
  GstPadQueryFunction		 queryfunc;
  gpointer                       querydata;
  GDestroyNotify                 querynotify;

  /* internal links */
  GstPadIterIntLinkFunction      iterintlinkfunc;
  gpointer                       iterintlinkdata;
  GDestroyNotify                 iterintlinknotify;

  /* counts number of probes attached. */
  gint				 num_probes;
  gint				 num_blocked;

  GstPadPrivate                 *priv;

  union {
    gpointer _gst_reserved[GST_PADDING];
    struct {
      GstFlowReturn last_flowret;
      GstPadEventFullFunction eventfullfunc;
    } abi;
  } ABI;
};

struct _GstPadClass {
  GstObjectClass	parent_class;

  /* signal callbacks */
  void		(*linked)		(GstPad *pad, GstPad *peer);
  void		(*unlinked)		(GstPad *pad, GstPad *peer);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};


/***** helper macros *****/
/* GstPad */

/**
 * GST_PAD_NAME:
 * @pad: a #GstPad
 *
 * Get name of the given pad.
 * No locking is performed in this function, use gst_pad_get_name() instead.
 */
#define GST_PAD_NAME(pad)		(GST_OBJECT_NAME(pad))
/**
 * GST_PAD_PARENT:
 * @pad: a #GstPad
 *
 * Get the @pad parent.
 * No locking is performed in this function, use gst_pad_get_parent() instead.
 */
#define GST_PAD_PARENT(pad)		(GST_ELEMENT_CAST(GST_OBJECT_PARENT(pad)))
/**
 * GST_PAD_ELEMENT_PRIVATE:
 * @pad: a #GstPad
 *
 * Get the private data of @pad, which is usually some pad- or stream-specific
 * structure created by the element and set on the pad when creating it.
 * No locking is performed in this function.
 */
#define GST_PAD_ELEMENT_PRIVATE(pad)    (GST_PAD_CAST(pad)->element_private)
/**
 * GST_PAD_PAD_TEMPLATE:
 * @pad: a #GstPad
 *
 * Get the @pad #GstPadTemplate. It describes the possible media types
 * a @pad or an element factory can handle.
 */
#define GST_PAD_PAD_TEMPLATE(pad)	(GST_PAD_CAST(pad)->padtemplate)
/**
 * GST_PAD_DIRECTION:
 * @pad: a #GstPad
 *
 * Get the #GstPadDirection of the given @pad. Accessor macro, use
 * gst_pad_get_direction() instead.
 */
#define GST_PAD_DIRECTION(pad)		(GST_PAD_CAST(pad)->direction)
/**
 * GST_PAD_TASK:
 * @pad: a #GstPad
 *
 * Get the #GstTask of @pad. Accessor macro used by GStreamer. Use the
 * gst_pad_start_task(), gst_pad_stop_task() and gst_pad_pause_task()
 * functions instead.
 */
#define GST_PAD_TASK(pad)		(GST_PAD_CAST(pad)->task)
/**
 * GST_PAD_MODE:
 * @pad: a #GstPad
 *
 * Get the #GstPadMode of pad, which will be GST_PAD_MODE_NONE if the pad
 * has not been activated yet, and otherwise either GST_PAD_MODE_PUSH or
 * GST_PAD_MODE_PULL depending on which mode the pad was activated in.
 */
#define GST_PAD_MODE(pad)	        (GST_PAD_CAST(pad)->mode)
/**
 * GST_PAD_ACTIVATEFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadActivateFunction from @pad.
 */
#define GST_PAD_ACTIVATEFUNC(pad)	(GST_PAD_CAST(pad)->activatefunc)
/**
 * GST_PAD_ACTIVATEMODEFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadActivateModeFunction from the given @pad.
 */
#define GST_PAD_ACTIVATEMODEFUNC(pad)	(GST_PAD_CAST(pad)->activatemodefunc)
/**
 * GST_PAD_CHAINFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadChainFunction from the given @pad.
 */
#define GST_PAD_CHAINFUNC(pad)		(GST_PAD_CAST(pad)->chainfunc)
/**
 * GST_PAD_CHAINLISTFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadChainListFunction from the given @pad.
 */
#define GST_PAD_CHAINLISTFUNC(pad)      (GST_PAD_CAST(pad)->chainlistfunc)
/**
 * GST_PAD_GETRANGEFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadGetRangeFunction from the given @pad.
 */
#define GST_PAD_GETRANGEFUNC(pad)	(GST_PAD_CAST(pad)->getrangefunc)
/**
 * GST_PAD_EVENTFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadEventFunction from the given @pad, which
 * is the function that handles events on the pad. You can
 * use this to set your own event handling function on a pad
 * after you create it.  If your element derives from a base
 * class, use the base class's virtual functions instead.
 */
#define GST_PAD_EVENTFUNC(pad)		(GST_PAD_CAST(pad)->eventfunc)
/**
 * GST_PAD_EVENTFULLFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadEventFullFunction from the given @pad, which
 * is the function that handles events on the pad. You can
 * use this to set your own event handling function on a pad
 * after you create it.  If your element derives from a base
 * class, use the base class's virtual functions instead.
 *
 * Since: 1.8
 */
#define GST_PAD_EVENTFULLFUNC(pad)	(GST_PAD_CAST(pad)->ABI.abi.eventfullfunc)
/**
 * GST_PAD_QUERYFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadQueryFunction from @pad, which is the function
 * that handles queries on the pad. You can  use this to set your
 * own query handling function on a pad after you create it. If your
 * element derives from a base class, use the base class's virtual
 * functions instead.
 */
#define GST_PAD_QUERYFUNC(pad)		(GST_PAD_CAST(pad)->queryfunc)
/**
 * GST_PAD_ITERINTLINKFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadIterIntLinkFunction from the given @pad.
 */
#define GST_PAD_ITERINTLINKFUNC(pad)    (GST_PAD_CAST(pad)->iterintlinkfunc)
/**
 * GST_PAD_PEER:
 * @pad: a #GstPad
 *
 * Return the pad's peer member. This member is a pointer to the linked @pad.
 * No locking is performed in this function, use gst_pad_get_peer() instead.
 */
#define GST_PAD_PEER(pad)		(GST_PAD_CAST(pad)->peer)
/**
 * GST_PAD_LINKFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadLinkFunction for the given @pad.
 */
#define GST_PAD_LINKFUNC(pad)		(GST_PAD_CAST(pad)->linkfunc)
/**
 * GST_PAD_UNLINKFUNC:
 * @pad: a #GstPad
 *
 * Get the #GstPadUnlinkFunction from the given @pad.
 */
#define GST_PAD_UNLINKFUNC(pad)		(GST_PAD_CAST(pad)->unlinkfunc)
/**
 * GST_PAD_IS_SRC:
 * @pad: a #GstPad
 *
 * Returns: %TRUE if the pad is a source pad (i.e. produces data).
 */
#define GST_PAD_IS_SRC(pad)		(GST_PAD_DIRECTION(pad) == GST_PAD_SRC)
/**
 * GST_PAD_IS_SINK:
 * @pad: a #GstPad
 *
 * Returns: %TRUE if the pad is a sink pad (i.e. consumes data).
 */
#define GST_PAD_IS_SINK(pad)		(GST_PAD_DIRECTION(pad) == GST_PAD_SINK)
/**
 * GST_PAD_IS_LINKED:
 * @pad: a #GstPad
 *
 * Returns: %TRUE if the pad is linked to another pad. Use gst_pad_is_linked()
 * instead.
 */
#define GST_PAD_IS_LINKED(pad)		(GST_PAD_PEER(pad) != NULL)
/**
 * GST_PAD_IS_ACTIVE:
 * @pad: a #GstPad
 *
 * Returns: %TRUE if the pad has been activated.
 */
#define GST_PAD_IS_ACTIVE(pad)          (GST_PAD_MODE(pad) != GST_PAD_MODE_NONE)
/**
 * GST_PAD_IS_BLOCKED:
 * @pad: a #GstPad
 *
 * Check if the dataflow on a @pad is blocked. Use gst_pad_is_blocked() instead.
 */
#define GST_PAD_IS_BLOCKED(pad)		(GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_BLOCKED))
/**
 * GST_PAD_IS_BLOCKING:
 * @pad: a #GstPad
 *
 * Check if the @pad is currently blocking on a buffer or event. Use
 * gst_pad_is_blocking() instead.
 */
#define GST_PAD_IS_BLOCKING(pad)	(GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_BLOCKING))
/**
 * GST_PAD_IS_FLUSHING:
 * @pad: a #GstPad
 *
 * Check if the given @pad is flushing.
 */
#define GST_PAD_IS_FLUSHING(pad)	(GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_FLUSHING))
/**
 * GST_PAD_SET_FLUSHING:
 * @pad: a #GstPad
 *
 * Set the given @pad to flushing state, which means it will not accept any
 * more events, queries or buffers, and return GST_FLOW_FLUSHING if any buffers
 * are pushed on it. This usually happens when the pad is shut down or when
 * a flushing seek happens. This is used inside GStreamer when flush start/stop
 * events pass through pads, or when an element state is changed and pads are
 * activated or deactivated.
 */
#define GST_PAD_SET_FLUSHING(pad)	(GST_OBJECT_FLAG_SET (pad, GST_PAD_FLAG_FLUSHING))
/**
 * GST_PAD_UNSET_FLUSHING:
 * @pad: a #GstPad
 *
 * Unset the flushing flag.
 */
#define GST_PAD_UNSET_FLUSHING(pad)	(GST_OBJECT_FLAG_UNSET (pad, GST_PAD_FLAG_FLUSHING))
/**
 * GST_PAD_IS_EOS:
 * @pad: a #GstPad
 *
 * Check if the @pad is in EOS state.
 */
#define GST_PAD_IS_EOS(pad)	        (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_EOS))
/**
 * GST_PAD_NEEDS_RECONFIGURE:
 * @pad: a #GstPad
 *
 * Check if the @pad should be reconfigured/renegotiated.
 * The flag has to be unset manually after reconfiguration happened.
 * Use gst_pad_needs_reconfigure() or gst_pad_check_reconfigure() instead.
 */
#define GST_PAD_NEEDS_RECONFIGURE(pad)  (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_NEED_RECONFIGURE))
/**
 * GST_PAD_HAS_PENDING_EVENTS:
 * @pad: a #GstPad
 *
 * Check if the given @pad has pending events. This is used internally by
 * GStreamer.
 */
#define GST_PAD_HAS_PENDING_EVENTS(pad) (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_PENDING_EVENTS))
/**
 * GST_PAD_IS_FIXED_CAPS:
 * @pad: a #GstPad
 *
 * Check if the given @pad is using fixed caps, which means that
 * once the caps are set on the @pad, the caps query function will
 * only return those caps. See gst_pad_use_fixed_caps().
 */
#define GST_PAD_IS_FIXED_CAPS(pad)	(GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_FIXED_CAPS))
/**
 * GST_PAD_NEEDS_PARENT:
 * @pad: a #GstPad
 *
 * Check if there is a parent object before calling into the @pad callbacks.
 * This is used internally by GStreamer.
 */
#define GST_PAD_NEEDS_PARENT(pad)       (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_NEED_PARENT))
/**
 * GST_PAD_IS_PROXY_CAPS:
 * @pad: a #GstPad
 *
 * Check if the given @pad is set to proxy caps. This means that the default
 * event and query handler will forward all events and queries to the
 * internally linked @pads instead of discarding them.
 */
#define GST_PAD_IS_PROXY_CAPS(pad)      (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_PROXY_CAPS))
/**
 * GST_PAD_SET_PROXY_CAPS:
 * @pad: a #GstPad
 *
 * Set @pad to proxy caps, so that all caps-related events and queries are
 * proxied down- or upstream to the other side of the element automatically.
 * Set this if the element always outputs data in the exact same format as it
 * receives as input. This is just for convenience to avoid implementing some
 * standard event and query handling code in an element.
 */
#define GST_PAD_SET_PROXY_CAPS(pad)     (GST_OBJECT_FLAG_SET (pad, GST_PAD_FLAG_PROXY_CAPS))
/**
 * GST_PAD_UNSET_PROXY_CAPS:
 * @pad: a #GstPad
 *
 * Unset proxy caps flag.
 */
#define GST_PAD_UNSET_PROXY_CAPS(pad)   (GST_OBJECT_FLAG_UNSET (pad, GST_PAD_FLAG_PROXY_CAPS))
/**
 * GST_PAD_IS_PROXY_ALLOCATION:
 * @pad: a #GstPad
 *
 * Check if the given @pad is set as proxy allocation which means
 * that the default query handler will forward allocation queries to the
 * internally linked @pads instead of discarding them.
 */
#define GST_PAD_IS_PROXY_ALLOCATION(pad)    (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_PROXY_ALLOCATION))
/**
 * GST_PAD_SET_PROXY_ALLOCATION:
 * @pad: a #GstPad
 *
 * Set @pad to proxy allocation queries, which means that the default query
 * handler will forward allocation queries to the internally linked @pads
 * instead of discarding them.
 * Set this if the element always outputs data in the exact same format as it
 * receives as input. This is just for convenience to avoid implementing some
 * standard query handling code in an element.
 */
#define GST_PAD_SET_PROXY_ALLOCATION(pad)   (GST_OBJECT_FLAG_SET (pad, GST_PAD_FLAG_PROXY_ALLOCATION))
/**
 * GST_PAD_UNSET_PROXY_ALLOCATION:
 * @pad: a #GstPad
 *
 * Unset proxy allocation flag.
 */
#define GST_PAD_UNSET_PROXY_ALLOCATION(pad) (GST_OBJECT_FLAG_UNSET (pad, GST_PAD_FLAG_PROXY_ALLOCATION))
/**
 * GST_PAD_IS_PROXY_SCHEDULING:
 * @pad: a #GstPad
 *
 * Check if the given @pad is set to proxy scheduling queries, which means that
 * the default query handler will forward scheduling queries to the internally
 * linked @pads instead of discarding them.
 */
#define GST_PAD_IS_PROXY_SCHEDULING(pad)    (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_PROXY_SCHEDULING))
/**
 * GST_PAD_SET_PROXY_SCHEDULING:
 * @pad: a #GstPad
 *
 * Set @pad to proxy scheduling queries, which means that the default query
 * handler will forward scheduling queries to the internally linked @pads
 * instead of discarding them. You will usually want to handle scheduling
 * queries explicitly if your element supports multiple scheduling modes.
 */
#define GST_PAD_SET_PROXY_SCHEDULING(pad)   (GST_OBJECT_FLAG_SET (pad, GST_PAD_FLAG_PROXY_SCHEDULING))
/**
 * GST_PAD_UNSET_PROXY_SCHEDULING:
 * @pad: a #GstPad
 *
 * Unset proxy scheduling flag.
 */
#define GST_PAD_UNSET_PROXY_SCHEDULING(pad) (GST_OBJECT_FLAG_UNSET (pad, GST_PAD_FLAG_PROXY_SCHEDULING))
/**
 * GST_PAD_IS_ACCEPT_INTERSECT:
 * @pad: a #GstPad
 *
 * Check if the pad's accept intersect flag is set. The default accept-caps
 * handler will check if the caps intersect the query-caps result instead of
 * checking for a subset. This is interesting for parser elements that can
 * accept incompletely specified caps.
 */
#define GST_PAD_IS_ACCEPT_INTERSECT(pad)    (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_ACCEPT_INTERSECT))
/**
 * GST_PAD_SET_ACCEPT_INTERSECT:
 * @pad: a #GstPad
 *
 * Set @pad to by default accept caps by intersecting the result instead of
 * checking for a subset. This is interesting for parser elements that can
 * accept incompletely specified caps.
 */
#define GST_PAD_SET_ACCEPT_INTERSECT(pad)   (GST_OBJECT_FLAG_SET (pad, GST_PAD_FLAG_ACCEPT_INTERSECT))
/**
 * GST_PAD_UNSET_ACCEPT_INTERSECT:
 * @pad: a #GstPad
 *
 * Unset accept intersect flag.
 */
#define GST_PAD_UNSET_ACCEPT_INTERSECT(pad) (GST_OBJECT_FLAG_UNSET (pad, GST_PAD_FLAG_ACCEPT_INTERSECT))
/**
 * GST_PAD_IS_ACCEPT_TEMPLATE:
 * @pad: a #GstPad
 *
 * Check if the pad's accept caps operation will use the pad template caps.
 * The default accept-caps will do a query caps to get the caps, which might
 * be querying downstream causing unnecessary overhead. It is recommended to
 * implement a proper accept-caps query handler or to use this flag to prevent
 * recursive accept-caps handling.
 *
 * Since: 1.6
 */
#define GST_PAD_IS_ACCEPT_TEMPLATE(pad)    (GST_OBJECT_FLAG_IS_SET (pad, GST_PAD_FLAG_ACCEPT_TEMPLATE))
/**
 * GST_PAD_SET_ACCEPT_TEMPLATE:
 * @pad: a #GstPad
 *
 * Set @pad to by default use the pad template caps to compare with
 * the accept caps instead of using a caps query result.
 *
 * Since: 1.6
 */
#define GST_PAD_SET_ACCEPT_TEMPLATE(pad)   (GST_OBJECT_FLAG_SET (pad, GST_PAD_FLAG_ACCEPT_TEMPLATE))
/**
 * GST_PAD_UNSET_ACCEPT_TEMPLATE:
 * @pad: a #GstPad
 *
 * Unset accept template flag.
 *
 * Since: 1.6
 */
#define GST_PAD_UNSET_ACCEPT_TEMPLATE(pad) (GST_OBJECT_FLAG_UNSET (pad, GST_PAD_FLAG_ACCEPT_TEMPLATE))
/**
 * GST_PAD_GET_STREAM_LOCK:
 * @pad: a #GstPad
 *
 * Get the stream lock of @pad. The stream lock is protecting the
 * resources used in the data processing functions of @pad. Accessor
 * macro, use GST_PAD_STREAM_LOCK() and GST_PAD_STREAM_UNLOCK() instead
 * to take/release the pad's stream lock.
 */
#define GST_PAD_GET_STREAM_LOCK(pad)    (&(GST_PAD_CAST(pad)->stream_rec_lock))
/**
 * GST_PAD_STREAM_LOCK:
 * @pad: a #GstPad
 *
 * Take the pad's stream lock. The stream lock is recursive and will be taken
 * when buffers or serialized downstream events are pushed on a pad.
 */
#define GST_PAD_STREAM_LOCK(pad)        g_rec_mutex_lock(GST_PAD_GET_STREAM_LOCK(pad))
/**
 * GST_PAD_STREAM_TRYLOCK:
 * @pad: a #GstPad
 *
 * Try to take the pad's stream lock, and return %TRUE if the lock could be
 * taken, and otherwise %FALSE.
 */
#define GST_PAD_STREAM_TRYLOCK(pad)     g_rec_mutex_trylock(GST_PAD_GET_STREAM_LOCK(pad))
/**
 * GST_PAD_STREAM_UNLOCK:
 * @pad: a #GstPad
 *
 * Release the pad's stream lock.
 */
#define GST_PAD_STREAM_UNLOCK(pad)      g_rec_mutex_unlock(GST_PAD_GET_STREAM_LOCK(pad))
/**
 * GST_PAD_LAST_FLOW_RETURN:
 * @pad: a #GstPad
 *
 * Gets the last flow return on this pad
 *
 * Since: 1.4
 */
#define GST_PAD_LAST_FLOW_RETURN(pad)   (GST_PAD_CAST(pad)->ABI.abi.last_flowret)

#define GST_PAD_BLOCK_GET_COND(pad)     (&GST_PAD_CAST(pad)->block_cond)
#define GST_PAD_BLOCK_WAIT(pad)         (g_cond_wait(GST_PAD_BLOCK_GET_COND (pad), GST_OBJECT_GET_LOCK (pad)))
#define GST_PAD_BLOCK_SIGNAL(pad)       (g_cond_signal(GST_PAD_BLOCK_GET_COND (pad)))
#define GST_PAD_BLOCK_BROADCAST(pad)    (g_cond_broadcast(GST_PAD_BLOCK_GET_COND (pad)))

GST_API
GType			gst_pad_get_type			(void);

/* creating pads */

GST_API
GstPad*			gst_pad_new				(const gchar *name, GstPadDirection direction);

GST_API
GstPad*			gst_pad_new_from_template		(GstPadTemplate *templ, const gchar *name);

GST_API
GstPad*			gst_pad_new_from_static_template	(GstStaticPadTemplate *templ, const gchar *name);


/**
 * gst_pad_get_name:
 * @pad: the pad to get the name from
 *
 * Get a copy of the name of the pad. g_free() after usage.
 *
 * MT safe.
 */
#define gst_pad_get_name(pad) gst_object_get_name (GST_OBJECT_CAST (pad))
/**
 * gst_pad_get_parent:
 * @pad: the pad to get the parent of
 *
 * Get the parent of @pad. This function increases the refcount
 * of the parent object so you should gst_object_unref() it after usage.
 * Can return %NULL if the pad did not have a parent.
 *
 * MT safe.
 *
 * Returns: (nullable): the parent
 */
#define gst_pad_get_parent(pad) gst_object_get_parent (GST_OBJECT_CAST (pad))

GST_API
GstPadDirection		gst_pad_get_direction			(GstPad *pad);

GST_API
gboolean		gst_pad_set_active			(GstPad *pad, gboolean active);

GST_API
gboolean		gst_pad_is_active			(GstPad *pad);

GST_API
gboolean		gst_pad_activate_mode			(GstPad *pad, GstPadMode mode,
                                                                 gboolean active);
GST_API
gulong                  gst_pad_add_probe                       (GstPad *pad,
								 GstPadProbeType mask,
								 GstPadProbeCallback callback,
                                                                 gpointer user_data,
                                                                 GDestroyNotify destroy_data);
GST_API
void                    gst_pad_remove_probe                    (GstPad *pad, gulong id);

GST_API
gboolean		gst_pad_is_blocked			(GstPad *pad);

GST_API
gboolean		gst_pad_is_blocking			(GstPad *pad);

GST_API
void                    gst_pad_mark_reconfigure                (GstPad *pad);

GST_API
gboolean		gst_pad_needs_reconfigure               (GstPad *pad);

GST_API
gboolean		gst_pad_check_reconfigure               (GstPad *pad);

GST_API
void			gst_pad_set_element_private		(GstPad *pad, gpointer priv);

GST_API
gpointer		gst_pad_get_element_private		(GstPad *pad);

GST_API
GstPadTemplate*		gst_pad_get_pad_template		(GstPad *pad);

GST_API
GstFlowReturn           gst_pad_store_sticky_event              (GstPad *pad, GstEvent *event);

GST_API
GstEvent*               gst_pad_get_sticky_event                (GstPad *pad, GstEventType event_type,
                                                                 guint idx);

GST_API
void                    gst_pad_sticky_events_foreach           (GstPad *pad, GstPadStickyEventsForeachFunction foreach_func, gpointer user_data);

/* data passing setup functions */

GST_API
void			gst_pad_set_activate_function_full	(GstPad *pad,
                                                                 GstPadActivateFunction activate,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
void			gst_pad_set_activatemode_function_full	(GstPad *pad,
                                                                 GstPadActivateModeFunction activatemode,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
/* data passing functions */

GST_API
void			gst_pad_set_chain_function_full		(GstPad *pad,
                                                                 GstPadChainFunction chain,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
void			gst_pad_set_chain_list_function_full	(GstPad *pad,
                                                                 GstPadChainListFunction chainlist,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
void			gst_pad_set_getrange_function_full	(GstPad *pad,
                                                                 GstPadGetRangeFunction get,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
void			gst_pad_set_event_function_full		(GstPad *pad,
                                                                 GstPadEventFunction event,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
void			gst_pad_set_event_full_function_full	(GstPad *pad,
                                                                 GstPadEventFullFunction event,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);

#define gst_pad_set_activate_function(p,f)      gst_pad_set_activate_function_full((p),(f),NULL,NULL)
#define gst_pad_set_activatemode_function(p,f)  gst_pad_set_activatemode_function_full((p),(f),NULL,NULL)
#define gst_pad_set_chain_function(p,f)         gst_pad_set_chain_function_full((p),(f),NULL,NULL)
#define gst_pad_set_chain_list_function(p,f)    gst_pad_set_chain_list_function_full((p),(f),NULL,NULL)
#define gst_pad_set_getrange_function(p,f)      gst_pad_set_getrange_function_full((p),(f),NULL,NULL)
#define gst_pad_set_event_function(p,f)         gst_pad_set_event_function_full((p),(f),NULL,NULL)
#define gst_pad_set_event_full_function(p,f)    gst_pad_set_event_full_function_full((p),(f),NULL,NULL)

/* pad links */

GST_API
void			gst_pad_set_link_function_full		(GstPad *pad,
                                                                 GstPadLinkFunction link,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
void			gst_pad_set_unlink_function_full        (GstPad *pad,
                                                                 GstPadUnlinkFunction unlink,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);

#define gst_pad_set_link_function(p,f)          gst_pad_set_link_function_full((p),(f),NULL,NULL)
#define gst_pad_set_unlink_function(p,f)        gst_pad_set_unlink_function_full((p),(f),NULL,NULL)

GST_API
gboolean                gst_pad_can_link                        (GstPad *srcpad, GstPad *sinkpad);

GST_API
GstPadLinkReturn        gst_pad_link				(GstPad *srcpad, GstPad *sinkpad);

GST_API
GstPadLinkReturn        gst_pad_link_full			(GstPad *srcpad, GstPad *sinkpad, GstPadLinkCheck flags);

GST_API
gboolean		gst_pad_unlink				(GstPad *srcpad, GstPad *sinkpad);

GST_API
gboolean		gst_pad_is_linked			(GstPad *pad);

GST_API
GstPad*			gst_pad_get_peer			(GstPad *pad);

GST_API
GstCaps*                gst_pad_get_pad_template_caps		(GstPad *pad);

/* capsnego function for linked/unlinked pads */

GST_API
GstCaps *		gst_pad_get_current_caps                (GstPad * pad);

GST_API
gboolean		gst_pad_has_current_caps                (GstPad * pad);

/* capsnego for linked pads */

GST_API
GstCaps *		gst_pad_get_allowed_caps		(GstPad * pad);

/* pad offsets */

GST_API
gint64                  gst_pad_get_offset                      (GstPad *pad);

GST_API
void                    gst_pad_set_offset                      (GstPad *pad, gint64 offset);

/* data passing functions to peer */

GST_API
GstFlowReturn		gst_pad_push				(GstPad *pad, GstBuffer *buffer);

GST_API
GstFlowReturn		gst_pad_push_list			(GstPad *pad, GstBufferList *list);

GST_API
GstFlowReturn		gst_pad_pull_range			(GstPad *pad, guint64 offset, guint size,
								 GstBuffer **buffer);
GST_API
gboolean		gst_pad_push_event			(GstPad *pad, GstEvent *event);

GST_API
gboolean		gst_pad_event_default			(GstPad *pad, GstObject *parent,
                                                                 GstEvent *event);
GST_API
GstFlowReturn           gst_pad_get_last_flow_return            (GstPad *pad);

/* data passing functions on pad */

GST_API
GstFlowReturn		gst_pad_chain				(GstPad *pad, GstBuffer *buffer);

GST_API
GstFlowReturn		gst_pad_chain_list                      (GstPad *pad, GstBufferList *list);

GST_API
GstFlowReturn		gst_pad_get_range			(GstPad *pad, guint64 offset, guint size,
								 GstBuffer **buffer);
GST_API
gboolean		gst_pad_send_event			(GstPad *pad, GstEvent *event);

/* pad tasks */

GST_API
gboolean		gst_pad_start_task			(GstPad *pad, GstTaskFunction func,
								 gpointer user_data, GDestroyNotify notify);
GST_API
gboolean		gst_pad_pause_task			(GstPad *pad);

GST_API
gboolean		gst_pad_stop_task			(GstPad *pad);

GST_API
GstTaskState	gst_pad_get_task_state		(GstPad *pad);

/* internal links */

GST_API
void                    gst_pad_set_iterate_internal_links_function_full (GstPad * pad,
                                                                 GstPadIterIntLinkFunction iterintlink,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);

GST_API
GstIterator *           gst_pad_iterate_internal_links          (GstPad * pad);

GST_API
GstIterator *           gst_pad_iterate_internal_links_default  (GstPad * pad, GstObject *parent);

#define gst_pad_set_iterate_internal_links_function(p,f) gst_pad_set_iterate_internal_links_function_full((p),(f),NULL,NULL)

GST_API
GstPad *                gst_pad_get_single_internal_link        (GstPad * pad);

/* generic query function */

GST_API
gboolean		gst_pad_query				(GstPad *pad, GstQuery *query);

GST_API
gboolean		gst_pad_peer_query			(GstPad *pad, GstQuery *query);

GST_API
void			gst_pad_set_query_function_full		(GstPad *pad, GstPadQueryFunction query,
                                                                 gpointer user_data,
                                                                 GDestroyNotify notify);
GST_API
gboolean		gst_pad_query_default			(GstPad *pad, GstObject *parent,
                                                                 GstQuery *query);

#define gst_pad_set_query_function(p,f)   gst_pad_set_query_function_full((p),(f),NULL,NULL)

/* misc helper functions */

GST_API
gboolean		gst_pad_forward                         (GstPad *pad, GstPadForwardFunction forward,
								 gpointer user_data);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstPad, gst_object_unref)

G_END_DECLS

#endif /* __GST_PAD_H__ */
