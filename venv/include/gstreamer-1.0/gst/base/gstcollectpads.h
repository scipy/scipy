/* GStreamer
 * Copyright (C) 2005 Wim Taymans <wim@fluendo.com>
 * Copyright (C) 2008 Mark Nauwelaerts <mnauw@users.sourceforge.net>
 *
 * gstcollectpads.h:
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

#ifndef __GST_COLLECT_PADS_H__
#define __GST_COLLECT_PADS_H__

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_COLLECT_PADS            (gst_collect_pads_get_type())
#define GST_COLLECT_PADS(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_COLLECT_PADS,GstCollectPads))
#define GST_COLLECT_PADS_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_COLLECT_PADS,GstCollectPadsClass))
#define GST_COLLECT_PADS_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_COLLECT_PADS,GstCollectPadsClass))
#define GST_IS_COLLECT_PADS(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_COLLECT_PADS))
#define GST_IS_COLLECT_PADS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_COLLECT_PADS))

typedef struct _GstCollectData GstCollectData;
typedef struct _GstCollectDataPrivate GstCollectDataPrivate;
typedef struct _GstCollectPads GstCollectPads;
typedef struct _GstCollectPadsPrivate GstCollectPadsPrivate;
typedef struct _GstCollectPadsClass GstCollectPadsClass;

/**
 * GstCollectDataDestroyNotify:
 * @data: the #GstCollectData that will be freed
 *
 * A function that will be called when the #GstCollectData will be freed.
 * It is passed the pointer to the structure and should free any custom
 * memory and resources allocated for it.
 */
typedef void (*GstCollectDataDestroyNotify) (GstCollectData *data);

/**
 * GstCollectPadsStateFlags:
 * @GST_COLLECT_PADS_STATE_EOS:         Set if collectdata's pad is EOS.
 * @GST_COLLECT_PADS_STATE_FLUSHING:    Set if collectdata's pad is flushing.
 * @GST_COLLECT_PADS_STATE_NEW_SEGMENT: Set if collectdata's pad received a
 *                                      new_segment event.
 * @GST_COLLECT_PADS_STATE_WAITING:     Set if collectdata's pad must be waited
 *                                      for when collecting.
 * @GST_COLLECT_PADS_STATE_LOCKED:      Set collectdata's pad WAITING state must
 *                                      not be changed.
 * #GstCollectPadsStateFlags indicate private state of a collectdata('s pad).
 */
typedef enum {
  GST_COLLECT_PADS_STATE_EOS = 1 << 0,
  GST_COLLECT_PADS_STATE_FLUSHING = 1 << 1,
  GST_COLLECT_PADS_STATE_NEW_SEGMENT = 1 << 2,
  GST_COLLECT_PADS_STATE_WAITING = 1 << 3,
  GST_COLLECT_PADS_STATE_LOCKED = 1 << 4
} GstCollectPadsStateFlags;

/**
 * GST_COLLECT_PADS_STATE:
 * @data: a #GstCollectData.
 *
 * A flags word containing #GstCollectPadsStateFlags flags set
 * on this collected pad.
 */
#define GST_COLLECT_PADS_STATE(data)                 (((GstCollectData *) data)->state)
/**
 * GST_COLLECT_PADS_STATE_IS_SET:
 * @data: a #GstCollectData.
 * @flag: the #GstCollectPadsStateFlags to check.
 *
 * Gives the status of a specific flag on a collected pad.
 */
#define GST_COLLECT_PADS_STATE_IS_SET(data,flag)     !!(GST_COLLECT_PADS_STATE (data) & flag)
/**
 * GST_COLLECT_PADS_STATE_SET:
 * @data: a #GstCollectData.
 * @flag: the #GstCollectPadsStateFlags to set.
 *
 * Sets a state flag on a collected pad.
 */
#define GST_COLLECT_PADS_STATE_SET(data,flag)        (GST_COLLECT_PADS_STATE (data) |= flag)
/**
 * GST_COLLECT_PADS_STATE_UNSET:
 * @data: a #GstCollectData.
 * @flag: the #GstCollectPadsStateFlags to clear.
 *
 * Clears a state flag on a collected pad.
 */
#define GST_COLLECT_PADS_STATE_UNSET(data,flag)      (GST_COLLECT_PADS_STATE (data) &= ~(flag))

/**
 * GST_COLLECT_PADS_DTS:
 * @data: A #GstCollectData.
 *
 * Returns the DTS that has been converted to running time when using
 * gst_collect_pads_clip_running_time(). Unlike the value saved into
 * the buffer, this value is of type gint64 and may be negative. This allow
 * properly handling streams with frame reordering where the first DTS may
 * be negative. If the initial DTS was not set, this value will be
 * set to %G_MININT64.
 *
 * Since: 1.6
 */
#define GST_COLLECT_PADS_DTS(data)                   (((GstCollectData *) data)->ABI.abi.dts)

/**
 * GST_COLLECT_PADS_DTS_IS_VALID:
 * @data: A #GstCollectData.
 *
 * Check if running DTS value store is valid.
 *
 * Since: 1.6
 */
#define GST_COLLECT_PADS_DTS_IS_VALID(data)          (GST_CLOCK_STIME_IS_VALID (GST_COLLECT_PADS_DTS (data)))

/**
 * GstCollectData:
 * @collect: owner #GstCollectPads
 * @pad: #GstPad managed by this data
 * @buffer: currently queued buffer.
 * @pos: position in the buffer
 * @segment: last segment received.
 * @dts: the signed version of the DTS converted to running time. To access
 *       this member, use %GST_COLLECT_PADS_DTS macro. (Since: 1.6)
 *
 * Structure used by the collect_pads.
 */
struct _GstCollectData
{
  /* with STREAM_LOCK of @collect */
  GstCollectPads        *collect;
  GstPad                *pad;
  GstBuffer             *buffer;
  guint                  pos;
  GstSegment             segment;

  /*< private >*/
  /* state: bitfield for easier extension;
   * eos, flushing, new_segment, waiting */
  GstCollectPadsStateFlags    state;

  GstCollectDataPrivate *priv;

  union {
    struct {
      /*< public >*/
      gint64 dts;
      /*< private >*/
    } abi;
    gpointer _gst_reserved[GST_PADDING];
  } ABI;
};

/**
 * GstCollectPadsFunction:
 * @pads: the #GstCollectPads that triggered the callback
 * @user_data: user data passed to gst_collect_pads_set_function()
 *
 * A function that will be called when all pads have received data.
 *
 * Returns: %GST_FLOW_OK for success
 */
typedef GstFlowReturn (*GstCollectPadsFunction) (GstCollectPads *pads, gpointer user_data);

/**
 * GstCollectPadsBufferFunction:
 * @pads: the #GstCollectPads that triggered the callback
 * @data: the #GstCollectData of pad that has received the buffer
 * @buffer: (transfer full): the #GstBuffer
 * @user_data: user data passed to gst_collect_pads_set_buffer_function()
 *
 * A function that will be called when a (considered oldest) buffer can be muxed.
 * If all pads have reached EOS, this function is called with %NULL @buffer
 * and %NULL @data.
 *
 * Returns: %GST_FLOW_OK for success
 */
typedef GstFlowReturn (*GstCollectPadsBufferFunction) (GstCollectPads *pads, GstCollectData *data,
                                                       GstBuffer *buffer, gpointer user_data);

/**
 * GstCollectPadsCompareFunction:
 * @pads: the #GstCollectPads that is comparing the timestamps
 * @data1: the first #GstCollectData
 * @timestamp1: the first timestamp
 * @data2: the second #GstCollectData
 * @timestamp2: the second timestamp
 * @user_data: user data passed to gst_collect_pads_set_compare_function()
 *
 * A function for comparing two timestamps of buffers or newsegments collected on one pad.
 *
 * Returns: Integer less than zero when first timestamp is deemed older than the second one.
 *          Zero if the timestamps are deemed equally old.
 *          Integer greater than zero when second timestamp is deemed older than the first one.
 */
typedef gint (*GstCollectPadsCompareFunction) (GstCollectPads *pads,
                                               GstCollectData * data1, GstClockTime timestamp1,
                                               GstCollectData * data2, GstClockTime timestamp2,
                                               gpointer user_data);

/**
 * GstCollectPadsEventFunction:
 * @pads: the #GstCollectPads that triggered the callback
 * @pad: the #GstPad that received an event
 * @event: the #GstEvent received
 * @user_data: user data passed to gst_collect_pads_set_event_function()
 *
 * A function that will be called while processing an event. It takes
 * ownership of the event and is responsible for chaining up (to
 * gst_collect_pads_event_default()) or dropping events (such typical cases
 * being handled by the default handler).
 *
 * Returns: %TRUE if the pad could handle the event
 */
typedef gboolean (*GstCollectPadsEventFunction)        (GstCollectPads *pads, GstCollectData * pad,
                                                        GstEvent * event, gpointer user_data);


/**
 * GstCollectPadsQueryFunction:
 * @pads: the #GstCollectPads that triggered the callback
 * @pad: the #GstPad that received an event
 * @query: the #GstEvent received
 * @user_data: user data passed to gst_collect_pads_set_query_function()
 *
 * A function that will be called while processing a query. It takes
 * ownership of the query and is responsible for chaining up (to
 * events downstream (with gst_pad_event_default()).
 *
 * Returns: %TRUE if the pad could handle the event
 */
typedef gboolean (*GstCollectPadsQueryFunction)        (GstCollectPads *pads, GstCollectData * pad,
                                                        GstQuery * query, gpointer user_data);

/**
 * GstCollectPadsClipFunction:
 * @pads: a #GstCollectPads
 * @data: a #GstCollectData
 * @inbuffer: (transfer full): the input #GstBuffer
 * @outbuffer: (out): the output #GstBuffer
 * @user_data: user data
 *
 * A function that will be called when @inbuffer is received on the pad managed
 * by @data in the collectpad object @pads.
 *
 * The function should use the segment of @data and the negotiated media type on
 * the pad to perform clipping of @inbuffer.
 *
 * This function takes ownership of @inbuffer and should output a buffer in
 * @outbuffer or return %NULL in @outbuffer if the buffer should be dropped.
 *
 * Returns: a #GstFlowReturn that corresponds to the result of clipping.
 */
typedef GstFlowReturn (*GstCollectPadsClipFunction) (GstCollectPads *pads, GstCollectData *data,
                                                     GstBuffer *inbuffer, GstBuffer **outbuffer,
                                                     gpointer user_data);


/**
 * GstCollectPadsFlushFunction:
 * @pads: a #GstCollectPads
 * @user_data: user data
 *
 * A function that will be called while processing a flushing seek event.
 *
 * The function should flush any internal state of the element and the state of
 * all the pads. It should clear only the state not directly managed by the
 * @pads object. It is therefore not necessary to call
 * gst_collect_pads_set_flushing nor gst_collect_pads_clear from this function.
 *
 * Since: 1.4
 */
typedef void (*GstCollectPadsFlushFunction) (GstCollectPads *pads, gpointer user_data);

/**
 * GST_COLLECT_PADS_GET_STREAM_LOCK:
 * @pads: a #GstCollectPads
 *
 * Get the stream lock of @pads. The stream lock is used to coordinate and
 * serialize execution among the various streams being collected, and in
 * protecting the resources used to accomplish this.
 */
#define GST_COLLECT_PADS_GET_STREAM_LOCK(pads) (&((GstCollectPads *)pads)->stream_lock)
/**
 * GST_COLLECT_PADS_STREAM_LOCK:
 * @pads: a #GstCollectPads
 *
 * Lock the stream lock of @pads.
 */
#define GST_COLLECT_PADS_STREAM_LOCK(pads)     g_rec_mutex_lock(GST_COLLECT_PADS_GET_STREAM_LOCK (pads))
/**
 * GST_COLLECT_PADS_STREAM_UNLOCK:
 * @pads: a #GstCollectPads
 *
 * Unlock the stream lock of @pads.
 */
#define GST_COLLECT_PADS_STREAM_UNLOCK(pads)   g_rec_mutex_unlock(GST_COLLECT_PADS_GET_STREAM_LOCK (pads))

/**
 * GstCollectPads:
 * @data: (element-type GstBase.CollectData): #GList of #GstCollectData managed
 *   by this #GstCollectPads.
 *
 * Collectpads object.
 */
struct _GstCollectPads {
  GstObject      object;

  /*< public >*/ /* with LOCK and/or STREAM_LOCK */
  GSList        *data;                  /* list of CollectData items */

  /*< private >*/
  GRecMutex      stream_lock;          /* used to serialize collection among several streams */

  GstCollectPadsPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

struct _GstCollectPadsClass {
  GstObjectClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_BASE_API
GType           gst_collect_pads_get_type (void);

/* creating the object */

GST_BASE_API
GstCollectPads* gst_collect_pads_new      (void);

/* set the callbacks */

GST_BASE_API
void            gst_collect_pads_set_function         (GstCollectPads *pads,
                                                       GstCollectPadsFunction func,
                                                       gpointer user_data);
GST_BASE_API
void            gst_collect_pads_set_buffer_function  (GstCollectPads *pads,
                                                       GstCollectPadsBufferFunction func,
                                                       gpointer user_data);
GST_BASE_API
void            gst_collect_pads_set_event_function   (GstCollectPads *pads,
                                                       GstCollectPadsEventFunction func,
                                                       gpointer user_data);
GST_BASE_API
void            gst_collect_pads_set_query_function   (GstCollectPads *pads,
                                                       GstCollectPadsQueryFunction func,
                                                       gpointer user_data);
GST_BASE_API
void            gst_collect_pads_set_compare_function (GstCollectPads *pads,
                                                       GstCollectPadsCompareFunction func,
                                                       gpointer user_data);
GST_BASE_API
void            gst_collect_pads_set_clip_function    (GstCollectPads *pads,
                                                       GstCollectPadsClipFunction clipfunc,
                                                       gpointer user_data);
GST_BASE_API
void            gst_collect_pads_set_flush_function    (GstCollectPads *pads,
                                                       GstCollectPadsFlushFunction func,
                                                       gpointer user_data);

/* pad management */

GST_BASE_API
GstCollectData* gst_collect_pads_add_pad       (GstCollectPads *pads, GstPad *pad, guint size,
                                                GstCollectDataDestroyNotify destroy_notify,
                                                gboolean lock);
GST_BASE_API
gboolean        gst_collect_pads_remove_pad    (GstCollectPads *pads, GstPad *pad);

/* start/stop collection */

GST_BASE_API
void            gst_collect_pads_start         (GstCollectPads *pads);

GST_BASE_API
void            gst_collect_pads_stop          (GstCollectPads *pads);

GST_BASE_API
void            gst_collect_pads_set_flushing  (GstCollectPads *pads, gboolean flushing);

/* get collected buffers */

GST_BASE_API
GstBuffer*      gst_collect_pads_peek          (GstCollectPads *pads, GstCollectData *data);

GST_BASE_API
GstBuffer*      gst_collect_pads_pop           (GstCollectPads *pads, GstCollectData *data);

/* get collected bytes */

GST_BASE_API
guint           gst_collect_pads_available     (GstCollectPads *pads);

GST_BASE_API
guint           gst_collect_pads_flush         (GstCollectPads *pads, GstCollectData *data,
                                                guint size);
GST_BASE_API
GstBuffer*      gst_collect_pads_read_buffer   (GstCollectPads * pads, GstCollectData * data,
                                                guint size);
GST_BASE_API
GstBuffer*      gst_collect_pads_take_buffer   (GstCollectPads * pads, GstCollectData * data,
                                                guint size);

/* setting and unsetting waiting mode */

GST_BASE_API
void            gst_collect_pads_set_waiting   (GstCollectPads *pads, GstCollectData *data,
                                                gboolean waiting);

/* convenience helper */

GST_BASE_API
GstFlowReturn	gst_collect_pads_clip_running_time (GstCollectPads * pads,
                                                    GstCollectData * cdata,
                                                    GstBuffer * buf, GstBuffer ** outbuf,
                                                    gpointer user_data);

/* default handlers */

GST_BASE_API
gboolean        gst_collect_pads_event_default (GstCollectPads * pads, GstCollectData * data,
                                                GstEvent * event, gboolean discard);
GST_BASE_API
gboolean        gst_collect_pads_src_event_default (GstCollectPads * pads, GstPad * pad,
                                                    GstEvent * event);
GST_BASE_API
gboolean        gst_collect_pads_query_default (GstCollectPads * pads, GstCollectData * data,
                                                GstQuery * query, gboolean discard);


G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstCollectPads, gst_object_unref)

G_END_DECLS

#endif /* __GST_COLLECT_PADS_H__ */
