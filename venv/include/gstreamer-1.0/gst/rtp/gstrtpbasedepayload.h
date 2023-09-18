/* GStreamer
 * Copyright (C) <2005> Philippe Khalaf <burger@speedy.org>
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

#ifndef __GST_RTP_BASE_DEPAYLOAD_H__
#define __GST_RTP_BASE_DEPAYLOAD_H__

#include <gst/gst.h>
#include <gst/rtp/gstrtpbuffer.h>

G_BEGIN_DECLS

#define GST_TYPE_RTP_BASE_DEPAYLOAD (gst_rtp_base_depayload_get_type())
#define GST_RTP_BASE_DEPAYLOAD(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_RTP_BASE_DEPAYLOAD,GstRTPBaseDepayload))
#define GST_RTP_BASE_DEPAYLOAD_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_RTP_BASE_DEPAYLOAD,GstRTPBaseDepayloadClass))
#define GST_RTP_BASE_DEPAYLOAD_GET_CLASS(obj) \
        (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_RTP_BASE_DEPAYLOAD,GstRTPBaseDepayloadClass))
#define GST_IS_RTP_BASE_DEPAYLOAD(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_RTP_BASE_DEPAYLOAD))
#define GST_IS_RTP_BASE_DEPAYLOAD_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_RTP_BASE_DEPAYLOAD))
#define GST_RTP_BASE_DEPAYLOAD_CAST(obj) ((GstRTPBaseDepayload *)(obj))

#define GST_RTP_BASE_DEPAYLOAD_SINKPAD(depayload) (GST_RTP_BASE_DEPAYLOAD_CAST (depayload)->sinkpad)
#define GST_RTP_BASE_DEPAYLOAD_SRCPAD(depayload)  (GST_RTP_BASE_DEPAYLOAD_CAST (depayload)->srcpad)

typedef struct _GstRTPBaseDepayload      GstRTPBaseDepayload;
typedef struct _GstRTPBaseDepayloadClass GstRTPBaseDepayloadClass;
typedef struct _GstRTPBaseDepayloadPrivate GstRTPBaseDepayloadPrivate;

struct _GstRTPBaseDepayload
{
  GstElement parent;

  GstPad *sinkpad, *srcpad;

  /* this attribute must be set by the child */
  guint clock_rate;

  GstSegment segment;
  gboolean need_newsegment;

  /*< private >*/
  GstRTPBaseDepayloadPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstRTPBaseDepayloadClass:
 * @parent_class: the parent class
 * @set_caps: configure the depayloader
 * @process: process incoming rtp packets. Subclass must implement either
 *   this method or @process_rtp_packet to process incoming rtp packets.
 *   If the child returns a buffer without a valid timestamp, the timestamp
 *   of the provided buffer will be applied to the result buffer and the
 *   buffer will be pushed. If this function returns %NULL, nothing is pushed.
 * @packet_lost: signal the depayloader about packet loss
 * @handle_event: custom event handling
 * @process_rtp_packet: Same as the process virtual function, but slightly more
 * efficient, since it is passed the rtp buffer structure that has already
 * been mapped (with GST_MAP_READ) by the base class and thus does not have
 * to be mapped again by the subclass. Can be used by the subclass to process
 * incoming rtp packets. If the subclass returns a buffer without a valid
 * timestamp, the timestamp of the input buffer will be applied to the result
 * buffer and the output buffer will be pushed out. If this function returns
 * %NULL, nothing is pushed out. Since: 1.6.
 *
 * Base class for RTP depayloaders.
 */
struct _GstRTPBaseDepayloadClass
{
  GstElementClass parent_class;

  /*< public >*/
  /* virtuals, inform the subclass of the caps. */
  gboolean (*set_caps) (GstRTPBaseDepayload *filter, GstCaps *caps);

  /* pure virtual function */
  GstBuffer * (*process) (GstRTPBaseDepayload *base, GstBuffer *in);

  /* non-pure function used to to signal the depayloader about packet loss. the
   * timestamp and duration are the estimated values of the lost packet.
   * The default implementation of this message pushes a segment update. */
  gboolean (*packet_lost) (GstRTPBaseDepayload *filter, GstEvent *event);

  /* the default implementation does the default actions for events but
   * implementation can override. */
  gboolean (*handle_event) (GstRTPBaseDepayload * filter, GstEvent * event);

  GstBuffer * (*process_rtp_packet) (GstRTPBaseDepayload *base, GstRTPBuffer * rtp_buffer);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING - 1];
};

GST_RTP_API
GType gst_rtp_base_depayload_get_type (void);

GST_RTP_API
GstFlowReturn   gst_rtp_base_depayload_push       (GstRTPBaseDepayload *filter, GstBuffer *out_buf);

GST_RTP_API
GstFlowReturn   gst_rtp_base_depayload_push_list  (GstRTPBaseDepayload *filter, GstBufferList *out_list);

GST_RTP_API
gboolean        gst_rtp_base_depayload_is_source_info_enabled  (GstRTPBaseDepayload * depayload);

GST_RTP_API
void            gst_rtp_base_depayload_set_source_info_enabled (GstRTPBaseDepayload * depayload,
                                                                gboolean enable);


G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstRTPBaseDepayload, gst_object_unref)

G_END_DECLS

#endif /* __GST_RTP_BASE_DEPAYLOAD_H__ */
