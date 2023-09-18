/* GStreamer
 * Copyright (C) 2015 Sebastian Dröge <sebastian@centricular.com>
 *
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

#ifndef __GST_PTP_CLOCK_H__
#define __GST_PTP_CLOCK_H__

#include <gst/gst.h>
#include <gst/gstsystemclock.h>
#include <gst/net/net-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_PTP_CLOCK \
  (gst_ptp_clock_get_type())
#define GST_PTP_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_PTP_CLOCK,GstPtpClock))
#define GST_PTP_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_PTP_CLOCK,GstPtpClockClass))
#define GST_IS_PTP_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_PTP_CLOCK))
#define GST_IS_PTP_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_PTP_CLOCK))

typedef struct _GstPtpClock GstPtpClock;
typedef struct _GstPtpClockClass GstPtpClockClass;
typedef struct _GstPtpClockPrivate GstPtpClockPrivate;

/**
 * GstPtpClock:
 *
 * Opaque #GstPtpClock structure.
 */
struct _GstPtpClock {
  GstSystemClock clock;

  /*< private >*/
  GstPtpClockPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstPtpClockClass:
 * @parent_class: parented to #GstSystemClockClass
 *
 * Opaque #GstPtpClockClass structure.
 */
struct _GstPtpClockClass {
  GstSystemClockClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GST_PTP_CLOCK_ID_NONE:
 * PTP clock identification that can be passed to gst_ptp_init() to
 * automatically select one based on the MAC address of interfaces
 */
#define GST_PTP_CLOCK_ID_NONE ((guint64) -1)

GST_NET_API
GType           gst_ptp_clock_get_type             (void);

GST_NET_API
gboolean        gst_ptp_is_supported               (void);

GST_NET_API
gboolean        gst_ptp_is_initialized             (void);

GST_NET_API
gboolean        gst_ptp_init                       (guint64 clock_id,
                                                    gchar ** interfaces);
GST_NET_API
void            gst_ptp_deinit                     (void);

#define GST_PTP_STATISTICS_NEW_DOMAIN_FOUND           "GstPtpStatisticsNewDomainFound"
#define GST_PTP_STATISTICS_BEST_MASTER_CLOCK_SELECTED "GstPtpStatisticsBestMasterClockSelected"
#define GST_PTP_STATISTICS_PATH_DELAY_MEASURED        "GstPtpStatisticsPathDelayMeasured"
#define GST_PTP_STATISTICS_TIME_UPDATED               "GstPtpStatisticsTimeUpdated"

/**
 * GstPtpStatisticsCallback:
 * @domain: PTP domain identifier
 * @stats: New statistics
 * @user_data: Data passed to gst_ptp_statistics_callback_add()
 *
 * The statistics can be the following structures:
 *
 * GST_PTP_STATISTICS_NEW_DOMAIN_FOUND:
 * "domain"                G_TYPE_UINT          The domain identifier of the domain
 * "clock"                 GST_TYPE_CLOCK       The internal clock that is slaved to the
 *                                              PTP domain
 *
 * GST_PTP_STATISTICS_BEST_MASTER_CLOCK_SELECTED:
 * "domain"                G_TYPE_UINT          The domain identifier of the domain
 * "master-clock-id"       G_TYPE_UINT64        PTP clock identifier of the selected master
 *                                              clock
 * "master-clock-port"     G_TYPE_UINT          PTP port number of the selected master clock
 * "grandmaster-clock-id"  G_TYPE_UINT64        PTP clock identifier of the grandmaster clock
 *
 * GST_PTP_STATISTICS_PATH_DELAY_MEASURED:
 * "domain"                G_TYPE_UINT          The domain identifier of the domain
 * "mean-path-delay-avg"   GST_TYPE_CLOCK_TIME  Average mean path delay
 * "mean-path-delay"       GST_TYPE_CLOCK_TIME  Latest mean path delay
 * "delay-request-delay"   GST_TYPE_CLOCK_TIME  Delay of DELAY_REQ / DELAY_RESP messages
 *
 * GST_PTP_STATISTICS_TIME_UPDATED:
 * "domain"                G_TYPE_UINT          The domain identifier of the domain
 * "mean-path-delay-avg"   GST_TYPE_CLOCK_TIME  Average mean path delay
 * "local-time"            GST_TYPE_CLOCK_TIME  Local time that corresponds to ptp-time
 * "ptp-time"              GST_TYPE_CLOCK_TIME  Newly measured PTP time at local-time
 * "estimated-ptp-time"    GST_TYPE_CLOCK_TIME  Estimated PTP time based on previous measurements
 * "discontinuity"         G_TYPE_INT64         Difference between estimated and measured PTP time
 * "synced"                G_TYPE_BOOLEAN       Currently synced to the remote clock
 * "r-squared"             G_TYPE_DOUBLE        R² of clock estimation regression
 * "internal-time"         GST_TYPE_CLOCK_TIME  Internal time clock parameter
 * "external-time"         GST_TYPE_CLOCK_TIME  External time clock parameter
 * "rate-num"              G_TYPE_UINT64        Internal/external rate numerator
 * "rate-den"              G_TYPE_UINT64        Internal/external rate denominator
 * "rate"                  G_TYPE_DOUBLE        Internal/external rate
 *
 * If %FALSE is returned, the callback is removed and never called again.
 *
 */
typedef gboolean  (*GstPtpStatisticsCallback)      (guint8 domain,
                                                    const GstStructure * stats,
                                                    gpointer user_data);
GST_NET_API
gulong          gst_ptp_statistics_callback_add    (GstPtpStatisticsCallback callback,
                                                    gpointer user_data, GDestroyNotify destroy_data);
GST_NET_API
void            gst_ptp_statistics_callback_remove (gulong id);

GST_NET_API
GstClock*       gst_ptp_clock_new                  (const gchar *name,
                                                    guint domain);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstPtpClock, gst_object_unref)

G_END_DECLS

#endif /* __GST_PTP_CLOCK_H__ */

