/* GstHarness - A test-harness for GStreamer testing
 *
 * Copyright (C) 2012-2015 Pexip <pexip.com>
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __GST_HARNESS_H__
#define __GST_HARNESS_H__

#include <gst/gst.h>
#include <gst/check/gsttestclock.h>
#include <gst/check/check-prelude.h>

G_BEGIN_DECLS

/**
 * GstHarnessThread:
 *
 * Opaque handle representing a GstHarness stress testing thread.
 *
 * Since: 1.6
 */
typedef struct _GstHarnessThread GstHarnessThread;

typedef struct _GstHarness GstHarness;
typedef struct _GstHarnessPrivate GstHarnessPrivate;

/**
 * GstHarness:
 * @element: the element inside the harness
 * @srcpad: the internal harness source pad
 * @sinkpad: the internal harness sink pad
 * @src_harness: the source (input) harness (if any)
 * @sink_harness: the sink (output) harness (if any)
 *
 * Since: 1.6
 */
struct _GstHarness {
  GstElement * element;

  GstPad * srcpad;
  GstPad * sinkpad;

  GstHarness * src_harness;
  GstHarness * sink_harness;

  /*< private >*/
  GstHarnessPrivate * priv;
};

/* Harness creation */

GST_CHECK_API
GstHarness * gst_harness_new_empty (void);

GST_CHECK_API
void         gst_harness_add_element_full (GstHarness           * h,
                                           GstElement           * element,
                                           GstStaticPadTemplate * hsrc,
                                           const gchar          * element_sinkpad_name,
                                           GstStaticPadTemplate * hsink,
                                           const gchar          * element_srcpad_name);

GST_CHECK_API
GstHarness * gst_harness_new_full (GstElement * element,
                                   GstStaticPadTemplate * hsrc,
                                   const gchar          * element_sinkpad_name,
                                   GstStaticPadTemplate * hsink,
                                   const gchar          * element_srcpad_name);

GST_CHECK_API
GstHarness * gst_harness_new_with_element  (GstElement  * element,
                                            const gchar * element_sinkpad_name,
                                            const gchar * element_srcpad_name);

GST_CHECK_API
GstHarness * gst_harness_new_with_padnames (const gchar * element_name,
                                            const gchar * element_sinkpad_name,
                                            const gchar * element_srcpad_name);

GST_CHECK_API
GstHarness * gst_harness_new_with_templates (const gchar * element_name,
                                             GstStaticPadTemplate * hsrc,
                                             GstStaticPadTemplate * hsink);

GST_CHECK_API
GstHarness * gst_harness_new (const gchar * element_name);

GST_CHECK_API
GstHarness * gst_harness_new_parse (const gchar * launchline);

GST_CHECK_API
void         gst_harness_add_parse (GstHarness * h, const gchar * launchline);

GST_CHECK_API
void         gst_harness_teardown (GstHarness * h);

GST_CHECK_API
void         gst_harness_add_element_src_pad  (GstHarness * h, GstPad * srcpad);

GST_CHECK_API
void         gst_harness_add_element_sink_pad (GstHarness * h, GstPad * sinkpad);

/* Caps Functions */

GST_CHECK_API
void         gst_harness_set_src_caps  (GstHarness * h, GstCaps * caps);

GST_CHECK_API
void         gst_harness_set_sink_caps (GstHarness * h, GstCaps * caps);

GST_CHECK_API
void         gst_harness_set_caps (GstHarness * h, GstCaps * in, GstCaps * out);

GST_CHECK_API
void         gst_harness_set_src_caps_str  (GstHarness * h, const gchar * str);

GST_CHECK_API
void         gst_harness_set_sink_caps_str (GstHarness * h, const gchar * str);

GST_CHECK_API
void         gst_harness_set_caps_str (GstHarness  * h,
                                       const gchar * in,
                                       const gchar * out);

/* Clock Functions */

GST_CHECK_API
void           gst_harness_use_systemclock (GstHarness * h);

GST_CHECK_API
void           gst_harness_use_testclock (GstHarness * h);

GST_CHECK_API
GstTestClock * gst_harness_get_testclock (GstHarness * h);

GST_CHECK_API
gboolean       gst_harness_set_time (GstHarness * h, GstClockTime time);

GST_CHECK_API
gboolean       gst_harness_wait_for_clock_id_waits (GstHarness * h,
                                                    guint waits,
                                                    guint timeout);

GST_CHECK_API
gboolean       gst_harness_crank_single_clock_wait (GstHarness * h);

GST_CHECK_API
gboolean       gst_harness_crank_multiple_clock_waits (GstHarness * h,
                                                       guint waits);

/* misc */

GST_CHECK_API
void           gst_harness_play (GstHarness * h);

GST_CHECK_API
void           gst_harness_set_blocking_push_mode (GstHarness * h);

GST_CHECK_API
void           gst_harness_set_forwarding (GstHarness * h, gboolean forwarding);

/* buffers */

GST_CHECK_API
GstBuffer *    gst_harness_create_buffer (GstHarness * h, gsize size);

GST_CHECK_API
GstFlowReturn  gst_harness_push (GstHarness * h, GstBuffer * buffer);

GST_CHECK_API
GstBuffer *    gst_harness_pull (GstHarness * h);

GST_CHECK_API
GstBuffer *    gst_harness_try_pull (GstHarness * h);

GST_CHECK_API
gboolean       gst_harness_pull_until_eos (GstHarness * h, GstBuffer ** buf);

GST_CHECK_API
GstBuffer *    gst_harness_push_and_pull (GstHarness * h, GstBuffer * buffer);

GST_CHECK_API
guint          gst_harness_buffers_received (GstHarness * h);

GST_CHECK_API
guint          gst_harness_buffers_in_queue (GstHarness * h);

GST_CHECK_API
void           gst_harness_set_drop_buffers (GstHarness * h, gboolean drop_buffers);

GST_CHECK_API
void           gst_harness_dump_to_file (GstHarness * h, const gchar * filename);

GST_CHECK_API
guint8 *       gst_harness_take_all_data (GstHarness * h, gsize * size);

GST_CHECK_API
GstBuffer *    gst_harness_take_all_data_as_buffer (GstHarness * h);

GST_CHECK_API
GBytes *       gst_harness_take_all_data_as_bytes (GstHarness * h);

GST_CHECK_API
GstClockTime   gst_harness_get_last_pushed_timestamp (GstHarness * h);

/* downstream events */

GST_CHECK_API
gboolean       gst_harness_push_event (GstHarness * h, GstEvent * event);

GST_CHECK_API
GstEvent *     gst_harness_pull_event (GstHarness * h);

GST_CHECK_API
GstEvent *     gst_harness_try_pull_event  (GstHarness * h);

GST_CHECK_API
guint          gst_harness_events_received (GstHarness * h);

GST_CHECK_API
guint          gst_harness_events_in_queue (GstHarness * h);

/* upstream events */

GST_CHECK_API
gboolean   gst_harness_push_upstream_event (GstHarness * h, GstEvent * event);

GST_CHECK_API
GstEvent * gst_harness_pull_upstream_event (GstHarness * h);

GST_CHECK_API
GstEvent * gst_harness_try_pull_upstream_event  (GstHarness * h);

GST_CHECK_API
guint      gst_harness_upstream_events_received (GstHarness * h);

GST_CHECK_API
guint      gst_harness_upstream_events_in_queue (GstHarness * h);

/* latency */

GST_CHECK_API
GstClockTime gst_harness_query_latency (GstHarness * h);

GST_CHECK_API
void         gst_harness_set_upstream_latency (GstHarness * h, GstClockTime latency);

GST_CHECK_API
void         gst_harness_set_live (GstHarness * h, gboolean is_live);

/* allocation query parameters */

GST_CHECK_API
void         gst_harness_set_propose_allocator (GstHarness                * h,
                                                GstAllocator              * allocator,
                                                const GstAllocationParams * params);

GST_CHECK_API
void         gst_harness_get_allocator         (GstHarness          * h,
                                                GstAllocator       ** allocator,
                                                GstAllocationParams * params);

GST_CHECK_API
void         gst_harness_add_propose_allocation_meta (GstHarness                * h,
                                                      GType                       api,
                                                      const GstStructure        * params);

/* src-harness */

GST_CHECK_API
void          gst_harness_add_src_harness (GstHarness * h,
                                           GstHarness * src_harness,
                                           gboolean has_clock_wait);

GST_CHECK_API
void          gst_harness_add_src (GstHarness  * h,
                                   const gchar * src_element_name,
                                   gboolean      has_clock_wait);

GST_CHECK_API
void          gst_harness_add_src_parse (GstHarness  * h,
                                         const gchar * launchline,
                                         gboolean      has_clock_wait);

GST_CHECK_API
GstFlowReturn gst_harness_push_from_src (GstHarness * h);

GST_CHECK_API
GstFlowReturn gst_harness_src_crank_and_push_many (GstHarness * h,
                                                   gint         cranks,
                                                   gint         pushes);

GST_CHECK_API
gboolean      gst_harness_src_push_event (GstHarness * h);

/* sink-harness */

GST_CHECK_API
void          gst_harness_add_sink_harness (GstHarness * h,
                                            GstHarness * sink_harness);

GST_CHECK_API
void          gst_harness_add_sink (GstHarness  * h,
                                    const gchar * sink_element_name);

GST_CHECK_API
void          gst_harness_add_sink_parse (GstHarness  * h,
                                          const gchar * launchline);

GST_CHECK_API
GstFlowReturn gst_harness_push_to_sink   (GstHarness * h);

GST_CHECK_API
GstFlowReturn gst_harness_sink_push_many (GstHarness * h, gint pushes);

/* convenience functions */

GST_CHECK_API
GstElement *  gst_harness_find_element (GstHarness * h,
                                       const gchar * element_name);

GST_CHECK_API
void          gst_harness_set (GstHarness  * h,
                               const gchar * element_name,
                               const gchar * first_property_name, ...) G_GNUC_NULL_TERMINATED;

GST_CHECK_API
void          gst_harness_get (GstHarness  * h,
                               const gchar * element_name,
                               const gchar * first_property_name, ...) G_GNUC_NULL_TERMINATED;

GST_CHECK_API
void          gst_harness_add_probe (GstHarness        * h,
                                     const gchar       * element_name,
                                     const gchar       * pad_name,
                                     GstPadProbeType     mask,
                                     GstPadProbeCallback callback,
                                     gpointer            user_data,
                                     GDestroyNotify      destroy_data);

/* Stress */

GST_CHECK_API
guint              gst_harness_stress_thread_stop  (GstHarnessThread * t);

GST_CHECK_API
GstHarnessThread * gst_harness_stress_custom_start (GstHarness * h,
                                                    GFunc        init,
                                                    GFunc        callback,
                                                    gpointer     data,
                                                    gulong       sleep);

#define gst_harness_stress_statechange_start(h)                                \
  gst_harness_stress_statechange_start_full (h, G_USEC_PER_SEC / 100)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_statechange_start_full (GstHarness * h,
                                                              gulong       sleep);

#define gst_harness_stress_push_buffer_start(h, c, s, b)                       \
  gst_harness_stress_push_buffer_start_full (h, c, s, b, 0)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_push_buffer_start_full (GstHarness * h,
                                                              GstCaps    * caps,
                                                              const GstSegment * segment,
                                                              GstBuffer  * buf,
                                                              gulong       sleep);

/**
 * GstHarnessPrepareBufferFunc:
 * @h: a #GstHarness
 * @data: user data
 *
 * Since: 1.6
 */
typedef GstBuffer * (*GstHarnessPrepareBufferFunc) (GstHarness * h, gpointer data);

#define gst_harness_stress_push_buffer_with_cb_start(h, c, s, f, d, n)         \
  gst_harness_stress_push_buffer_with_cb_start_full (h, c, s, f, d, n, 0)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_push_buffer_with_cb_start_full (GstHarness   * h,
                                                                      GstCaps      * caps,
                                                                      const GstSegment * segment,
                                                                      GstHarnessPrepareBufferFunc func,
                                                                      gpointer       data,
                                                                      GDestroyNotify notify,
                                                                      gulong         sleep);

#define gst_harness_stress_push_event_start(h, e)                              \
  gst_harness_stress_push_event_start_full (h, e, 0)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_push_event_start_full (GstHarness * h,
                                                             GstEvent   * event,
                                                             gulong       sleep);

/**
 * GstHarnessPrepareEventFunc:
 * @h: a #GstHarness
 * @data: user data
 *
 * Since: 1.8
 */
typedef GstEvent * (*GstHarnessPrepareEventFunc) (GstHarness * h, gpointer data);

#define gst_harness_stress_push_event_with_cb_start(h, f, d, n)                \
  gst_harness_stress_push_event_with_cb_start_full (h, f, d, n, 0)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_push_event_with_cb_start_full (GstHarness   * h,
                                                                     GstHarnessPrepareEventFunc func,
                                                                     gpointer       data,
                                                                     GDestroyNotify notify,
                                                                     gulong         sleep);

#define gst_harness_stress_send_upstream_event_start(h, e)                     \
  gst_harness_stress_push_upstream_event_start_full (h, e, 0)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_push_upstream_event_start_full (GstHarness * h,
                                                                      GstEvent   * event,
                                                                      gulong       sleep);

#define gst_harness_stress_send_upstream_event_with_cb_start(h, f, d, n)       \
  gst_harness_stress_push_upstream_event_with_cb_start_full (h, f, d, n, 0)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_push_upstream_event_with_cb_start_full (GstHarness   * h,
                                                                              GstHarnessPrepareEventFunc func,
                                                                              gpointer       data,
                                                                              GDestroyNotify notify,
                                                                              gulong         sleep);


#define gst_harness_stress_property_start(h, n, v)                             \
  gst_harness_stress_property_start_full (h, n, v, G_USEC_PER_SEC / 1000)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_property_start_full (GstHarness   * h,
                                                           const gchar  * name,
                                                           const GValue * value,
                                                           gulong         sleep);

#define gst_harness_stress_requestpad_start(h, t, n, c, r)                     \
  gst_harness_stress_requestpad_start_full (h, t, n, c, r, G_USEC_PER_SEC / 100)

GST_CHECK_API
GstHarnessThread * gst_harness_stress_requestpad_start_full (GstHarness     * h,
                                                             GstPadTemplate * templ,
                                                             const gchar    * name,
                                                             GstCaps        * caps,
                                                             gboolean         release,
                                                             gulong           sleep);

G_END_DECLS

#endif /* __GST_HARNESS_H__ */
