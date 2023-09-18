/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wim.taymans@chello.be>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *                    2011 Wim Taymans <wim.taymans@gmail.com>
 *
 * gstquery.h: GstQuery API declaration
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


#ifndef __GST_QUERY_H__
#define __GST_QUERY_H__

#include <glib.h>
#include <glib-object.h>
#include <gst/gstconfig.h>

G_BEGIN_DECLS

typedef struct _GstQuery GstQuery;

#include <gst/gstminiobject.h>

/**
 * GstQueryTypeFlags:
 * @GST_QUERY_TYPE_UPSTREAM:     Set if the query can travel upstream.
 * @GST_QUERY_TYPE_DOWNSTREAM:   Set if the query can travel downstream.
 * @GST_QUERY_TYPE_SERIALIZED:   Set if the query should be serialized with data
 *                               flow.
 *
 * #GstQueryTypeFlags indicate the aspects of the different #GstQueryType
 * values. You can get the type flags of a #GstQueryType with the
 * gst_query_type_get_flags() function.
 */
typedef enum {
  GST_QUERY_TYPE_UPSTREAM       = 1 << 0,
  GST_QUERY_TYPE_DOWNSTREAM     = 1 << 1,
  GST_QUERY_TYPE_SERIALIZED     = 1 << 2
} GstQueryTypeFlags;

/**
 * GST_QUERY_TYPE_BOTH: (value 3) (type GstQueryTypeFlags)
 *
 * The same thing as #GST_QUERY_TYPE_UPSTREAM | #GST_QUERY_TYPE_DOWNSTREAM.
 */
#define GST_QUERY_TYPE_BOTH \
    ((GstQueryTypeFlags)(GST_QUERY_TYPE_UPSTREAM | GST_QUERY_TYPE_DOWNSTREAM))

#define GST_QUERY_NUM_SHIFT     (8)

/**
 * GST_QUERY_MAKE_TYPE:
 * @num: the query number to create
 * @flags: the query flags
 *
 * when making custom query types, use this macro with the num and
 * the given flags
 */
#define GST_QUERY_MAKE_TYPE(num,flags) \
    (((num) << GST_QUERY_NUM_SHIFT) | (flags))

#define _FLAG(name) GST_QUERY_TYPE_##name


/**
 * GstQueryType:
 * @GST_QUERY_UNKNOWN: unknown query type
 * @GST_QUERY_POSITION: current position in stream
 * @GST_QUERY_DURATION: total duration of the stream
 * @GST_QUERY_LATENCY: latency of stream
 * @GST_QUERY_JITTER: current jitter of stream
 * @GST_QUERY_RATE: current rate of the stream
 * @GST_QUERY_SEEKING: seeking capabilities
 * @GST_QUERY_SEGMENT: segment start/stop positions
 * @GST_QUERY_CONVERT: convert values between formats
 * @GST_QUERY_FORMATS: query supported formats for convert
 * @GST_QUERY_BUFFERING: query available media for efficient seeking.
 * @GST_QUERY_CUSTOM: a custom application or element defined query.
 * @GST_QUERY_URI: query the URI of the source or sink.
 * @GST_QUERY_ALLOCATION: the buffer allocation properties
 * @GST_QUERY_SCHEDULING: the scheduling properties
 * @GST_QUERY_ACCEPT_CAPS: the accept caps query
 * @GST_QUERY_CAPS: the caps query
 * @GST_QUERY_DRAIN: wait till all serialized data is consumed downstream
 * @GST_QUERY_CONTEXT: query the pipeline-local context from
 *     downstream or upstream (since 1.2)
 * @GST_QUERY_BITRATE: the bitrate query (since 1.16)
 * @GST_QUERY_SELECTABLE: Query stream selection capability (Since: 1.22)
 *
 * Standard predefined Query types
 */
/* NOTE: don't forget to update the table in gstquery.c when changing
 * this enum */
typedef enum {
  GST_QUERY_UNKNOWN      = GST_QUERY_MAKE_TYPE (0, 0),
  GST_QUERY_POSITION     = GST_QUERY_MAKE_TYPE (10, _FLAG(BOTH)),
  GST_QUERY_DURATION     = GST_QUERY_MAKE_TYPE (20, _FLAG(BOTH)),
  GST_QUERY_LATENCY      = GST_QUERY_MAKE_TYPE (30, _FLAG(BOTH)),
  GST_QUERY_JITTER       = GST_QUERY_MAKE_TYPE (40, _FLAG(BOTH)),
  GST_QUERY_RATE         = GST_QUERY_MAKE_TYPE (50, _FLAG(BOTH)),
  GST_QUERY_SEEKING      = GST_QUERY_MAKE_TYPE (60, _FLAG(BOTH)),
  GST_QUERY_SEGMENT      = GST_QUERY_MAKE_TYPE (70, _FLAG(BOTH)),
  GST_QUERY_CONVERT      = GST_QUERY_MAKE_TYPE (80, _FLAG(BOTH)),
  GST_QUERY_FORMATS      = GST_QUERY_MAKE_TYPE (90, _FLAG(BOTH)),
  GST_QUERY_BUFFERING    = GST_QUERY_MAKE_TYPE (110, _FLAG(BOTH)),
  GST_QUERY_CUSTOM       = GST_QUERY_MAKE_TYPE (120, _FLAG(BOTH)),
  GST_QUERY_URI          = GST_QUERY_MAKE_TYPE (130, _FLAG(BOTH)),
  GST_QUERY_ALLOCATION   = GST_QUERY_MAKE_TYPE (140, _FLAG(DOWNSTREAM) | _FLAG(SERIALIZED)),
  GST_QUERY_SCHEDULING   = GST_QUERY_MAKE_TYPE (150, _FLAG(UPSTREAM)),
  GST_QUERY_ACCEPT_CAPS  = GST_QUERY_MAKE_TYPE (160, _FLAG(BOTH)),
  GST_QUERY_CAPS         = GST_QUERY_MAKE_TYPE (170, _FLAG(BOTH)),
  GST_QUERY_DRAIN        = GST_QUERY_MAKE_TYPE (180, _FLAG(DOWNSTREAM) | _FLAG(SERIALIZED)),
  GST_QUERY_CONTEXT      = GST_QUERY_MAKE_TYPE (190, _FLAG(BOTH)),
  GST_QUERY_BITRATE      = GST_QUERY_MAKE_TYPE (200, _FLAG(DOWNSTREAM)),

  /**
   * GST_QUERY_SELECTABLE:
   *
   * Query stream selection capability.
   *
   * Since: 1.22
   */
  GST_QUERY_SELECTABLE   = GST_QUERY_MAKE_TYPE (210, _FLAG(BOTH)),
} GstQueryType;
#undef _FLAG

GST_API GType _gst_query_type;

#define GST_TYPE_QUERY                         (_gst_query_type)
#define GST_IS_QUERY(obj)                      (GST_IS_MINI_OBJECT_TYPE (obj, GST_TYPE_QUERY))
#define GST_QUERY_CAST(obj)                    ((GstQuery*)(obj))
#define GST_QUERY(obj)                         (GST_QUERY_CAST(obj))

/**
 * GST_QUERY_TYPE:
 * @query: the query to query
 *
 * Get the #GstQueryType of the query.
 */
#define GST_QUERY_TYPE(query)  (((GstQuery*)(query))->type)

/**
 * GST_QUERY_TYPE_NAME:
 * @query: the query to query
 *
 * Get a constant string representation of the #GstQueryType of the query.
 */
#define GST_QUERY_TYPE_NAME(query) (gst_query_type_get_name(GST_QUERY_TYPE(query)))

/**
 * GST_QUERY_IS_UPSTREAM:
 * @ev: the query to query
 *
 * Check if an query can travel upstream.
 */
#define GST_QUERY_IS_UPSTREAM(ev)       !!(GST_QUERY_TYPE (ev) & GST_QUERY_TYPE_UPSTREAM)
/**
 * GST_QUERY_IS_DOWNSTREAM:
 * @ev: the query to query
 *
 * Check if an query can travel downstream.
 */
#define GST_QUERY_IS_DOWNSTREAM(ev)     !!(GST_QUERY_TYPE (ev) & GST_QUERY_TYPE_DOWNSTREAM)
/**
 * GST_QUERY_IS_SERIALIZED:
 * @ev: the query to query
 *
 * Check if an query is serialized with the data stream.
 */
#define GST_QUERY_IS_SERIALIZED(ev)     !!(GST_QUERY_TYPE (ev) & GST_QUERY_TYPE_SERIALIZED)


/**
 * GstQuery:
 * @mini_object: The parent #GstMiniObject type
 * @type: the #GstQueryType
 *
 * The #GstQuery structure.
 */
struct _GstQuery
{
  GstMiniObject mini_object;

  /*< public > *//* with COW */
  GstQueryType type;
};

/**
 * GstBufferingMode:
 * @GST_BUFFERING_STREAM: a small amount of data is buffered
 * @GST_BUFFERING_DOWNLOAD: the stream is being downloaded
 * @GST_BUFFERING_TIMESHIFT: the stream is being downloaded in a ringbuffer
 * @GST_BUFFERING_LIVE: the stream is a live stream
 *
 * The different types of buffering methods.
 */
typedef enum {
  GST_BUFFERING_STREAM,
  GST_BUFFERING_DOWNLOAD,
  GST_BUFFERING_TIMESHIFT,
  GST_BUFFERING_LIVE
} GstBufferingMode;

#include <gst/gstiterator.h>
#include <gst/gststructure.h>
#include <gst/gstformat.h>
#include <gst/gstpad.h>
#include <gst/gstallocator.h>
#include <gst/gsttoc.h>
#include <gst/gstcontext.h>

GST_API
const gchar*    gst_query_type_get_name        (GstQueryType type);

GST_API
GQuark          gst_query_type_to_quark        (GstQueryType type);

GST_API
GstQueryTypeFlags
                gst_query_type_get_flags       (GstQueryType type);


GST_API
GType           gst_query_get_type             (void);

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* refcounting */
static inline GstQuery *
gst_query_ref (GstQuery * q)
{
  return GST_QUERY_CAST (gst_mini_object_ref (GST_MINI_OBJECT_CAST (q)));
}

static inline void
gst_query_unref (GstQuery * q)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (q));
}

static inline void
gst_clear_query (GstQuery ** query_ptr)
{
  gst_clear_mini_object ((GstMiniObject **) query_ptr);
}

/* copy query */
static inline GstQuery *
gst_query_copy (const GstQuery * q)
{
  return GST_QUERY_CAST (gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (q)));
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstQuery *  gst_query_ref   (GstQuery * q);

GST_API
void        gst_query_unref (GstQuery * q);

GST_API
void        gst_clear_query (GstQuery ** query_ptr);

GST_API
GstQuery *  gst_query_copy  (const GstQuery * q);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/**
 * gst_query_is_writable:
 * @q: a #GstQuery
 *
 * Tests if you can safely write data into a query's structure.
 */
#define         gst_query_is_writable(q)     gst_mini_object_is_writable (GST_MINI_OBJECT_CAST (q))
/**
 * gst_query_make_writable:
 * @q: (transfer full): a #GstQuery to make writable
 *
 * Makes a writable query from the given query.
 *
 * Returns: (transfer full): a new writable query (possibly same as @q)
 */
#define         gst_query_make_writable(q)      GST_QUERY_CAST (gst_mini_object_make_writable (GST_MINI_OBJECT_CAST (q)))

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
static inline gboolean
gst_query_replace (GstQuery **old_query, GstQuery *new_query)
{
  return gst_mini_object_replace ((GstMiniObject **) old_query, (GstMiniObject *) new_query);
}

static inline gboolean
gst_query_take (GstQuery **old_query, GstQuery *new_query)
{
  return gst_mini_object_take ((GstMiniObject **) old_query,
      (GstMiniObject *) new_query);
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
gboolean        gst_query_replace               (GstQuery ** old_query,
                                                 GstQuery * new_query);

GST_API
gboolean        gst_query_take                  (GstQuery ** old_query,
                                                 GstQuery * new_query);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/* application specific query */

GST_API
GstQuery *      gst_query_new_custom            (GstQueryType type, GstStructure *structure) G_GNUC_MALLOC;

GST_API
const GstStructure *
                gst_query_get_structure         (GstQuery *query);

GST_API
GstStructure *  gst_query_writable_structure    (GstQuery *query);

/* position query */

GST_API
GstQuery*       gst_query_new_position          (GstFormat format) G_GNUC_MALLOC;

GST_API
void            gst_query_set_position          (GstQuery *query, GstFormat format, gint64 cur);

GST_API
void            gst_query_parse_position        (GstQuery *query, GstFormat *format, gint64 *cur);

/* duration query */

GST_API
GstQuery*       gst_query_new_duration          (GstFormat format) G_GNUC_MALLOC;

GST_API
void            gst_query_set_duration          (GstQuery *query, GstFormat format, gint64 duration);

GST_API
void            gst_query_parse_duration        (GstQuery *query, GstFormat *format, gint64 *duration);

/* latency query */

GST_API
GstQuery*       gst_query_new_latency           (void) G_GNUC_MALLOC;

GST_API
void            gst_query_set_latency           (GstQuery *query, gboolean live, GstClockTime min_latency,
                                                 GstClockTime max_latency);

GST_API
void            gst_query_parse_latency         (GstQuery *query, gboolean *live, GstClockTime *min_latency,
                                                 GstClockTime *max_latency);

/* convert query */

GST_API
GstQuery*       gst_query_new_convert           (GstFormat src_format, gint64 value, GstFormat dest_format) G_GNUC_MALLOC;

GST_API
void            gst_query_set_convert           (GstQuery *query, GstFormat src_format, gint64 src_value,
                                                 GstFormat dest_format, gint64 dest_value);

GST_API
void            gst_query_parse_convert         (GstQuery *query, GstFormat *src_format, gint64 *src_value,
                                                 GstFormat *dest_format, gint64 *dest_value);
/* segment query */

GST_API
GstQuery*       gst_query_new_segment           (GstFormat format) G_GNUC_MALLOC;

GST_API
void            gst_query_set_segment           (GstQuery *query, gdouble rate, GstFormat format,
                                                 gint64 start_value, gint64 stop_value);

GST_API
void            gst_query_parse_segment         (GstQuery *query, gdouble *rate, GstFormat *format,
                                                 gint64 *start_value, gint64 *stop_value);

/* seeking query */

GST_API
GstQuery*       gst_query_new_seeking           (GstFormat format) G_GNUC_MALLOC;

GST_API
void            gst_query_set_seeking           (GstQuery *query, GstFormat format,
                                                 gboolean seekable,
                                                 gint64 segment_start,
                                                 gint64 segment_end);

GST_API
void            gst_query_parse_seeking         (GstQuery *query, GstFormat *format,
                                                 gboolean *seekable,
                                                 gint64 *segment_start,
                                                 gint64 *segment_end);
/* formats query */

GST_API
GstQuery*       gst_query_new_formats           (void) G_GNUC_MALLOC;

GST_API
void            gst_query_set_formats           (GstQuery *query, gint n_formats, ...);

GST_API
void            gst_query_set_formatsv          (GstQuery *query, gint n_formats, const GstFormat *formats);

GST_API
void            gst_query_parse_n_formats       (GstQuery *query, guint *n_formats);

GST_API
void            gst_query_parse_nth_format      (GstQuery *query, guint nth, GstFormat *format);

/* buffering query */

GST_API
GstQuery*       gst_query_new_buffering           (GstFormat format) G_GNUC_MALLOC;

GST_API
void            gst_query_set_buffering_percent   (GstQuery *query, gboolean busy, gint percent);

GST_API
void            gst_query_parse_buffering_percent (GstQuery *query, gboolean *busy, gint *percent);

GST_API
void            gst_query_set_buffering_stats     (GstQuery *query, GstBufferingMode mode,
                                                   gint avg_in, gint avg_out,
                                                   gint64 buffering_left);

GST_API
void            gst_query_parse_buffering_stats    (GstQuery *query, GstBufferingMode *mode,
                                                   gint *avg_in, gint *avg_out,
                                                   gint64 *buffering_left);

GST_API
void            gst_query_set_buffering_range     (GstQuery *query, GstFormat format,
                                                   gint64 start, gint64 stop,
                                                   gint64 estimated_total);

GST_API
void            gst_query_parse_buffering_range   (GstQuery *query, GstFormat *format,
                                                   gint64 *start, gint64 *stop,
                                                   gint64 *estimated_total);

GST_API
gboolean        gst_query_add_buffering_range       (GstQuery *query,
                                                     gint64 start, gint64 stop);

GST_API
guint           gst_query_get_n_buffering_ranges    (GstQuery *query);

GST_API
gboolean        gst_query_parse_nth_buffering_range (GstQuery *query,
                                                     guint index, gint64 *start,
                                                     gint64 *stop);

/* URI query */

GST_API
GstQuery *      gst_query_new_uri                    (void) G_GNUC_MALLOC;

GST_API
void            gst_query_parse_uri                  (GstQuery *query, gchar **uri);

GST_API
void            gst_query_set_uri                    (GstQuery *query, const gchar *uri);

GST_API
void            gst_query_parse_uri_redirection      (GstQuery *query, gchar **uri);

GST_API
void            gst_query_set_uri_redirection        (GstQuery *query, const gchar *uri);

GST_API
void            gst_query_parse_uri_redirection_permanent (GstQuery *query, gboolean * permanent);

GST_API
void            gst_query_set_uri_redirection_permanent (GstQuery *query, gboolean permanent);

/* allocation query */

GST_API
GstQuery *      gst_query_new_allocation             (GstCaps *caps, gboolean need_pool) G_GNUC_MALLOC;

GST_API
void            gst_query_parse_allocation           (GstQuery *query, GstCaps **caps, gboolean *need_pool);

/* pools */

GST_API
void            gst_query_add_allocation_pool        (GstQuery *query, GstBufferPool *pool,
                                                      guint size, guint min_buffers,
                                                      guint max_buffers);

GST_API
guint           gst_query_get_n_allocation_pools     (GstQuery *query);

GST_API
void            gst_query_parse_nth_allocation_pool  (GstQuery *query, guint index,
                                                      GstBufferPool **pool,
                                                      guint *size, guint *min_buffers,
                                                      guint *max_buffers);

GST_API
void            gst_query_set_nth_allocation_pool    (GstQuery *query, guint index,
                                                      GstBufferPool *pool,
                                                      guint size, guint min_buffers,
                                                      guint max_buffers);

GST_API
void            gst_query_remove_nth_allocation_pool (GstQuery *query, guint index);

/* allocators */

GST_API
void            gst_query_add_allocation_param       (GstQuery *query, GstAllocator *allocator,
                                                      const GstAllocationParams *params);

GST_API
guint           gst_query_get_n_allocation_params    (GstQuery *query);

GST_API
void            gst_query_parse_nth_allocation_param (GstQuery *query, guint index,
                                                      GstAllocator **allocator,
                                                      GstAllocationParams *params);

GST_API
void            gst_query_set_nth_allocation_param   (GstQuery *query, guint index,
                                                      GstAllocator *allocator,
                                                      const GstAllocationParams *params);

GST_API
void            gst_query_remove_nth_allocation_param (GstQuery *query, guint index);

/* metadata */

GST_API
void            gst_query_add_allocation_meta        (GstQuery *query, GType api, const GstStructure *params);

GST_API
guint           gst_query_get_n_allocation_metas     (GstQuery *query);

GST_API
GType           gst_query_parse_nth_allocation_meta  (GstQuery *query, guint index,
                                                      const GstStructure **params);

GST_API
void            gst_query_remove_nth_allocation_meta (GstQuery *query, guint index);

GST_API
gboolean        gst_query_find_allocation_meta       (GstQuery *query, GType api, guint *index);


/* scheduling query */
/**
 * GstSchedulingFlags:
 * @GST_SCHEDULING_FLAG_SEEKABLE: if seeking is possible
 * @GST_SCHEDULING_FLAG_SEQUENTIAL: if sequential access is recommended
 * @GST_SCHEDULING_FLAG_BANDWIDTH_LIMITED: if bandwidth is limited and buffering possible (since 1.2)
 *
 * The different scheduling flags.
 */
typedef enum {
  GST_SCHEDULING_FLAG_SEEKABLE          = (1 << 0),
  GST_SCHEDULING_FLAG_SEQUENTIAL        = (1 << 1),
  GST_SCHEDULING_FLAG_BANDWIDTH_LIMITED = (1 << 2)
} GstSchedulingFlags;

GST_API
GstQuery *      gst_query_new_scheduling          (void) G_GNUC_MALLOC;

GST_API
void            gst_query_set_scheduling          (GstQuery *query, GstSchedulingFlags flags,
                                                   gint minsize, gint maxsize, gint align);

GST_API
void            gst_query_parse_scheduling        (GstQuery *query, GstSchedulingFlags *flags,
                                                   gint *minsize, gint *maxsize, gint *align);

GST_API
void            gst_query_add_scheduling_mode       (GstQuery *query, GstPadMode mode);

GST_API
guint           gst_query_get_n_scheduling_modes    (GstQuery *query);

GST_API
GstPadMode      gst_query_parse_nth_scheduling_mode (GstQuery *query, guint index);

GST_API
gboolean        gst_query_has_scheduling_mode       (GstQuery *query, GstPadMode mode);

GST_API
gboolean        gst_query_has_scheduling_mode_with_flags (GstQuery * query, GstPadMode mode,
                                                    GstSchedulingFlags flags);

/* accept-caps query */

GST_API
GstQuery *      gst_query_new_accept_caps          (GstCaps *caps) G_GNUC_MALLOC;

GST_API
void            gst_query_parse_accept_caps        (GstQuery *query, GstCaps **caps);

GST_API
void            gst_query_set_accept_caps_result   (GstQuery *query, gboolean result);

GST_API
void            gst_query_parse_accept_caps_result (GstQuery *query, gboolean *result);

/* caps query */

GST_API
GstQuery *      gst_query_new_caps                 (GstCaps *filter) G_GNUC_MALLOC;

GST_API
void            gst_query_parse_caps               (GstQuery *query, GstCaps **filter);

GST_API
void            gst_query_set_caps_result          (GstQuery *query, GstCaps *caps);

GST_API
void            gst_query_parse_caps_result        (GstQuery *query, GstCaps **caps);

/* drain query */

GST_API
GstQuery *      gst_query_new_drain                (void) G_GNUC_MALLOC;

/* context query */

GST_API
GstQuery *      gst_query_new_context              (const gchar * context_type) G_GNUC_MALLOC;

GST_API
gboolean        gst_query_parse_context_type       (GstQuery * query, const gchar ** context_type);

GST_API
void            gst_query_set_context              (GstQuery *query, GstContext *context);

GST_API
void            gst_query_parse_context            (GstQuery *query, GstContext **context);

/* bitrate query */

GST_API
GstQuery *      gst_query_new_bitrate              (void) G_GNUC_MALLOC;

GST_API
void            gst_query_set_bitrate              (GstQuery * query, guint nominal_bitrate);

GST_API
void            gst_query_parse_bitrate            (GstQuery * query, guint * nominal_bitrate);

/* selectable query */

GST_API
GstQuery *      gst_query_new_selectable           (void) G_GNUC_MALLOC;

GST_API
void            gst_query_set_selectable           (GstQuery *query, gboolean selectable);

GST_API
void            gst_query_parse_selectable         (GstQuery *query, gboolean * selectable);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstQuery, gst_query_unref)

G_END_DECLS

#endif /* __GST_QUERY_H__ */

