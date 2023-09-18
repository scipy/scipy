/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstclock.h: Header for clock subsystem
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

#ifndef __GST_CLOCK_H__
#define __GST_CLOCK_H__

#include <gst/gstconfig.h>
#include <glib.h>

G_BEGIN_DECLS

/* --- standard type macros --- */
#define GST_TYPE_CLOCK                  (gst_clock_get_type ())
#define GST_CLOCK(clock)                (G_TYPE_CHECK_INSTANCE_CAST ((clock), GST_TYPE_CLOCK, GstClock))
#define GST_IS_CLOCK(clock)             (G_TYPE_CHECK_INSTANCE_TYPE ((clock), GST_TYPE_CLOCK))
#define GST_CLOCK_CLASS(cclass)         (G_TYPE_CHECK_CLASS_CAST ((cclass), GST_TYPE_CLOCK, GstClockClass))
#define GST_IS_CLOCK_CLASS(cclass)      (G_TYPE_CHECK_CLASS_TYPE ((cclass), GST_TYPE_CLOCK))
#define GST_CLOCK_GET_CLASS(clock)      (G_TYPE_INSTANCE_GET_CLASS ((clock), GST_TYPE_CLOCK, GstClockClass))
#define GST_CLOCK_CAST(clock)           ((GstClock*)(clock))

/**
 * GstClockTime:
 *
 * A datatype to hold a time, measured in nanoseconds.
 */
typedef guint64 GstClockTime;

/**
 * GST_TYPE_CLOCK_TIME:
 *
 * The #GType of a #GstClockTime.
 */
#define GST_TYPE_CLOCK_TIME G_TYPE_UINT64

/**
 * GstClockTimeDiff:
 *
 * A datatype to hold a time difference, measured in nanoseconds.
 */
typedef gint64 GstClockTimeDiff;
/**
 * GstClockID:
 *
 * A datatype to hold the handle to an outstanding sync or async clock callback.
 */
typedef gpointer GstClockID;

/**
 * GST_CLOCK_TIME_NONE: (value 18446744073709551615) (type GstClockTime)
 *
 * Constant to define an undefined clock time.
 */
#define GST_CLOCK_TIME_NONE             ((GstClockTime) -1)
/**
 * GST_CLOCK_TIME_IS_VALID:
 * @time: clock time to validate
 *
 * Tests if a given #GstClockTime represents a valid defined time.
 */
#define GST_CLOCK_TIME_IS_VALID(time)   (((GstClockTime)(time)) != GST_CLOCK_TIME_NONE)

/**
 * GST_CLOCK_STIME_NONE: (value -9223372036854775808) (type GstClockTimeDiff)
 *
 * Constant to define an undefined clock time.
 */
#define GST_CLOCK_STIME_NONE             ((GstClockTimeDiff)G_MININT64)
/**
 * GST_CLOCK_STIME_IS_VALID:
 * @time: signed clock time to validate
 *
 * Tests if a given #GstClockTimeDiff of #gint64 represents a valid defined time.
 *
 * Since: 1.6
 */
#define GST_CLOCK_STIME_IS_VALID(time)   (((GstClockTimeDiff)(time)) != GST_CLOCK_STIME_NONE)

/**
 * GST_SECOND: (value 1000000000) (type GstClockTimeDiff)
 *
 * Constant that defines one GStreamer second.
 */
#define GST_SECOND  ((GstClockTimeDiff)(G_USEC_PER_SEC * G_GINT64_CONSTANT (1000)))
/**
 * GST_MSECOND: (value 1000000) (type GstClockTimeDiff)
 *
 * Constant that defines one GStreamer millisecond.
 */
#define GST_MSECOND ((GstClockTimeDiff)(GST_SECOND / G_GINT64_CONSTANT (1000)))
/**
 * GST_USECOND: (value 1000) (type GstClockTimeDiff)
 *
 * Constant that defines one GStreamer microsecond.
 */
#define GST_USECOND ((GstClockTimeDiff)(GST_SECOND / G_GINT64_CONSTANT (1000000)))
/**
 * GST_NSECOND: (value 1) (type GstClockTimeDiff)
 *
 * Constant that defines one GStreamer nanosecond
 */
#define GST_NSECOND ((GstClockTimeDiff)(GST_SECOND / G_GINT64_CONSTANT (1000000000)))


/**
 * GST_TIME_AS_SECONDS:
 * @time: the time
 *
 * Converts a #GstClockTime to seconds.
 */
#define GST_TIME_AS_SECONDS(time)  ((time) / GST_SECOND)
/**
 * GST_TIME_AS_MSECONDS:
 * @time: the time
 *
 * Converts a #GstClockTime to milliseconds (1/1000 of a second).
 */
#define GST_TIME_AS_MSECONDS(time) ((time) / G_GINT64_CONSTANT (1000000))
/**
 * GST_TIME_AS_USECONDS:
 * @time: the time
 *
 * Converts a #GstClockTime to microseconds (1/1000000 of a second).
 */
#define GST_TIME_AS_USECONDS(time) ((time) / G_GINT64_CONSTANT (1000))
/**
 * GST_TIME_AS_NSECONDS:
 * @time: the time
 *
 * Converts a #GstClockTime to nanoseconds (1/1000000000 of a second).
 */
#define GST_TIME_AS_NSECONDS(time) (time)

/**
 * GST_CLOCK_DIFF:
 * @s: the first time
 * @e: the second time
 *
 * Calculates a difference between two clock times as a #GstClockTimeDiff.
 * The difference is calculated as @e - @s.
 */
#define GST_CLOCK_DIFF(s, e)            (GstClockTimeDiff)((e) - (s))

/**
 * GST_TIMEVAL_TO_TIME:
 * @tv: the timeval to convert
 *
 * Converts a GTimeVal to a #GstClockTime.
 */
#define GST_TIMEVAL_TO_TIME(tv)         (GstClockTime)((tv).tv_sec * GST_SECOND + (tv).tv_usec * GST_USECOND)

/**
 * GST_TIME_TO_TIMEVAL:
 * @t: The #GstClockTime to convert
 * @tv: The target timeval
 *
 * Converts a #GstClockTime to a GTimeVal
 *
 * > on 32-bit systems, a timeval has a range of only 2^32 - 1 seconds,
 * > which is about 68 years.  Expect trouble if you want to schedule stuff
 * > in your pipeline for 2038.
 */
#define GST_TIME_TO_TIMEVAL(t,tv)                               \
G_STMT_START {                                                  \
  g_assert ("Value of time " #t " is out of timeval's range" && \
      ((t) / GST_SECOND) < G_MAXLONG);                          \
  (tv).tv_sec  = (glong) (((GstClockTime) (t)) / GST_SECOND);   \
  (tv).tv_usec = (glong) ((((GstClockTime) (t)) -               \
                  ((GstClockTime) (tv).tv_sec) * GST_SECOND)    \
                 / GST_USECOND);                                \
} G_STMT_END

/**
 * GST_TIMESPEC_TO_TIME:
 * @ts: the timespec to convert
 *
 * Converts a struct timespec (see `man pselect`) to a #GstClockTime.
 */
#define GST_TIMESPEC_TO_TIME(ts)        (GstClockTime)((ts).tv_sec * GST_SECOND + (ts).tv_nsec * GST_NSECOND)
/**
 * GST_TIME_TO_TIMESPEC:
 * @t: The #GstClockTime to convert
 * @ts: The target timespec
 *
 * Converts a #GstClockTime to a struct timespec (see `man pselect`)
 */
#define GST_TIME_TO_TIMESPEC(t,ts)                                \
G_STMT_START {                                                    \
  g_assert ("Value of time " #t " is out of timespec's range" &&  \
      ((t) / GST_SECOND) < G_MAXLONG);                            \
  (ts).tv_sec  =  (glong) ((t) / GST_SECOND);                     \
  (ts).tv_nsec = (glong) (((t) - (ts).tv_sec * GST_SECOND) / GST_NSECOND);        \
} G_STMT_END

/* timestamp debugging macros */
/**
 * GST_TIME_FORMAT: (skip):
 *
 * A string that can be used in printf-like format strings to display a
 * #GstClockTime value in `h:m:s` format.  Use GST_TIME_ARGS() to construct
 * the matching arguments.
 *
 * Example:
 *
 * ``` C
 * printf("%" GST_TIME_FORMAT "\n", GST_TIME_ARGS(ts));
 * ``` 
 */
#define GST_TIME_FORMAT "u:%02u:%02u.%09u"
/**
 * GST_TIME_ARGS: (skip):
 * @t: a #GstClockTime
 *
 * Formats @t for the #GST_TIME_FORMAT format string. Note: @t will be
 * evaluated more than once.
 */
#define GST_TIME_ARGS(t) \
        GST_CLOCK_TIME_IS_VALID (t) ? \
        (guint) (((GstClockTime)(t)) / (GST_SECOND * 60 * 60)) : 99, \
        GST_CLOCK_TIME_IS_VALID (t) ? \
        (guint) ((((GstClockTime)(t)) / (GST_SECOND * 60)) % 60) : 99, \
        GST_CLOCK_TIME_IS_VALID (t) ? \
        (guint) ((((GstClockTime)(t)) / GST_SECOND) % 60) : 99, \
        GST_CLOCK_TIME_IS_VALID (t) ? \
        (guint) (((GstClockTime)(t)) % GST_SECOND) : 999999999
/**
 * GST_STIME_FORMAT: (skip):
 *
 * A string that can be used in printf-like format strings to display a signed
 * #GstClockTimeDiff or #gint64 value in `h:m:s` format.  Use GST_TIME_ARGS() to
 * construct the matching arguments.
 *
 * Example:
 *
 * ``` C
 * printf("%" GST_STIME_FORMAT "\n", GST_STIME_ARGS(ts));
 * ```
 *
 * Since: 1.6
 */
#define GST_STIME_FORMAT "c%" GST_TIME_FORMAT
/**
 * GST_STIME_ARGS: (skip):
 * @t: a #GstClockTimeDiff or #gint64
 *
 * Formats @t for the #GST_STIME_FORMAT format string. Note: @t will be
 * evaluated more than once.
 *
 * Since: 1.6
 */
#define GST_STIME_ARGS(t)						\
  ((t) == GST_CLOCK_STIME_NONE || (t) >= 0) ? '+' : '-',		\
    GST_CLOCK_STIME_IS_VALID (t) ?					\
    (guint) (((GstClockTime)(ABS(t))) / (GST_SECOND * 60 * 60)) : 99,	\
    GST_CLOCK_STIME_IS_VALID (t) ?					\
    (guint) ((((GstClockTime)(ABS(t))) / (GST_SECOND * 60)) % 60) : 99,	\
    GST_CLOCK_STIME_IS_VALID (t) ?					\
    (guint) ((((GstClockTime)(ABS(t))) / GST_SECOND) % 60) : 99,	\
    GST_CLOCK_STIME_IS_VALID (t) ?					\
    (guint) (((GstClockTime)(ABS(t))) % GST_SECOND) : 999999999

typedef struct _GstClockEntry   GstClockEntry;
typedef struct _GstClock        GstClock;
typedef struct _GstClockClass   GstClockClass;
typedef struct _GstClockPrivate GstClockPrivate;

/* --- prototype for async callbacks --- */
/**
 * GstClockCallback:
 * @clock: The clock that triggered the callback
 * @time: The time it was triggered
 * @id: The #GstClockID that expired
 * @user_data: user data passed in the gst_clock_id_wait_async() function
 *
 * The function prototype of the callback.
 *
 * Returns: %TRUE or %FALSE (currently unused)
 */
typedef gboolean        (*GstClockCallback)     (GstClock *clock, GstClockTime time,
                                                 GstClockID id, gpointer user_data);
/**
 * GstClockReturn:
 * @GST_CLOCK_OK: The operation succeeded.
 * @GST_CLOCK_EARLY: The operation was scheduled too late.
 * @GST_CLOCK_UNSCHEDULED: The clockID was unscheduled
 * @GST_CLOCK_BUSY: The ClockID is busy
 * @GST_CLOCK_BADTIME: A bad time was provided to a function.
 * @GST_CLOCK_ERROR: An error occurred
 * @GST_CLOCK_UNSUPPORTED: Operation is not supported
 * @GST_CLOCK_DONE: The ClockID is done waiting
 *
 * The return value of a clock operation.
 */
typedef enum
{
  GST_CLOCK_OK          =  0,
  GST_CLOCK_EARLY       =  1,
  GST_CLOCK_UNSCHEDULED =  2,
  GST_CLOCK_BUSY        =  3,
  GST_CLOCK_BADTIME     =  4,
  GST_CLOCK_ERROR       =  5,
  GST_CLOCK_UNSUPPORTED =  6,
  GST_CLOCK_DONE        =  7
} GstClockReturn;

/**
 * GstClockEntryType:
 * @GST_CLOCK_ENTRY_SINGLE: a single shot timeout
 * @GST_CLOCK_ENTRY_PERIODIC: a periodic timeout request
 *
 * The type of the clock entry
 */
typedef enum {
  GST_CLOCK_ENTRY_SINGLE,
  GST_CLOCK_ENTRY_PERIODIC
} GstClockEntryType;

/**
 * GST_CLOCK_ENTRY:
 * @entry: the entry to cast
 *
 * Casts to a clock entry
 */
#define GST_CLOCK_ENTRY(entry)          ((GstClockEntry *)(entry))

#ifndef GST_DISABLE_DEPRECATED
/**
 * GST_CLOCK_ENTRY_CLOCK:
 * @entry: the entry to query
 *
 * Gets the owner clock of the entry
 *
 * Deprecated: Use gst_clock_id_get_clock() instead.
 */
#define GST_CLOCK_ENTRY_CLOCK(entry)    ((entry)->clock)
#endif
/**
 * GST_CLOCK_ENTRY_TYPE:
 * @entry: the entry to query
 *
 * Gets the type of the clock entry
 */
#define GST_CLOCK_ENTRY_TYPE(entry)     ((entry)->type)
/**
 * GST_CLOCK_ENTRY_TIME:
 * @entry: the entry to query
 *
 * Gets the requested time of this entry
 */
#define GST_CLOCK_ENTRY_TIME(entry)     ((entry)->time)
/**
 * GST_CLOCK_ENTRY_INTERVAL:
 * @entry: the entry to query
 *
 * Gets the interval of this periodic entry
 */
#define GST_CLOCK_ENTRY_INTERVAL(entry) ((entry)->interval)
/**
 * GST_CLOCK_ENTRY_STATUS:
 * @entry: the entry to query
 *
 * The status of the entry
 */
#define GST_CLOCK_ENTRY_STATUS(entry)   ((entry)->status)

/**
 * GstClockEntry:
 * @refcount: reference counter (read-only)
 *
 * All pending timeouts or periodic notifies are converted into
 * an entry.
 * Note that GstClockEntry should be treated as an opaque structure. It must
 * not be extended or allocated using a custom allocator.
 */
struct _GstClockEntry {
  gint                  refcount;
  /*< private >*/
#ifndef GST_REMOVE_DEPRECATED
#ifndef GST_DISABLE_DEPRECATED
  GstClock              *clock;
#else
  gpointer               _clock;
#endif
#endif
  GstClockEntryType      type;
  GstClockTime           time;
  GstClockTime           interval;
  GstClockReturn         status;
  GstClockCallback       func;
  gpointer               user_data;
  GDestroyNotify         destroy_data;
  gboolean               unscheduled;
  gboolean               woken_up;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

#include <gst/gstobject.h>

/**
 * GstClockFlags:
 * @GST_CLOCK_FLAG_CAN_DO_SINGLE_SYNC: clock can do a single sync timeout request
 * @GST_CLOCK_FLAG_CAN_DO_SINGLE_ASYNC: clock can do a single async timeout request
 * @GST_CLOCK_FLAG_CAN_DO_PERIODIC_SYNC: clock can do sync periodic timeout requests
 * @GST_CLOCK_FLAG_CAN_DO_PERIODIC_ASYNC: clock can do async periodic timeout callbacks
 * @GST_CLOCK_FLAG_CAN_SET_RESOLUTION: clock's resolution can be changed
 * @GST_CLOCK_FLAG_CAN_SET_MASTER: clock can be slaved to a master clock
 * @GST_CLOCK_FLAG_LAST: subclasses can add additional flags starting from this flag
 *
 * The capabilities of this clock
 */
typedef enum {
  GST_CLOCK_FLAG_CAN_DO_SINGLE_SYNC     = (GST_OBJECT_FLAG_LAST << 0),
  GST_CLOCK_FLAG_CAN_DO_SINGLE_ASYNC    = (GST_OBJECT_FLAG_LAST << 1),
  GST_CLOCK_FLAG_CAN_DO_PERIODIC_SYNC   = (GST_OBJECT_FLAG_LAST << 2),
  GST_CLOCK_FLAG_CAN_DO_PERIODIC_ASYNC  = (GST_OBJECT_FLAG_LAST << 3),
  GST_CLOCK_FLAG_CAN_SET_RESOLUTION     = (GST_OBJECT_FLAG_LAST << 4),
  GST_CLOCK_FLAG_CAN_SET_MASTER         = (GST_OBJECT_FLAG_LAST << 5),

  /**
   * GST_CLOCK_FLAG_NEEDS_STARTUP_SYNC:
   *
   * clock needs to be synced before it can be used
   *
   * Since: 1.6
   */
  GST_CLOCK_FLAG_NEEDS_STARTUP_SYNC     = (GST_OBJECT_FLAG_LAST << 6),
  /* padding */
  GST_CLOCK_FLAG_LAST                   = (GST_OBJECT_FLAG_LAST << 8)
} GstClockFlags;

/**
 * GST_CLOCK_FLAGS:
 * @clock: the clock to query
 *
 * Gets the #GstClockFlags clock flags.
 */
#define GST_CLOCK_FLAGS(clock)  GST_OBJECT_FLAGS(clock)

/**
 * GstClock:
 * @object: the parent structure
 *
 * #GstClock base structure.
 */
struct _GstClock {
  GstObject      object;

  /*< private >*/
  GstClockPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstClockClass:
 * @parent_class: the parent class structure
 *
 * GStreamer clock class. Override the vmethods to implement the clock
 * functionality.
 */
struct _GstClockClass {
  GstObjectClass        parent_class;

  /*< public >*/
  /* vtable */

  /**
   * GstClockClass::change_resolution:
   * @clock: the #GstClock
   * @old_resolution: the previous resolution
   * @new_resolution: the new resolution
   *
   * Change the resolution of the clock. Not all values might
   * be acceptable.
   *
   * Returns: the new resolution
   */
  GstClockTime          (*change_resolution)    (GstClock *clock,
                                                 GstClockTime old_resolution,
                                                 GstClockTime new_resolution);

  /**
   * GstClockClass::get_resolution:
   * @clock: the #GstClock
   *
   * Get the resolution of the clock.
   *
   * Returns: the current resolution
   */
  GstClockTime          (*get_resolution)       (GstClock *clock);

  /**
   * GstClockClass::get_internal_time:
   * @clock: the #GstClock
   *
   * Get the internal unadjusted time of the clock.
   *
   * Implement #GstClockClass::wait instead.
   *
   * Returns: the internal time
   */
  GstClockTime          (*get_internal_time)    (GstClock *clock);

  /* waiting on an ID */

  /**
   * GstClockClass::wait:
   * @clock: the #GstClock
   * @entry: the entry to wait on
   * @jitter: (out) (allow-none): a pointer that will contain the jitter
   *
   * Perform a blocking wait on the given #GstClockEntry and return
   * the jitter.
   *
   * Returns: the result of the blocking wait. #GST_CLOCK_EARLY will be returned
   * if the current clock time is past the time of @id, #GST_CLOCK_OK if
   * @id was scheduled in time. #GST_CLOCK_UNSCHEDULED if @id was
   * unscheduled with gst_clock_id_unschedule().
   */
  GstClockReturn        (*wait)                 (GstClock *clock, GstClockEntry *entry,
                                                 GstClockTimeDiff *jitter);

  /**
   * GstClockClass::wait_async:
   * @clock: the #GstClock
   * @entry: the entry to wait on
   *
   * Perform an asynchronous wait on the given #GstClockEntry.
   *
   * Returns: the result of the non blocking wait.
   */
  GstClockReturn        (*wait_async)           (GstClock *clock, GstClockEntry *entry);

  /**
   * GstClockClass::unschedule:
   * @clock: the #GstClock
   * @entry: the entry to unschedule
   *
   * Unblock a blocking or async wait operation.
   */
  void                  (*unschedule)           (GstClock *clock, GstClockEntry *entry);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType                   gst_clock_get_type              (void);

GST_API
GstClockTime            gst_clock_set_resolution        (GstClock *clock,
                                                         GstClockTime resolution);
GST_API
GstClockTime            gst_clock_get_resolution        (GstClock *clock);

GST_API
GstClockTime            gst_clock_get_time              (GstClock *clock);

GST_API
void                    gst_clock_set_calibration       (GstClock *clock, GstClockTime internal,
                                                         GstClockTime external,
                                                         GstClockTime rate_num,
                                                         GstClockTime rate_denom);
GST_API
void                    gst_clock_get_calibration       (GstClock *clock, GstClockTime *internal,
                                                         GstClockTime *external,
                                                         GstClockTime *rate_num,
                                                         GstClockTime *rate_denom);

/* master/slave clocks */

GST_API
gboolean                gst_clock_set_master            (GstClock *clock, GstClock *master);

GST_API
GstClock*               gst_clock_get_master            (GstClock *clock);

GST_API
void                    gst_clock_set_timeout           (GstClock *clock,
                                                         GstClockTime timeout);
GST_API
GstClockTime            gst_clock_get_timeout           (GstClock *clock);

GST_API
gboolean                gst_clock_add_observation       (GstClock *clock, GstClockTime slave,
                                                         GstClockTime master, gdouble *r_squared);
GST_API
gboolean                gst_clock_add_observation_unapplied (GstClock *clock, GstClockTime slave,
                                                         GstClockTime master, gdouble *r_squared,
                                                         GstClockTime *internal,
                                                         GstClockTime *external,
                                                         GstClockTime *rate_num,
                                                         GstClockTime *rate_denom);

/* getting and adjusting internal/external time */

GST_API
GstClockTime            gst_clock_get_internal_time     (GstClock *clock);

GST_API
GstClockTime            gst_clock_adjust_unlocked       (GstClock *clock, GstClockTime internal);

GST_API
GstClockTime            gst_clock_adjust_with_calibration (GstClock *clock,
                                                         GstClockTime internal_target,
                                                         GstClockTime cinternal,
                                                         GstClockTime cexternal,
                                                         GstClockTime cnum,
                                                         GstClockTime cdenom);
GST_API
GstClockTime            gst_clock_unadjust_with_calibration (GstClock *clock,
                                                         GstClockTime external_target,
                                                         GstClockTime cinternal,
                                                         GstClockTime cexternal,
                                                         GstClockTime cnum,
                                                         GstClockTime cdenom);
GST_API
GstClockTime            gst_clock_unadjust_unlocked     (GstClock * clock, GstClockTime external);

/* waiting for, signalling and checking for synchronization */

GST_API
gboolean                gst_clock_wait_for_sync         (GstClock * clock, GstClockTime timeout);

GST_API
gboolean                gst_clock_is_synced             (GstClock * clock);

/* to be used by subclasses only */

GST_API
void                    gst_clock_set_synced            (GstClock * clock, gboolean synced);

/* creating IDs that can be used to get notifications */

GST_API
GstClockID              gst_clock_new_single_shot_id    (GstClock *clock,
                                                         GstClockTime time);
GST_API
GstClockID              gst_clock_new_periodic_id       (GstClock *clock,
                                                         GstClockTime start_time,
                                                         GstClockTime interval);

/* reference counting */

GST_API
GstClockID              gst_clock_id_ref                (GstClockID id);

GST_API
void                    gst_clock_id_unref              (GstClockID id);

/* operations on IDs */

GST_API
gint                    gst_clock_id_compare_func       (gconstpointer id1, gconstpointer id2);

GST_API
GstClock *              gst_clock_id_get_clock          (GstClockID id);

GST_API
gboolean                gst_clock_id_uses_clock         (GstClockID id, GstClock * clock);

GST_API
GstClockTime            gst_clock_id_get_time           (GstClockID id);

GST_API
GstClockReturn          gst_clock_id_wait               (GstClockID id,
                                                         GstClockTimeDiff *jitter);
GST_API
GstClockReturn          gst_clock_id_wait_async         (GstClockID id,
                                                         GstClockCallback func,
                                                         gpointer user_data,
                                                         GDestroyNotify destroy_data);
GST_API
void                    gst_clock_id_unschedule         (GstClockID id);

GST_API
gboolean                gst_clock_single_shot_id_reinit (GstClock * clock,
                                                         GstClockID id,
                                                         GstClockTime time);
GST_API
gboolean                gst_clock_periodic_id_reinit    (GstClock * clock,
                                                         GstClockID id,
                                                         GstClockTime start_time,
                                                         GstClockTime interval);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstClock, gst_object_unref)
G_DEFINE_AUTO_CLEANUP_FREE_FUNC(GstClockID, gst_clock_id_unref, 0)

G_END_DECLS

#endif /* __GST_CLOCK_H__ */
