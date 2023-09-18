/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstaudioclock.h: Clock for use by audio plugins
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

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef __GST_AUDIO_CLOCK_H__
#define __GST_AUDIO_CLOCK_H__

#include <gst/gst.h>
#include <gst/gstsystemclock.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_CLOCK \
  (gst_audio_clock_get_type())
#define GST_AUDIO_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_CLOCK,GstAudioClock))
#define GST_AUDIO_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_CLOCK,GstAudioClockClass))
#define GST_IS_AUDIO_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_CLOCK))
#define GST_IS_AUDIO_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_CLOCK))
#define GST_AUDIO_CLOCK_CAST(obj) \
  ((GstAudioClock*)(obj))

typedef struct _GstAudioClock GstAudioClock;
typedef struct _GstAudioClockClass GstAudioClockClass;

/**
 * GstAudioClockGetTimeFunc:
 * @clock: the #GstAudioClock
 * @user_data: user data
 *
 * This function will be called whenever the current clock time needs to be
 * calculated. If this function returns #GST_CLOCK_TIME_NONE, the last reported
 * time will be returned by the clock.
 *
 * Returns: the current time or #GST_CLOCK_TIME_NONE if the previous time should
 * be used.
 */
typedef GstClockTime (*GstAudioClockGetTimeFunc) (GstClock *clock, gpointer user_data);

/**
 * GstAudioClock:
 *
 * Opaque #GstAudioClock.
 */
struct _GstAudioClock {
  GstSystemClock clock;

  /*< protected >*/
  GstAudioClockGetTimeFunc func;
  gpointer                 user_data;
  GDestroyNotify           destroy_notify;

  /*< private >*/
  GstClockTime             last_time;
  GstClockTimeDiff         time_offset;

  gpointer _gst_reserved[GST_PADDING];
};

struct _GstAudioClockClass {
  GstSystemClockClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType           gst_audio_clock_get_type        (void);

GST_AUDIO_API
GstClock*       gst_audio_clock_new             (const gchar *name, GstAudioClockGetTimeFunc func,
                                                 gpointer user_data, GDestroyNotify destroy_notify);

GST_AUDIO_API
void            gst_audio_clock_reset           (GstAudioClock *clock, GstClockTime time);

GST_AUDIO_API
GstClockTime    gst_audio_clock_get_time        (GstAudioClock * clock);

GST_AUDIO_API
GstClockTime    gst_audio_clock_adjust          (GstAudioClock * clock, GstClockTime time);

GST_AUDIO_API
void            gst_audio_clock_invalidate      (GstAudioClock * clock);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioClock, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_CLOCK_H__ */
