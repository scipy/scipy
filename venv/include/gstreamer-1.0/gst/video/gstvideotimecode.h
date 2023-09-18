/* GStreamer
 * Copyright (C) <2016> Vivia Nikolaidou <vivia@toolsonair.com>
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

#ifndef __GST_VIDEO_TIME_CODE_H__
#define __GST_VIDEO_TIME_CODE_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

typedef struct _GstVideoTimeCodeConfig GstVideoTimeCodeConfig;
typedef struct _GstVideoTimeCode GstVideoTimeCode;
typedef struct _GstVideoTimeCodeInterval GstVideoTimeCodeInterval;

/**
 * GstVideoTimeCodeFlags:
 * @GST_VIDEO_TIME_CODE_FLAGS_NONE: No flags
 * @GST_VIDEO_TIME_CODE_FLAGS_DROP_FRAME: Whether we have drop frame rate
 * @GST_VIDEO_TIME_CODE_FLAGS_INTERLACED: Whether we have interlaced video
 *
 * Flags related to the time code information.
 * For drop frame, only 30000/1001 and 60000/1001 frame rates are supported.
 *
 * Since: 1.10
 */
typedef enum
{
  GST_VIDEO_TIME_CODE_FLAGS_NONE = 0,
  GST_VIDEO_TIME_CODE_FLAGS_DROP_FRAME = (1<<0),
  GST_VIDEO_TIME_CODE_FLAGS_INTERLACED = (1<<1)
  /* Not supported yet:
   * GST_VIDEO_TIME_CODE_ALLOW_MORE_THAN_24H = (1<<2)
   * GST_VIDEO_TIME_CODE_ALLOW_NEGATIVE = (1<<3)
   */
} GstVideoTimeCodeFlags;

/**
 * GstVideoTimeCodeConfig:
 * @fps_n: Numerator of the frame rate
 * @fps_d: Denominator of the frame rate
 * @flags: the corresponding #GstVideoTimeCodeFlags
 * @latest_daily_jam: The latest daily jam information, if present, or NULL
 *
 * Supported frame rates: 30000/1001, 60000/1001 (both with and without drop
 * frame), and integer frame rates e.g. 25/1, 30/1, 50/1, 60/1.
 *
 * The configuration of the time code.
 *
 * Since: 1.10
 */
struct _GstVideoTimeCodeConfig {
  guint fps_n;
  guint fps_d;
  GstVideoTimeCodeFlags flags;
  GDateTime *latest_daily_jam;
};

/**
 * GstVideoTimeCode:
 * @hours: the hours field of #GstVideoTimeCode
 * @minutes: the minutes field of #GstVideoTimeCode
 * @seconds: the seconds field of #GstVideoTimeCode
 * @frames: the frames field of #GstVideoTimeCode
 * @field_count: Interlaced video field count
 * @config: the corresponding #GstVideoTimeCodeConfig
 *
 * @field_count must be 0 for progressive video and 1 or 2 for interlaced.
 *
 * A representation of a SMPTE time code.
 *
 * @hours must be positive and less than 24. Will wrap around otherwise.
 * @minutes and @seconds must be positive and less than 60.
 * @frames must be less than or equal to @config.fps_n / @config.fps_d
 * These values are *NOT* automatically normalized.
 *
 * Since: 1.10
 */
struct _GstVideoTimeCode {
  GstVideoTimeCodeConfig config;

  guint hours;
  guint minutes;
  guint seconds;
  guint frames;
  guint field_count;
};

/**
 * GstVideoTimeCodeInterval:
 * @hours: the hours field of #GstVideoTimeCodeInterval
 * @minutes: the minutes field of #GstVideoTimeCodeInterval
 * @seconds: the seconds field of #GstVideoTimeCodeInterval
 * @frames: the frames field of #GstVideoTimeCodeInterval
 *
 * A representation of a difference between two #GstVideoTimeCode instances.
 * Will not necessarily correspond to a real timecode (e.g. 00:00:10;00)
 *
 * Since: 1.12
 */
struct _GstVideoTimeCodeInterval {
  guint hours;
  guint minutes;
  guint seconds;
  guint frames;
};

#define GST_VIDEO_TIME_CODE_INIT { {0, 0, 0, NULL}, 0, 0, 0, 0, 0 }

#define GST_TYPE_VIDEO_TIME_CODE (gst_video_time_code_get_type())
GST_VIDEO_API
GType gst_video_time_code_get_type (void);

GST_VIDEO_API
GstVideoTimeCode * gst_video_time_code_new          (guint                    fps_n,
                                                     guint                    fps_d,
                                                     GDateTime              * latest_daily_jam,
                                                     GstVideoTimeCodeFlags    flags,
                                                     guint                    hours,
                                                     guint                    minutes,
                                                     guint                    seconds,
                                                     guint                    frames,
                                                     guint                    field_count);

GST_VIDEO_API
GstVideoTimeCode * gst_video_time_code_new_empty    (void);

GST_VIDEO_API
GstVideoTimeCode * gst_video_time_code_new_from_string    (const gchar * tc_str);

GST_VIDEO_DEPRECATED_FOR(gst_video_time_code_new_from_date_time_full)
GstVideoTimeCode * gst_video_time_code_new_from_date_time (guint                    fps_n,
                                                           guint                    fps_d,
                                                           GDateTime              * dt,
                                                           GstVideoTimeCodeFlags    flags,
                                                           guint                    field_count);

GST_VIDEO_API
GstVideoTimeCode * gst_video_time_code_new_from_date_time_full (guint                    fps_n,
                                                                guint                    fps_d,
                                                                GDateTime              * dt,
                                                                GstVideoTimeCodeFlags    flags,
                                                                guint                    field_count);

GST_VIDEO_API
void gst_video_time_code_free                       (GstVideoTimeCode       * tc);

GST_VIDEO_API
GstVideoTimeCode * gst_video_time_code_copy         (const GstVideoTimeCode * tc);

GST_VIDEO_API
void gst_video_time_code_init                       (GstVideoTimeCode       * tc,
                                                     guint                    fps_n,
                                                     guint                    fps_d,
                                                     GDateTime              * latest_daily_jam,
                                                     GstVideoTimeCodeFlags    flags,
                                                     guint                    hours,
                                                     guint                    minutes,
                                                     guint                    seconds,
                                                     guint                    frames,
                                                     guint                    field_count);

GST_VIDEO_DEPRECATED_FOR(gst_video_time_code_init_from_date_time_full)
void gst_video_time_code_init_from_date_time        (GstVideoTimeCode       * tc,
                                                     guint                    fps_n,
                                                     guint                    fps_d,
                                                     GDateTime              * dt,
                                                     GstVideoTimeCodeFlags    flags,
                                                     guint                    field_count);
GST_VIDEO_API
gboolean gst_video_time_code_init_from_date_time_full (GstVideoTimeCode       * tc,
                                                       guint                    fps_n,
                                                       guint                    fps_d,
                                                       GDateTime              * dt,
                                                       GstVideoTimeCodeFlags    flags,
                                                       guint                    field_count);

GST_VIDEO_API
void gst_video_time_code_clear                      (GstVideoTimeCode       * tc);

GST_VIDEO_API
gboolean gst_video_time_code_is_valid               (const GstVideoTimeCode * tc);

GST_VIDEO_API
gint gst_video_time_code_compare                    (const GstVideoTimeCode * tc1,
                                                     const GstVideoTimeCode * tc2);

GST_VIDEO_API
void gst_video_time_code_increment_frame            (GstVideoTimeCode       * tc);

GST_VIDEO_API
void gst_video_time_code_add_frames                 (GstVideoTimeCode       * tc,
                                                     gint64                   frames);

GST_VIDEO_API
gchar *gst_video_time_code_to_string                (const GstVideoTimeCode * tc);

GST_VIDEO_API
GDateTime *gst_video_time_code_to_date_time         (const GstVideoTimeCode * tc);

GST_VIDEO_API
guint64 gst_video_time_code_nsec_since_daily_jam    (const GstVideoTimeCode * tc);

GST_VIDEO_API
guint64 gst_video_time_code_frames_since_daily_jam  (const GstVideoTimeCode * tc);

GST_VIDEO_API
GstVideoTimeCode * gst_video_time_code_add_interval (const GstVideoTimeCode * tc, const GstVideoTimeCodeInterval * tc_inter);

#define GST_TYPE_VIDEO_TIME_CODE_INTERVAL (gst_video_time_code_interval_get_type())
GST_VIDEO_API
GType gst_video_time_code_interval_get_type (void);

GST_VIDEO_API
GstVideoTimeCodeInterval * gst_video_time_code_interval_new  (guint                    hours,
                                                     guint                    minutes,
                                                     guint                    seconds,
                                                     guint                    frames);

GST_VIDEO_API
GstVideoTimeCodeInterval * gst_video_time_code_interval_new_from_string    (const gchar * tc_inter_str);

GST_VIDEO_API
void gst_video_time_code_interval_free                   (GstVideoTimeCodeInterval       * tc);

GST_VIDEO_API
GstVideoTimeCodeInterval * gst_video_time_code_interval_copy (const GstVideoTimeCodeInterval * tc);

GST_VIDEO_API
void gst_video_time_code_interval_init                   (GstVideoTimeCodeInterval       * tc,
                                                     guint                    hours,
                                                     guint                    minutes,
                                                     guint                    seconds,
                                                     guint                    frames);

GST_VIDEO_API
void gst_video_time_code_interval_clear                  (GstVideoTimeCodeInterval       * tc);

G_END_DECLS

#endif /* __GST_VIDEO_TIME_CODE_H__ */
