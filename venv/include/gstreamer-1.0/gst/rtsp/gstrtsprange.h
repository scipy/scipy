/* GStreamer
 * Copyright (C) <2005,2006> Wim Taymans <wim@fluendo.com>
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
/*
 * Unless otherwise indicated, Source Code is licensed under MIT license.
 * See further explanation attached in License Statement (distributed in the file
 * LICENSE).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __GST_RTSP_RANGE_H__
#define __GST_RTSP_RANGE_H__

#include <glib.h>
#include <gst/gst.h>

#include <gst/rtsp/gstrtspdefs.h>

G_BEGIN_DECLS

/**
 * GstRTSPRangeUnit:
 * @GST_RTSP_RANGE_SMPTE: SMPTE timecode
 * @GST_RTSP_RANGE_SMPTE_30_DROP: 29.97 frames per second
 * @GST_RTSP_RANGE_SMPTE_25: 25 frames per second
 * @GST_RTSP_RANGE_NPT: Normal play time
 * @GST_RTSP_RANGE_CLOCK: Absolute time expressed as ISO 8601 timestamps
 *
 * Different possible time range units.
 */
typedef enum
{
  GST_RTSP_RANGE_SMPTE,
  GST_RTSP_RANGE_SMPTE_30_DROP,
  GST_RTSP_RANGE_SMPTE_25,
  GST_RTSP_RANGE_NPT,
  GST_RTSP_RANGE_CLOCK
} GstRTSPRangeUnit;

typedef struct _GstRTSPTimeRange GstRTSPTimeRange;
typedef struct _GstRTSPTime GstRTSPTime;
typedef struct _GstRTSPTime2 GstRTSPTime2;

/**
 * GstRTSPTimeType:
 * @GST_RTSP_TIME_SECONDS: seconds
 * @GST_RTSP_TIME_NOW: now
 * @GST_RTSP_TIME_END: end
 * @GST_RTSP_TIME_FRAMES: frames and subframes
 * @GST_RTSP_TIME_UTC: UTC time
 *
 * Possible time types.
 */
typedef enum {
  GST_RTSP_TIME_SECONDS,
  GST_RTSP_TIME_NOW,
  GST_RTSP_TIME_END,
  GST_RTSP_TIME_FRAMES,
  GST_RTSP_TIME_UTC
} GstRTSPTimeType;

/**
 * GstRTSPTime:
 * @type: the time of the time
 * @seconds: seconds when @type is GST_RTSP_TIME_SECONDS,
 *           GST_RTSP_TIME_UTC and GST_RTSP_TIME_FRAMES
 *
 * A time indication.
 */
struct _GstRTSPTime {
  GstRTSPTimeType type;
  gdouble         seconds;
};

/**
 * GstRTSPTime2:
 * @frames: frames and subframes when type in GstRTSPTime is
 *          GST_RTSP_TIME_FRAMES
 * @year: year when type is GST_RTSP_TIME_UTC
 * @month: month when type is GST_RTSP_TIME_UTC
 * @day: day when type is GST_RTSP_TIME_UTC
 *
 * Extra fields for a time indication.
 *
 * Since: 1.2
 */
struct _GstRTSPTime2 {
  gdouble         frames;
  guint           year;
  guint           month;
  guint           day;
};

/**
 * GstRTSPTimeRange:
 * @unit: the time units used
 * @min: the minimum interval
 * @max: the maximum interval
 * @min2: extra fields in the minimum interval (Since: 1.2)
 * @max2: extra fields in the maximum interval (Since: 1.2)
 *
 * A time range.
 */
struct _GstRTSPTimeRange {
  GstRTSPRangeUnit unit;

  GstRTSPTime  min;
  GstRTSPTime  max;
  GstRTSPTime2 min2;
  GstRTSPTime2 max2;
};

GST_RTSP_API
GstRTSPResult   gst_rtsp_range_parse        (const gchar *rangestr, GstRTSPTimeRange **range);

GST_RTSP_API
gchar *         gst_rtsp_range_to_string    (const GstRTSPTimeRange *range);

GST_RTSP_API
void            gst_rtsp_range_free         (GstRTSPTimeRange *range);

GST_RTSP_API
gboolean        gst_rtsp_range_get_times     (const GstRTSPTimeRange *range,
                                              GstClockTime *min, GstClockTime *max);

GST_RTSP_API
gboolean        gst_rtsp_range_convert_units (GstRTSPTimeRange * range,
                                              GstRTSPRangeUnit unit);

G_END_DECLS

#endif /* __GST_RTSP_RANGE_H__ */
