/* GStreamer
 * Copyright (C) 2016 Stefan Sauer <ensonic@users.sf.net>
 *
 * gsttracerrecord.h: tracer log record class
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

#ifndef __GST_TRACER_RECORD_H__
#define __GST_TRACER_RECORD_H__

#include <gst/gstobject.h>

G_BEGIN_DECLS

/**
 * GstTracerRecord:
 *
 * The opaque GstTracerRecord instance structure
 *
 * Since: 1.8
 */
typedef struct _GstTracerRecord GstTracerRecord;
typedef struct _GstTracerRecordClass GstTracerRecordClass;

#define GST_TYPE_TRACER_RECORD            (gst_tracer_record_get_type())
#define GST_TRACER_RECORD(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_TRACER_RECORD,GstTracerRecord))
#define GST_TRACER_RECORD_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_TRACER_RECORD,GstTracerRecordClass))
#define GST_IS_TRACER_RECORD(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_TRACER_RECORD))
#define GST_IS_TRACER_RECORD_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_TRACER_RECORD))
#define GST_TRACER_RECORD_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_TRACER_RECORD,GstTracerRecordClass))
#define GST_TRACER_RECORD_CAST(obj)       ((GstTracerRecord *)(obj))

GST_API
GType gst_tracer_record_get_type          (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTracerRecord, gst_object_unref)

/**
 * GstTracerValueScope:
 * @GST_TRACER_VALUE_SCOPE_PROCESS: the value is related to the process
 * @GST_TRACER_VALUE_SCOPE_THREAD: the value is related to a thread
 * @GST_TRACER_VALUE_SCOPE_ELEMENT: the value is related to an #GstElement
 * @GST_TRACER_VALUE_SCOPE_PAD: the value is related to a #GstPad
 *
 * Tracing record will contain fields that contain a measured value or extra
 * meta-data. One such meta data are values that tell where a measurement was
 * taken. This enumerating declares to which scope such a meta data field
 * relates to. If it is e.g. %GST_TRACER_VALUE_SCOPE_PAD, then each of the log
 * events may contain values for different #GstPads.
 *
 * Since: 1.8
 */
typedef enum
{
  GST_TRACER_VALUE_SCOPE_PROCESS,
  GST_TRACER_VALUE_SCOPE_THREAD,
  GST_TRACER_VALUE_SCOPE_ELEMENT,
  GST_TRACER_VALUE_SCOPE_PAD
} GstTracerValueScope;

/**
 * GstTracerValueFlags:
 * @GST_TRACER_VALUE_FLAGS_NONE: no flags
 * @GST_TRACER_VALUE_FLAGS_OPTIONAL: the value is optional. When using this flag
 *   one need to have an additional boolean arg before this value in the
 *   var-args list passed to  gst_tracer_record_log().
 * @GST_TRACER_VALUE_FLAGS_AGGREGATED: the value is a combined figure, since the
 *   start of tracing. Examples are averages or timestamps.
 *
 * Flag that describe the value. These flags help applications processing the
 * logs to understand the values.
 */
typedef enum
{
  GST_TRACER_VALUE_FLAGS_NONE = 0,
  GST_TRACER_VALUE_FLAGS_OPTIONAL = (1 << 0),
  GST_TRACER_VALUE_FLAGS_AGGREGATED = (1 << 1),
} GstTracerValueFlags;

GST_API
GstTracerRecord * gst_tracer_record_new (const gchar * name, const gchar * firstfield, ...) G_GNUC_NULL_TERMINATED;

#ifndef GST_DISABLE_GST_DEBUG
GST_API
void              gst_tracer_record_log (GstTracerRecord *self, ...);
#else
#define gst_tracer_record_log(...) G_STMT_START {} G_STMT_END
#endif

G_END_DECLS

#endif /* __GST_TRACER_RECORD_H__ */
