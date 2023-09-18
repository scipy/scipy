/* GStreamer
 * Copyright (C) 2013 Stefan Sauer <ensonic@users.sf.net>
 *
 * gsttracer.h: tracer base class
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

#ifndef __GST_TRACER_H__
#define __GST_TRACER_H__

#include <glib.h>
#include <glib-object.h>
#include <gst/gstobject.h>
#include <gst/gstconfig.h>

G_BEGIN_DECLS

typedef struct _GstTracer GstTracer;
typedef struct _GstTracerPrivate GstTracerPrivate;
typedef struct _GstTracerClass GstTracerClass;

#define GST_TYPE_TRACER            (gst_tracer_get_type())
#define GST_TRACER(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_TRACER,GstTracer))
#define GST_TRACER_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_TRACER,GstTracerClass))
#define GST_IS_TRACER(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_TRACER))
#define GST_IS_TRACER_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_TRACER))
#define GST_TRACER_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_TRACER,GstTracerClass))
#define GST_TRACER_CAST(obj)       ((GstTracer *)(obj))

/**
 * GstTracer:
 *
 * The opaque GstTracer instance structure
 */
struct _GstTracer {
  GstObject        parent;
  /*< private >*/
  GstTracerPrivate *priv;
  gpointer _gst_reserved[GST_PADDING];
};

struct _GstTracerClass {
  GstObjectClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType gst_tracer_get_type          (void);

GST_API
void gst_tracing_register_hook (GstTracer *tracer, const gchar *detail,
  GCallback func);

/* tracing modules */

GST_API
gboolean gst_tracer_register (GstPlugin * plugin, const gchar * name, GType type);

GST_API
GList* gst_tracing_get_active_tracers (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTracer, gst_object_unref)

G_END_DECLS

#endif /* __GST_TRACER_H__ */

