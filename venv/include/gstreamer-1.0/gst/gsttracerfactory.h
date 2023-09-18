/* GStreamer
 * Copyright (C) 2013 Stefan Sauer <ensonic@users.sf.net>
 *
 * gsttracerfactory.h: tracing subsystem
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

#ifndef __GST_TRACER_FACTORY_H__
#define __GST_TRACER_FACTORY_H__

#include <gst/gstcaps.h>
#include <gst/gstplugin.h>
#include <gst/gstpluginfeature.h>

G_BEGIN_DECLS

#define GST_TYPE_TRACER_FACTORY                 (gst_tracer_factory_get_type())
#define GST_TRACER_FACTORY(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_TRACER_FACTORY, GstTracerFactory))
#define GST_IS_TRACER_FACTORY(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_TRACER_FACTORY))
#define GST_TRACER_FACTORY_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_TRACER_FACTORY, GstTracerFactoryClass))
#define GST_IS_TRACER_FACTORY_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_TRACER_FACTORY))
#define GST_TRACER_FACTORY_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_TRACER_FACTORY, GstTracerFactoryClass))
#define GST_TRACER_FACTORY_CAST(obj)            ((GstTracerFactory *)(obj))

/**
 * GstTracerFactory:
 *
 * Opaque object that stores information about a tracer function.
 *
 * Since: 1.8
 */
typedef struct _GstTracerFactory GstTracerFactory;
typedef struct _GstTracerFactoryClass GstTracerFactoryClass;

/* tracering interface */

GST_API
GType           gst_tracer_factory_get_type          (void);

GST_API
GList *         gst_tracer_factory_get_list          (void);

GST_API
GType           gst_tracer_factory_get_tracer_type   (GstTracerFactory * factory);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTracerFactory, gst_object_unref)

G_END_DECLS

#endif /* __GST_TRACER_FACTORY_H__ */
