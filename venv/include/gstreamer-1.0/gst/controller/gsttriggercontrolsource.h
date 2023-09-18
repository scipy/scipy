/* GStreamer
 *
 * Copyright (C) 2007 Sebastian Dröge <slomo@circular-chaos.org>
 *               2011 Stefan Sauer <ensonic@users.sf.net>
 *
 * gsttriggercontrolsource.h: Control source that provides some values at time-
 *                            stamps
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


#ifndef __GST_TRIGGER_CONTROL_SOURCE_H__
#define __GST_TRIGGER_CONTROL_SOURCE_H__

#include <glib-object.h>
#include <gst/gst.h>

#include <gst/controller/gsttimedvaluecontrolsource.h>

G_BEGIN_DECLS

#define GST_TYPE_TRIGGER_CONTROL_SOURCE \
  (gst_trigger_control_source_get_type ())
#define GST_TRIGGER_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_TRIGGER_CONTROL_SOURCE, GstTriggerControlSource))
#define GST_TRIGGER_CONTROL_SOURCE_CLASS(vtable) \
  (G_TYPE_CHECK_CLASS_CAST ((vtable), GST_TYPE_TRIGGER_CONTROL_SOURCE, GstTriggerControlSourceClass))
#define GST_IS_TRIGGER_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_TRIGGER_CONTROL_SOURCE))
#define GST_IS_TRIGGER_CONTROL_SOURCE_CLASS(vtable) \
  (G_TYPE_CHECK_CLASS_TYPE ((vtable), GST_TYPE_TRIGGER_CONTROL_SOURCE))
#define GST_TRIGGER_CONTROL_SOURCE_GET_CLASS(inst) \
  (G_TYPE_INSTANCE_GET_CLASS ((inst), GST_TYPE_TRIGGER_CONTROL_SOURCE, GstTriggerControlSourceClass))

#define GST_TYPE_TRIGGER_WAVEFORM (gst_trigger_waveform_get_type ())

typedef struct _GstTriggerControlSource GstTriggerControlSource;
typedef struct _GstTriggerControlSourceClass GstTriggerControlSourceClass;
typedef struct _GstTriggerControlSourcePrivate GstTriggerControlSourcePrivate;

/**
 * GstTriggerControlSource:
 *
 * The instance structure of #GstControlSource.
 */
struct _GstTriggerControlSource {
  GstTimedValueControlSource parent;

  /*< private >*/
  GstTriggerControlSourcePrivate *priv;
  gpointer _gst_reserved[GST_PADDING];
};

struct _GstTriggerControlSourceClass {
  GstTimedValueControlSourceClass parent_class;
  
  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_CONTROLLER_API
GType gst_trigger_control_source_get_type (void);

/* Functions */

GST_CONTROLLER_API
GstControlSource *gst_trigger_control_source_new (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTriggerControlSource, gst_object_unref)

G_END_DECLS

#endif /* __GST_TRIGGER_CONTROL_SOURCE_H__ */
