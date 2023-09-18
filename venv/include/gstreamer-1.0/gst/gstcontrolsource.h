/* GStreamer
 *
 * Copyright (C) 2007 Sebastian Dröge <slomo@circular-chaos.org>
 *
 * gstcontrolsource.h: Interface declaration for control sources
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

#ifndef __GST_CONTROL_SOURCE_H__
#define __GST_CONTROL_SOURCE_H__

#include <gst/gstconfig.h>

#include <glib-object.h>

#include <gst/gstclock.h>

G_BEGIN_DECLS

#define GST_TYPE_CONTROL_SOURCE \
  (gst_control_source_get_type())
#define GST_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_CONTROL_SOURCE,GstControlSource))
#define GST_CONTROL_SOURCE_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_CONTROL_SOURCE,GstControlSourceClass))
#define GST_IS_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_CONTROL_SOURCE))
#define GST_IS_CONTROL_SOURCE_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_CONTROL_SOURCE))
#define GST_CONTROL_SOURCE_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_CONTOL_SOURCE, GstControlSourceClass))

typedef struct _GstControlSource GstControlSource;
typedef struct _GstControlSourceClass GstControlSourceClass;
typedef struct _GstTimedValue GstTimedValue;
typedef struct _GstValueArray GstValueArray;

/**
 * GstTimedValue:
 * @timestamp: timestamp of the value change
 * @value: the corresponding value
 *
 * Structure for storing a timestamp and a value.
 */
struct _GstTimedValue
{
  GstClockTime timestamp;
  gdouble      value;
};

/**
 * GstControlSourceGetValue:
 * @self: the #GstControlSource instance
 * @timestamp: timestamp for which a value should be calculated
 * @value: a value which will be set to the result.
 *
 * Function for returning a value for a given timestamp.
 *
 * Returns: %TRUE if the value was successfully calculated.
 *
 */
typedef gboolean (* GstControlSourceGetValue) (GstControlSource *self, 
    GstClockTime timestamp, gdouble *value);

/**
 * GstControlSourceGetValueArray:
 * @self: the #GstControlSource instance
 * @timestamp: timestamp for which a value should be calculated
 * @interval: the time spacing between subsequent values
 * @n_values: the number of values
 * @values: array to put control-values in
 *
 * Function for returning an array of values starting at a given timestamp.
 *
 * Returns: %TRUE if the values were successfully calculated.
 *
 */
typedef gboolean (* GstControlSourceGetValueArray) (GstControlSource *self, 
    GstClockTime timestamp, GstClockTime interval, guint n_values, gdouble *values);

/**
 * GstControlSource:
 * @parent: the parent structure
 * @get_value: Function for returning a value for a given timestamp
 * @get_value_array: Function for returning a values array for a given timestamp
 *
 * The instance structure of #GstControlSource.
 */
struct _GstControlSource {
  GstObject parent;

  /*< public >*/
  GstControlSourceGetValue get_value;             /* Returns the value for a property at a given timestamp */
  GstControlSourceGetValueArray get_value_array;  /* Returns values for a property in a given timespan */

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstControlSourceClass:
 * @parent_class: Parent class
 *
 * The class structure of #GstControlSource.
 */

struct _GstControlSourceClass
{
  GstObjectClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType          gst_control_source_get_type (void);

/* Functions */

GST_API
gboolean       gst_control_source_get_value             (GstControlSource *self, GstClockTime timestamp,
                                                         gdouble *value);
GST_API
gboolean       gst_control_source_get_value_array       (GstControlSource *self, GstClockTime timestamp,
                                                         GstClockTime interval, guint n_values,
                                                         gdouble *values);
G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstControlSource, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstValueArray, gst_object_unref)

G_END_DECLS

#endif /* __GST_CONTROL_SOURCE_H__ */
