/* GStreamer
 *
 * Copyright (C) 2007 Sebastian Dröge <slomo@circular-chaos.org>
 *
 * gstinterpolationcontrolsource.h: Control source that provides several
 *                                  interpolation methods
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

#ifndef __GST_INTERPOLATION_CONTROL_SOURCE_H__
#define __GST_INTERPOLATION_CONTROL_SOURCE_H__

#include <glib-object.h>
#include <gst/gst.h>

#include <gst/controller/gsttimedvaluecontrolsource.h>
#include <gst/controller/controller-enumtypes.h>

G_BEGIN_DECLS

#define GST_TYPE_INTERPOLATION_CONTROL_SOURCE \
  (gst_interpolation_control_source_get_type ())
#define GST_INTERPOLATION_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_INTERPOLATION_CONTROL_SOURCE, GstInterpolationControlSource))
#define GST_INTERPOLATION_CONTROL_SOURCE_CLASS(vtable) \
  (G_TYPE_CHECK_CLASS_CAST ((vtable), GST_TYPE_INTERPOLATION_CONTROL_SOURCE, GstInterpolationControlSourceClass))
#define GST_IS_INTERPOLATION_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_INTERPOLATION_CONTROL_SOURCE))
#define GST_IS_INTERPOLATION_CONTROL_SOURCE_CLASS(vtable) \
  (G_TYPE_CHECK_CLASS_TYPE ((vtable), GST_TYPE_INTERPOLATION_CONTROL_SOURCE))
#define GST_INTERPOLATION_CONTROL_SOURCE_GET_CLASS(inst) \
  (G_TYPE_INSTANCE_GET_CLASS ((inst), GST_TYPE_INTERPOLATION_CONTROL_SOURCE, GstInterpolationControlSourceClass))

typedef struct _GstInterpolationControlSource GstInterpolationControlSource;
typedef struct _GstInterpolationControlSourceClass GstInterpolationControlSourceClass;
typedef struct _GstInterpolationControlSourcePrivate GstInterpolationControlSourcePrivate;

/**
 * GstInterpolationMode:
 * @GST_INTERPOLATION_MODE_NONE: steps-like interpolation, default
 * @GST_INTERPOLATION_MODE_LINEAR: linear interpolation
 * @GST_INTERPOLATION_MODE_CUBIC: cubic interpolation (natural), may overshoot
 *   the min or max values set by the control point, but is more 'curvy'
 * @GST_INTERPOLATION_MODE_CUBIC_MONOTONIC: monotonic cubic interpolation, will not
 *   produce any values outside of the min-max range set by the control points
 *   (Since: 1.8)
 *
 * The various interpolation modes available.
 */
typedef enum
{
  GST_INTERPOLATION_MODE_NONE,
  GST_INTERPOLATION_MODE_LINEAR,
  GST_INTERPOLATION_MODE_CUBIC,
  GST_INTERPOLATION_MODE_CUBIC_MONOTONIC,
} GstInterpolationMode;

/**
 * GstInterpolationControlSource:
 *
 * The instance structure of #GstControlSource.
 */
struct _GstInterpolationControlSource {
  GstTimedValueControlSource parent;

  /*< private >*/
  GstInterpolationControlSourcePrivate *priv;
  gpointer _gst_reserved[GST_PADDING];
};

struct _GstInterpolationControlSourceClass {
  GstTimedValueControlSourceClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_CONTROLLER_API
GType gst_interpolation_control_source_get_type (void);

/* Functions */

GST_CONTROLLER_API
GstControlSource * gst_interpolation_control_source_new (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstInterpolationControlSource, gst_object_unref)

G_END_DECLS

#endif /* __GST_INTERPOLATION_CONTROL_SOURCE_H__ */
