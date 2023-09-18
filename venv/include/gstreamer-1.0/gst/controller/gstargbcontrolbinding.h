/* GStreamer
 *
 * Copyright (C) 2011 Stefan Sauer <ensonic@users.sf.net>
 *
 * gstargbcontrolbinding.h: Attachment for multiple control sources to gargb
 *                            properties
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

#ifndef __GST_ARGB_CONTROL_BINDING_H__
#define __GST_ARGB_CONTROL_BINDING_H__

#include <gst/gstconfig.h>

#include <glib-object.h>

#include <gst/gstcontrolsource.h>
#include <gst/controller/controller-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_ARGB_CONTROL_BINDING \
  (gst_argb_control_binding_get_type())
#define GST_ARGB_CONTROL_BINDING(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_ARGB_CONTROL_BINDING,GstARGBControlBinding))
#define GST_ARGB_CONTROL_BINDING_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_ARGB_CONTROL_BINDING,GstARGBControlBindingClass))
#define GST_IS_ARGB_CONTROL_BINDING(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_ARGB_CONTROL_BINDING))
#define GST_IS_ARGB_CONTROL_BINDING_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_ARGB_CONTROL_BINDING))
#define GST_ARGB_CONTROL_BINDING_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_CONTOL_SOURCE, GstARGBControlBindingClass))

typedef struct _GstARGBControlBinding GstARGBControlBinding;
typedef struct _GstARGBControlBindingClass GstARGBControlBindingClass;

/**
 * GstARGBControlBinding:
 * @name: name of the property of this binding
 *
 * The instance structure of #GstARGBControlBinding.
 */
struct _GstARGBControlBinding {
  GstControlBinding parent;
  
  /*< private >*/
  GstControlSource *cs_a;       /* GstControlSources for this property */
  GstControlSource *cs_r;
  GstControlSource *cs_g;
  GstControlSource *cs_b;

  GValue cur_value;
  guint32 last_value;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstARGBControlBindingClass:
 * @parent_class: Parent class
 * @convert: Class method to convert control-values
 *
 * The class structure of #GstARGBControlBinding.
 */

struct _GstARGBControlBindingClass
{
  GstControlBindingClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_CONTROLLER_API
GType gst_argb_control_binding_get_type (void);

/* Functions */

GST_CONTROLLER_API
GstControlBinding * gst_argb_control_binding_new   (GstObject * object, const gchar * property_name,
                                                            GstControlSource * cs_a, GstControlSource * cs_r,
                                                            GstControlSource * cs_g, GstControlSource * cs_b);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstARGBControlBinding, gst_object_unref)

G_END_DECLS

#endif /* __GST_ARGB_CONTROL_BINDING_H__ */
