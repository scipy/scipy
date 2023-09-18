/*
 * GStreamer
 * Copyright (C) 2016 Matthew Waters <matthew@centricular.com>
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

#ifndef __GST_PROXY_CONTROL_BINDING_H__
#define __GST_PROXY_CONTROL_BINDING_H__

#include <gst/gst.h>
#include <gst/controller/controller-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_PROXY_CONTROL_BINDING  (gst_proxy_control_binding_get_type())
#define GST_PROXY_CONTROL_BINDING(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_PROXY_CONTROL_BINDING,GstProxyControlBinding))
#define GST_PROXY_CONTROL_BINDING_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_PROXY_CONTROL_BINDING,GstProxyControlBindingClass))
#define GST_IS_PROXY_CONTROL_BINDING(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_PROXY_CONTROL_BINDING))
#define GST_IS_PROXY_CONTROL_BINDING_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_PROXY_CONTROL_BINDING))
#define GST_PROXY_CONTROL_BINDING_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_CONTOL_SOURCE, GstProxyControlBindingClass))

typedef struct _GstProxyControlBinding GstProxyControlBinding;
typedef struct _GstProxyControlBindingClass GstProxyControlBindingClass;

/**
 * GstProxyControlBinding:
 *
 * Opaque #GstProxyControlBinding struct
 */
struct _GstProxyControlBinding
{
  /* <private> */
  GstControlBinding parent;

  GWeakRef ref_object;
  gchar *property_name;

  gpointer _padding[GST_PADDING];
};

/**
 * GstProxyControlBindingClass:
 *
 * Opaque #GstProxyControlBindingClass struct
 */
struct _GstProxyControlBindingClass
{
  /* <private> */
  GstControlBindingClass parent_class;

  gpointer _padding[GST_PADDING];
};

GST_CONTROLLER_API
GType                   gst_proxy_control_binding_get_type (void);

GST_CONTROLLER_API
GstControlBinding *     gst_proxy_control_binding_new (GstObject * object,
                                                       const gchar * property_name,
                                                       GstObject * ref_object,
                                                       const gchar * ref_property_name);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstProxyControlBinding, gst_object_unref)
G_END_DECLS

#endif /* __GST_PROXY_CONTROL_BINDING_H__ */
