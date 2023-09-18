/* GStreamer
 * Copyright (C) 2005 Stefan Kost <ensonic@users.sf.net>
 *
 * gstchildproxy.h: interface header for multi child elements
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

#ifndef __GST_CHILD_PROXY_H__
#define __GST_CHILD_PROXY_H__

#include <glib-object.h>
#include <gst/gst.h>

G_BEGIN_DECLS


#define GST_TYPE_CHILD_PROXY               (gst_child_proxy_get_type ())
#define GST_CHILD_PROXY(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_CHILD_PROXY, GstChildProxy))
#define GST_IS_CHILD_PROXY(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_CHILD_PROXY))
#define GST_CHILD_PROXY_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), GST_TYPE_CHILD_PROXY, GstChildProxyInterface))

/**
 * GstChildProxy:
 *
 * Opaque #GstChildProxy data structure.
 */
typedef struct _GstChildProxy GstChildProxy;    /* dummy object */
typedef struct _GstChildProxyInterface GstChildProxyInterface;

/**
 * GstChildProxyInterface:
 * @parent: parent interface type.
 *
 * #GstChildProxy interface.
 */
struct _GstChildProxyInterface
{
  GTypeInterface parent;

  /* methods */

  /**
   * GstChildProxyInterface.get_child_by_name:
   * @parent: the #GstChildProxy
   * @name: the name of the child to fetch
   *
   * Fetch a child object by name
   *
   * Returns: (transfer full) (nullable): the child object
   */
  GObject * (*get_child_by_name)  (GstChildProxy * parent, const gchar * name);

  /**
   * GstChildProxyInterface.get_child_by_index:
   * @parent: the #GstChildProxy
   * @index: the index of the child to fetch
   *
   * Fetch a child object by index
   *
   * Returns: (transfer full) (nullable): the child object
   */
  GObject * (*get_child_by_index) (GstChildProxy * parent, guint index);

  /**
   * GstChildProxyInterface.get_children_count:
   * @parent: the #GstChildProxy
   *
   * Get the number of children in @parent
   *
   * Returns: the number of children
   */
  guint     (*get_children_count) (GstChildProxy * parent);

  /*< private >*/
  /* signals */

  /**
   * GstChildProxyInterface.child_added:
   * @parent: the #GstChildProxy
   * @child: the child object
   * @name: the name of the child object
   *
   * Called when @child is added to @parent
   */
  void      (*child_added)        (GstChildProxy * parent, GObject * child, const gchar * name);

  /**
   * GstChildProxyInterface.child_removed:
   * @parent: the #GstChildProxy
   * @child: the child object
   * @name: the name of the child object
   *
   * Called when @child is removed from @parent
   */
  void      (*child_removed)      (GstChildProxy * parent, GObject * child, const gchar * name);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType     gst_child_proxy_get_type (void);

GST_API
GObject * gst_child_proxy_get_child_by_name  (GstChildProxy * parent, const gchar * name);

GST_API
GObject * gst_child_proxy_get_child_by_name_recurse (GstChildProxy * child_proxy,
                                                     const gchar *name);

GST_API
guint     gst_child_proxy_get_children_count (GstChildProxy * parent);

GST_API
GObject * gst_child_proxy_get_child_by_index (GstChildProxy * parent, guint index);

GST_API
gboolean  gst_child_proxy_lookup             (GstChildProxy *object, const gchar *name,
                                              GObject **target, GParamSpec **pspec);
GST_API
void      gst_child_proxy_get_property       (GstChildProxy * object, const gchar *name,
                                              GValue *value);
GST_API
void      gst_child_proxy_get_valist         (GstChildProxy * object,
                                              const gchar * first_property_name,
                                              va_list var_args);
GST_API
void      gst_child_proxy_get                (GstChildProxy * object,
                                              const gchar * first_property_name,
                                              ...) G_GNUC_NULL_TERMINATED;
GST_API
void      gst_child_proxy_set_property       (GstChildProxy * object, const gchar *name,
                                              const GValue *value);

GST_API
void      gst_child_proxy_set_valist         (GstChildProxy* object,
                                              const gchar * first_property_name,
                                              va_list var_args);
GST_API
void      gst_child_proxy_set                (GstChildProxy * object,
                                              const gchar * first_property_name,
                                              ...) G_GNUC_NULL_TERMINATED;
GST_API
void      gst_child_proxy_child_added        (GstChildProxy * parent, GObject * child,
                                              const gchar *name);
GST_API
void      gst_child_proxy_child_removed      (GstChildProxy * parent, GObject * child,
                                              const gchar *name);

G_END_DECLS

#endif /* __GST_CHILD_PROXY_H__ */
