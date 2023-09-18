/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *
 * gstpluginfeature.h: Header for base GstPluginFeature
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


#ifndef __GST_PLUGIN_FEATURE_H__
#define __GST_PLUGIN_FEATURE_H__

#include <glib-object.h>
#include <gst/gstobject.h>
#include <gst/gstplugin.h>

G_BEGIN_DECLS

#define GST_TYPE_PLUGIN_FEATURE                 (gst_plugin_feature_get_type())
#define GST_PLUGIN_FEATURE(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_PLUGIN_FEATURE, GstPluginFeature))
#define GST_IS_PLUGIN_FEATURE(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_PLUGIN_FEATURE))
#define GST_PLUGIN_FEATURE_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_PLUGIN_FEATURE, GstPluginFeatureClass))
#define GST_IS_PLUGIN_FEATURE_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_PLUGIN_FEATURE))
#define GST_PLUGIN_FEATURE_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_PLUGIN_FEATURE, GstPluginFeatureClass))
#define GST_PLUGIN_FEATURE_CAST(obj)            ((GstPluginFeature*)(obj))

/**
 * GstPluginFeature:
 *
 * Opaque #GstPluginFeature structure.
 */
typedef struct _GstPluginFeature GstPluginFeature;
typedef struct _GstPluginFeatureClass GstPluginFeatureClass;

/**
 * GstRank:
 * @GST_RANK_NONE: will be chosen last or not at all
 * @GST_RANK_MARGINAL: unlikely to be chosen
 * @GST_RANK_SECONDARY: likely to be chosen
 * @GST_RANK_PRIMARY: will be chosen first
 *
 * Element priority ranks. Defines the order in which the autoplugger (or
 * similar rank-picking mechanisms, such as e.g. gst_element_make_from_uri())
 * will choose this element over an alternative one with the same function.
 *
 * These constants serve as a rough guidance for defining the rank of a
 * #GstPluginFeature. Any value is valid, including values bigger than
 * @GST_RANK_PRIMARY.
 */
typedef enum {
  GST_RANK_NONE                 = 0,
  GST_RANK_MARGINAL             = 64,
  GST_RANK_SECONDARY            = 128,
  GST_RANK_PRIMARY              = 256
} GstRank;

/**
 * gst_plugin_feature_get_name:
 * @feature: a #GstPluginFeature to get the name of @feature.
 *
 * Returns the name of @feature.
 * For a nameless plugin feature, this returns %NULL.
 *
 * Returns: (transfer none) (nullable): the name of @feature. MT safe.
 *
 */
#define                 gst_plugin_feature_get_name(feature)      GST_OBJECT_NAME(feature)

/**
 * gst_plugin_feature_set_name:
 * @feature: a #GstPluginFeature to set the name of.
 * @name: the new name
 *
 * Sets the name of the plugin feature, getting rid of the old name if there was one.
 */
#define                 gst_plugin_feature_set_name(feature,name) gst_object_set_name(GST_OBJECT_CAST(feature),name)

/**
 * GstPluginFeatureFilter:
 * @feature: the pluginfeature to check
 * @user_data: the user_data that has been passed on e.g.
 *  gst_registry_feature_filter()
 *
 * A function that can be used with e.g. gst_registry_feature_filter()
 * to get a list of pluginfeature that match certain criteria.
 *
 * Returns: %TRUE for a positive match, %FALSE otherwise
 */
typedef gboolean        (*GstPluginFeatureFilter)       (GstPluginFeature *feature,
                                                         gpointer user_data);

/* normal GObject stuff */

GST_API
GType           gst_plugin_feature_get_type             (void);

GST_API
GstPluginFeature *
                gst_plugin_feature_load                 (GstPluginFeature *feature);

GST_API
void            gst_plugin_feature_set_rank             (GstPluginFeature *feature, guint rank);

GST_API
guint           gst_plugin_feature_get_rank             (GstPluginFeature *feature);

GST_API
GstPlugin     * gst_plugin_feature_get_plugin           (GstPluginFeature *feature);

GST_API
const gchar   * gst_plugin_feature_get_plugin_name      (GstPluginFeature *feature);

GST_API
void            gst_plugin_feature_list_free            (GList *list);

GST_API
GList          *gst_plugin_feature_list_copy            (GList *list) G_GNUC_MALLOC;

GST_API
void            gst_plugin_feature_list_debug           (GList *list);

/**
 * GST_PLUGIN_FEATURE_LIST_DEBUG:
 * @list: (transfer none) (element-type Gst.PluginFeature): a #GList of
 *     plugin features
 *
 * Debug the plugin feature names in @list.
 */
#ifndef GST_DISABLE_GST_DEBUG
#define GST_PLUGIN_FEATURE_LIST_DEBUG(list) gst_plugin_feature_list_debug(list)
#else
#define GST_PLUGIN_FEATURE_LIST_DEBUG(list)
#endif

GST_API
gboolean        gst_plugin_feature_check_version        (GstPluginFeature *feature,
                                                         guint             min_major,
                                                         guint             min_minor,
                                                         guint             min_micro);
GST_API
gint            gst_plugin_feature_rank_compare_func    (gconstpointer p1,
							 gconstpointer p2);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstPluginFeature, gst_object_unref)

G_END_DECLS


#endif /* __GST_PLUGIN_FEATURE_H__ */

