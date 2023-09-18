/* GStreamer
 * Copyright (C) 2013 Collabora Ltd.
 *   Author: Sebastian Dröge <sebastian.droege@collabora.co.uk>
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

#ifndef __GST_CAPS_FEATURES_H__
#define __GST_CAPS_FEATURES_H__

#include <gst/gstconfig.h>
#include <gst/glib-compat.h>
#include <glib-object.h>
#include <glib.h>

G_BEGIN_DECLS

typedef struct _GstCapsFeatures GstCapsFeatures;

GST_API GType _gst_caps_features_type;

#define GST_TYPE_CAPS_FEATURES (_gst_caps_features_type)
#define GST_IS_CAPS_FEATURES(object)       (gst_is_caps_features(object))
#define GST_CAPS_FEATURES_CAST(object)     ((GstCapsFeatures *)(object))
#define GST_CAPS_FEATURES(object)          (GST_CAPS_FEATURES_CAST(object))

#define GST_CAPS_FEATURE_MEMORY_SYSTEM_MEMORY "memory:SystemMemory"

GST_API GstCapsFeatures *_gst_caps_features_any;
#define GST_CAPS_FEATURES_ANY (_gst_caps_features_any)

GST_API GstCapsFeatures *_gst_caps_features_memory_system_memory;
#define GST_CAPS_FEATURES_MEMORY_SYSTEM_MEMORY (_gst_caps_features_memory_system_memory)

GST_API
GType             gst_caps_features_get_type (void);

GST_API
gboolean          gst_is_caps_features (gconstpointer obj);

GST_API
GstCapsFeatures * gst_caps_features_new_empty (void);

GST_API
GstCapsFeatures * gst_caps_features_new_any (void);

GST_API
GstCapsFeatures * gst_caps_features_new_single (const gchar *feature) G_GNUC_MALLOC;

GST_API
GstCapsFeatures * gst_caps_features_new (const gchar *feature1, ...) G_GNUC_NULL_TERMINATED;

GST_API
GstCapsFeatures * gst_caps_features_new_valist (const gchar *feature1, va_list varargs);

GST_API
GstCapsFeatures * gst_caps_features_new_id (GQuark feature1, ...);

GST_API
GstCapsFeatures * gst_caps_features_new_id_valist (GQuark feature1, va_list varargs);

GST_API
gboolean          gst_caps_features_set_parent_refcount  (GstCapsFeatures *features, gint * refcount);

GST_API
GstCapsFeatures * gst_caps_features_copy (const GstCapsFeatures * features);

GST_API
void              gst_caps_features_free (GstCapsFeatures * features);

GST_API
gchar *           gst_caps_features_to_string (const GstCapsFeatures * features);

GST_API
GstCapsFeatures * gst_caps_features_from_string (const gchar * features);

GST_API
guint             gst_caps_features_get_size (const GstCapsFeatures * features);

GST_API
const gchar *     gst_caps_features_get_nth (const GstCapsFeatures * features, guint i);

GST_API
GQuark            gst_caps_features_get_nth_id (const GstCapsFeatures * features, guint i);

GST_API
gboolean          gst_caps_features_contains (const GstCapsFeatures * features, const gchar * feature);

GST_API
gboolean          gst_caps_features_contains_id (const GstCapsFeatures * features, GQuark feature);

GST_API
gboolean          gst_caps_features_is_equal (const GstCapsFeatures * features1, const GstCapsFeatures * features2);

GST_API
gboolean          gst_caps_features_is_any (const GstCapsFeatures * features);

GST_API
void              gst_caps_features_add (GstCapsFeatures * features, const gchar * feature);

GST_API
void              gst_caps_features_add_id ( GstCapsFeatures * features, GQuark feature);

GST_API
void              gst_caps_features_remove (GstCapsFeatures * features, const gchar * feature);

GST_API
void              gst_caps_features_remove_id (GstCapsFeatures * features, GQuark feature);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstCapsFeatures, gst_caps_features_free)

G_END_DECLS

#endif /* __GST_CAPS_FEATURES_H__ */
