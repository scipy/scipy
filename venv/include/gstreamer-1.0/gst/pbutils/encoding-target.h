/* GStreamer encoding profile registry
 * Copyright (C) 2010 Edward Hervey <edward.hervey@collabora.co.uk>
 *           (C) 2010 Nokia Corporation
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

#ifndef __GST_PROFILE_REGISTRY_H__
#define __GST_PROFILE_REGISTRY_H__

#include <gst/pbutils/encoding-profile.h>

G_BEGIN_DECLS


/* FIXME/UNKNOWNS
 *
 * Should encoding categories be well-known strings/quarks ?
 *
 */

/**
 * GST_ENCODING_CATEGORY_DEVICE:
 *
 * #GstEncodingTarget category for device-specific targets.
 * The name of the target will usually be the constructor and model of the device,
 * and that target will contain #GstEncodingProfiles suitable for that device.
 */
#define GST_ENCODING_CATEGORY_DEVICE            "device"

/**
 * GST_ENCODING_CATEGORY_ONLINE_SERVICE:
 *
 * #GstEncodingTarget category for online-services.
 * The name of the target will usually be the name of the online service
 * and that target will contain #GstEncodingProfiles suitable for that online
 * service.
 */

#define GST_ENCODING_CATEGORY_ONLINE_SERVICE    "online-service"

/**
 * GST_ENCODING_CATEGORY_STORAGE_EDITING:
 *
 * #GstEncodingTarget category for storage, archiving and editing targets.
 * Those targets can be lossless and/or provide very fast random access content.
 * The name of the target will usually be the container type or editing target,
 * and that target will contain #GstEncodingProfiles suitable for editing or
 * storage.
 */
#define GST_ENCODING_CATEGORY_STORAGE_EDITING   "storage-editing"

/**
 * GST_ENCODING_CATEGORY_CAPTURE:
 *
 * #GstEncodingTarget category for recording and capture.
 * Targets within this category are optimized for low latency encoding.
 */
#define GST_ENCODING_CATEGORY_CAPTURE           "capture"

/**
 * GST_ENCODING_CATEGORY_FILE_EXTENSION:
 *
 * #GstEncodingTarget category for file extensions.
 * The name of the target will be the name of the file extensions possible
 * for a particular target. Those targets are defining like 'default' formats
 * usually used for a particular file extension.
 */

#define GST_ENCODING_CATEGORY_FILE_EXTENSION    "file-extension"

/**
 * GstEncodingTarget:
 *
 * Collection of #GstEncodingProfile for a specific target or use-case.
 *
 * When being stored/loaded, targets come from a specific category, like
 * #GST_ENCODING_CATEGORY_DEVICE.
 */
#define GST_TYPE_ENCODING_TARGET                        \
  (gst_encoding_target_get_type ())
#define GST_ENCODING_TARGET(obj)                        \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_ENCODING_TARGET, GstEncodingTarget))
#define GST_IS_ENCODING_TARGET(obj)                     \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_ENCODING_TARGET))

typedef struct _GstEncodingTarget GstEncodingTarget;
typedef GObjectClass GstEncodingTargetClass;

GST_PBUTILS_API
GType gst_encoding_target_get_type (void);

/**
 * gst_encoding_target_unref:
 * @target: a #GstEncodingTarget
 *
 * Decreases the reference count of the @target, possibly freeing it.
 */
#define gst_encoding_target_unref(target) \
  (g_object_unref ((GObject*) target))

/**
 * gst_encoding_target_ref:
 * @target: a #GstEncodingTarget
 *
 * Increases the reference count of the @target.
 */
#define gst_encoding_target_ref(target) \
  (g_object_ref ((GObject*) target))

GST_PBUTILS_API
GstEncodingTarget *     gst_encoding_target_new                 (const gchar *name,
                                                                 const gchar *category,
                                                                 const gchar *description,
                                                                 const GList *profiles);

GST_PBUTILS_API
const gchar *           gst_encoding_target_get_name            (GstEncodingTarget *target);

GST_PBUTILS_API
const gchar *           gst_encoding_target_get_category        (GstEncodingTarget *target);

GST_PBUTILS_API
const gchar *           gst_encoding_target_get_description     (GstEncodingTarget *target);

GST_PBUTILS_API
const gchar *           gst_encoding_target_get_path            (GstEncodingTarget *target);

GST_PBUTILS_API
const GList *           gst_encoding_target_get_profiles        (GstEncodingTarget *target);

GST_PBUTILS_API
GstEncodingProfile *    gst_encoding_target_get_profile         (GstEncodingTarget *target,
                                                                 const gchar *name);

GST_PBUTILS_API
gboolean                gst_encoding_target_add_profile         (GstEncodingTarget *target,
                                                                 GstEncodingProfile *profile);

GST_PBUTILS_API
gboolean                gst_encoding_target_save                (GstEncodingTarget *target,
                                                                 GError **error);

GST_PBUTILS_API
gboolean                gst_encoding_target_save_to_file        (GstEncodingTarget *target,
                                                                 const gchar *filepath,
                                                                 GError **error);

GST_PBUTILS_API
GstEncodingTarget *     gst_encoding_target_load                (const gchar *name,
                                                                 const gchar *category,
                                                                 GError **error);

GST_PBUTILS_API
GstEncodingTarget *     gst_encoding_target_load_from_file      (const gchar *filepath,
                                                                 GError **error);

GST_PBUTILS_API
GList *                 gst_encoding_list_available_categories  (void);

GST_PBUTILS_API
GList *                 gst_encoding_list_all_targets           (const gchar * categoryname);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstEncodingTarget, gst_object_unref)

G_END_DECLS

#endif  /* __GST_PROFILE_REGISTRY_H__ */
