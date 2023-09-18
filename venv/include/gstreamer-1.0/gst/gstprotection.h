/* GStreamer
 * Copyright (C) <2015> YouView TV Ltd.
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

#ifndef __GST_PROTECTION_H__
#define __GST_PROTECTION_H__

#include <gst/gst.h>

G_BEGIN_DECLS

/**
 * GST_PROTECTION_SYSTEM_ID_CAPS_FIELD:
 *
 * The field name in a GstCaps that is used to signal the UUID of the protection
 * system.
 *
 * Since: 1.6
 */
#define GST_PROTECTION_SYSTEM_ID_CAPS_FIELD "protection-system"

/**
 * GST_PROTECTION_UNSPECIFIED_SYSTEM_ID:
 *
 * The protection system value of the unspecified UUID.
 * In some cases the system protection ID is not present in the contents or in their
 * metadata, as encrypted WebM.
 * This define is used to set the value of the "system_id" field in GstProtectionEvent,
 * with this value, the application will use an external information to choose which
 * protection system to use.
 *
 * Example: The matroskademux uses this value in the case of encrypted WebM,
 * the application will choose the appropriate protection system based on the information
 * received through EME API.
 *
 * Since: 1.16
 */
#define GST_PROTECTION_UNSPECIFIED_SYSTEM_ID "unspecified-system-id"

typedef struct _GstProtectionMeta GstProtectionMeta;
/**
 * GstProtectionMeta:
 * @meta: the parent #GstMeta.
 * @info: the cryptographic information needed to decrypt the sample.
 *
 * Metadata type that holds information about a sample from a protection-protected
 * track, including the information needed to decrypt it (if it is encrypted).
 *
 * Since: 1.6
 */
struct _GstProtectionMeta
{
  GstMeta meta;

  GstStructure *info;
};

/**
 * gst_protection_meta_api_get_type: (attributes doc.skip=true)
 */
GST_API
GType gst_protection_meta_api_get_type (void);

#define GST_PROTECTION_META_API_TYPE (gst_protection_meta_api_get_type())

#define gst_buffer_get_protection_meta(b) \
    ((GstProtectionMeta*)gst_buffer_get_meta ((b), GST_PROTECTION_META_API_TYPE))

#define GST_PROTECTION_META_INFO (gst_protection_meta_get_info())

GST_API
const GstMetaInfo * gst_protection_meta_get_info (void);

GST_API
GstProtectionMeta * gst_buffer_add_protection_meta (GstBuffer * buffer,
                                                    GstStructure * info);
GST_API
const gchar * gst_protection_select_system (const gchar ** system_identifiers);

GST_API
gchar ** gst_protection_filter_systems_by_available_decryptors (
    const gchar ** system_identifiers);

G_END_DECLS
#endif /* __GST_PROTECTION_META_H__ */
