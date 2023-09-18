/* GStreamer
 * Copyright (C) 2006 Stefan Kost <ensonic@users.sf.net>
 *
 * gstpreset.h: helper interface header for element presets
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

#ifndef __GST_PRESET_H__
#define __GST_PRESET_H__

#include <glib-object.h>
#include <gst/gstconfig.h>

G_BEGIN_DECLS

#define GST_TYPE_PRESET               (gst_preset_get_type())
#define GST_PRESET(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_PRESET, GstPreset))
#define GST_IS_PRESET(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_PRESET))
#define GST_PRESET_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), GST_TYPE_PRESET, GstPresetInterface))

/**
 * GstPreset:
 *
 * Opaque #GstPreset data structure.
 */
typedef struct _GstPreset GstPreset; /* dummy object */
typedef struct _GstPresetInterface GstPresetInterface;

/**
 * GstPresetInterface:
 * @parent: parent interface type.
 * @get_preset_names: virtual method to get list of presets
 * @get_property_names: virtual methods to get properties that are persistent
 * @load_preset: virtual methods to load a preset into properties
 * @save_preset: virtual methods to save properties into a preset
 * @rename_preset: virtual methods to rename a preset
 * @delete_preset: virtual methods to remove a preset
 * @set_meta: virtual methods to set textual meta data to a preset
 * @get_meta: virtual methods to get textual meta data from a preset
 *
 * #GstPreset interface.
 */
struct _GstPresetInterface
{
  GTypeInterface parent;

  /* methods */
  gchar**      (*get_preset_names)    (GstPreset *preset);

  gchar**      (*get_property_names)  (GstPreset *preset);

  gboolean     (*load_preset)         (GstPreset *preset, const gchar *name);
  gboolean     (*save_preset)         (GstPreset *preset, const gchar *name);
  gboolean     (*rename_preset)       (GstPreset *preset, const gchar *old_name,
                                       const gchar *new_name);
  gboolean     (*delete_preset)       (GstPreset *preset, const gchar *name);

  gboolean     (*set_meta)            (GstPreset *preset, const gchar *name,
                                       const gchar *tag, const gchar *value);
  gboolean     (*get_meta)            (GstPreset *preset, const gchar *name,
                                       const gchar *tag, gchar **value);
  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType        gst_preset_get_type (void);

GST_API
gchar**      gst_preset_get_preset_names   (GstPreset *preset) G_GNUC_MALLOC;

GST_API
gchar**      gst_preset_get_property_names (GstPreset *preset) G_GNUC_MALLOC;

GST_API
gboolean     gst_preset_load_preset        (GstPreset *preset, const gchar *name);

GST_API
gboolean     gst_preset_save_preset        (GstPreset *preset, const gchar *name);

GST_API
gboolean     gst_preset_rename_preset      (GstPreset *preset, const gchar *old_name,
                                            const gchar *new_name);
GST_API
gboolean     gst_preset_delete_preset      (GstPreset *preset, const gchar *name);

GST_API
gboolean     gst_preset_set_meta           (GstPreset *preset, const gchar *name,
                                            const gchar *tag, const gchar *value);
GST_API
gboolean     gst_preset_get_meta           (GstPreset *preset, const gchar *name,
                                            const gchar *tag, gchar **value);
GST_API
gboolean     gst_preset_set_app_dir        (const gchar *app_dir);

GST_API
const gchar *gst_preset_get_app_dir        (void);

GST_API
gboolean     gst_preset_is_editable        (GstPreset *preset);

G_END_DECLS

#endif /* __GST_PRESET_H__ */
