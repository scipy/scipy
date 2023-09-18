/* GStreamer
 * Copyright (C) 2003 David A. Schleef <ds@schleef.org>
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

#ifndef __GST_STRUCTURE_H__
#define __GST_STRUCTURE_H__

#include <gst/gstconfig.h>
#include <glib-object.h>
#include <gst/gstclock.h>
#include <gst/gstdatetime.h>
#include <gst/glib-compat.h>

G_BEGIN_DECLS

GST_API GType _gst_structure_type;

typedef struct _GstStructure GstStructure;

/**
 * GstSerializeFlags:
 * @GST_SERIALIZE_FLAG_NONE: No special flags specified.
 * @GST_SERIALIZE_FLAG_BACKWARD_COMPAT: Serialize using the old format for
 *                                      nested structures.
 *
 * Since: 1.20
 */
typedef enum
{
  GST_SERIALIZE_FLAG_NONE = 0,
  GST_SERIALIZE_FLAG_BACKWARD_COMPAT = (1 << 0),
} GstSerializeFlags;

#define GST_TYPE_STRUCTURE             (_gst_structure_type)
#define GST_IS_STRUCTURE(object)       ((object) && (GST_STRUCTURE(object)->type == GST_TYPE_STRUCTURE))
#define GST_STRUCTURE_CAST(object)     ((GstStructure *)(object))
#define GST_STRUCTURE(object)          (GST_STRUCTURE_CAST(object))


/**
 * GstStructureForeachFunc:
 * @field_id: the #GQuark of the field name
 * @value: the #GValue of the field
 * @user_data: user data
 *
 * A function that will be called in gst_structure_foreach(). The function may
 * not modify @value.
 *
 * Returns: %TRUE if the foreach operation should continue, %FALSE if
 * the foreach operation should stop with %FALSE.
 */
typedef gboolean (*GstStructureForeachFunc) (GQuark   field_id,
                                             const GValue * value,
                                             gpointer user_data);

/**
 * GstStructureMapFunc:
 * @field_id: the #GQuark of the field name
 * @value: the #GValue of the field
 * @user_data: user data
 *
 * A function that will be called in gst_structure_map_in_place(). The function
 * may modify @value.
 *
 * Returns: %TRUE if the map operation should continue, %FALSE if
 * the map operation should stop with %FALSE.
 */
typedef gboolean (*GstStructureMapFunc)     (GQuark   field_id,
                                             GValue * value,
                                             gpointer user_data);

/**
 * GstStructureFilterMapFunc:
 * @field_id: the #GQuark of the field name
 * @value: the #GValue of the field
 * @user_data: user data
 *
 * A function that will be called in gst_structure_filter_and_map_in_place().
 * The function may modify @value, and the value will be removed from
 * the structure if %FALSE is returned.
 *
 * Returns: %TRUE if the field should be preserved, %FALSE if it
 * should be removed.
 */
typedef gboolean (*GstStructureFilterMapFunc) (GQuark   field_id,
                                               GValue * value,
                                               gpointer user_data);

/**
 * GstStructure:
 * @type: the GType of a structure
 *
 * The GstStructure object. Most fields are private.
 */
struct _GstStructure {
  GType type;

  /*< private >*/
  GQuark name;
};

GST_API
GType                 gst_structure_get_type             (void);

GST_API
GstStructure *        gst_structure_new_empty            (const gchar * name) G_GNUC_MALLOC;

GST_API
GstStructure *        gst_structure_new_id_empty         (GQuark quark) G_GNUC_MALLOC;

GST_API
GstStructure *        gst_structure_new                  (const gchar * name,
                                                          const gchar * firstfield,
                                                          ...) G_GNUC_NULL_TERMINATED  G_GNUC_MALLOC;
GST_API
GstStructure *        gst_structure_new_valist           (const gchar * name,
                                                          const gchar * firstfield,
                                                          va_list       varargs) G_GNUC_MALLOC;
GST_API
GstStructure *        gst_structure_new_id               (GQuark name_quark,
                                                          GQuark field_quark,
                                                          ...) G_GNUC_MALLOC;
GST_API
GstStructure *        gst_structure_new_from_string      (const gchar * string);

GST_API
GstStructure *        gst_structure_copy                 (const GstStructure  * structure) G_GNUC_MALLOC;

GST_API
gboolean              gst_structure_set_parent_refcount  (GstStructure        * structure,
                                                          gint                * refcount);
GST_API
void                  gst_structure_free                 (GstStructure        * structure);

GST_API
void                  gst_clear_structure                (GstStructure **structure_ptr);
#define               gst_clear_structure(structure_ptr) g_clear_pointer ((structure_ptr), gst_structure_free)

GST_API
gboolean              gst_structure_take                 (GstStructure ** oldstr_ptr,
                                                          GstStructure * newstr);
GST_API
const gchar *         gst_structure_get_name             (const GstStructure  * structure);

GST_API
GQuark                gst_structure_get_name_id          (const GstStructure  * structure);

GST_API
gboolean              gst_structure_has_name             (const GstStructure  * structure,
                                                          const gchar         * name);
GST_API
void                  gst_structure_set_name             (GstStructure        * structure,
                                                          const gchar         * name);
GST_API
void                  gst_structure_id_set_value         (GstStructure        * structure,
                                                          GQuark                field,
                                                          const GValue        * value);
GST_API
void                  gst_structure_set_value            (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          const GValue        * value);
GST_API
void                  gst_structure_set_array            (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          const GValueArray   * array);
GST_API
void                  gst_structure_set_list             (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          const GValueArray   * array);
GST_API
void                  gst_structure_id_take_value        (GstStructure        * structure,
                                                          GQuark                field,
                                                          GValue              * value);
GST_API
void                  gst_structure_take_value           (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          GValue              * value);
GST_API
void                  gst_structure_set                  (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          ...) G_GNUC_NULL_TERMINATED;
GST_API
void                  gst_structure_set_valist           (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          va_list varargs);
GST_API
void                  gst_structure_id_set               (GstStructure        * structure,
                                                          GQuark                fieldname,
                                                          ...) G_GNUC_NULL_TERMINATED;
GST_API
void                  gst_structure_id_set_valist        (GstStructure        * structure,
                                                          GQuark                fieldname,
                                                          va_list varargs);
GST_API
gboolean              gst_structure_get_valist           (const GstStructure  * structure,
                                                          const char          * first_fieldname,
                                                          va_list              args);
GST_API
gboolean              gst_structure_get                  (const GstStructure  * structure,
                                                          const char          * first_fieldname,
                                                          ...) G_GNUC_NULL_TERMINATED;
GST_API
gboolean              gst_structure_id_get_valist        (const GstStructure  * structure,
                                                          GQuark                first_field_id,
                                                          va_list               args);
GST_API
gboolean              gst_structure_id_get               (const GstStructure  * structure,
                                                          GQuark                first_field_id,
                                                          ...) G_GNUC_NULL_TERMINATED;
GST_API
const GValue *        gst_structure_id_get_value         (const GstStructure  * structure,
                                                          GQuark                field);
GST_API
const GValue *        gst_structure_get_value            (const GstStructure  * structure,
                                                          const gchar         * fieldname);
GST_API
void                  gst_structure_remove_field         (GstStructure        * structure,
                                                          const gchar         * fieldname);
GST_API
void                  gst_structure_remove_fields        (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          ...) G_GNUC_NULL_TERMINATED;
GST_API
void                  gst_structure_remove_fields_valist (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          va_list               varargs);
GST_API
void                  gst_structure_remove_all_fields    (GstStructure        * structure);

GST_API
GType                 gst_structure_get_field_type       (const GstStructure  * structure,
                                                          const gchar         * fieldname);
GST_API
gboolean              gst_structure_foreach              (const GstStructure  * structure,
                                                          GstStructureForeachFunc   func,
                                                          gpointer              user_data);
GST_API
gboolean              gst_structure_map_in_place         (GstStructure        * structure,
                                                          GstStructureMapFunc   func,
                                                          gpointer              user_data);
GST_API
void                  gst_structure_filter_and_map_in_place (GstStructure        * structure,
                                                          GstStructureFilterMapFunc   func,
                                                          gpointer              user_data);
GST_API
gint                  gst_structure_n_fields             (const GstStructure  * structure);

GST_API
const gchar *         gst_structure_nth_field_name       (const GstStructure  * structure,
                                                          guint                 index);
GST_API
gboolean              gst_structure_id_has_field         (const GstStructure  * structure,
                                                          GQuark                field);
GST_API
gboolean              gst_structure_id_has_field_typed   (const GstStructure  * structure,
                                                          GQuark                field,
                                                          GType                 type);
GST_API
gboolean              gst_structure_has_field            (const GstStructure  * structure,
                                                          const gchar         * fieldname);
GST_API
gboolean              gst_structure_has_field_typed      (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          GType                 type);

/* utility functions */

GST_API
gboolean              gst_structure_get_boolean          (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          gboolean            * value);
GST_API
gboolean              gst_structure_get_int              (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          gint                * value);
GST_API
gboolean              gst_structure_get_uint             (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          guint               * value);
GST_API
gboolean              gst_structure_get_int64            (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          gint64              * value);
GST_API
gboolean              gst_structure_get_uint64           (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          guint64             * value);
GST_API
gboolean              gst_structure_get_double           (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          gdouble             * value);
GST_API
gboolean              gst_structure_get_date             (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          GDate              ** value);
GST_API
gboolean              gst_structure_get_date_time        (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          GstDateTime        ** value);
GST_API
gboolean              gst_structure_get_clock_time       (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          GstClockTime        * value);
GST_API
const gchar *         gst_structure_get_string           (const GstStructure  * structure,
                                                          const gchar         * fieldname);
GST_API
gboolean              gst_structure_get_enum             (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          GType                 enumtype,
                                                          gint                * value);
GST_API
gboolean              gst_structure_get_fraction         (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          gint                * value_numerator,
                                                          gint                * value_denominator);
GST_API
gboolean              gst_structure_get_flagset          (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          guint               * value_flags,
                                                          guint               * value_mask);
GST_API
gboolean              gst_structure_get_array            (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          GValueArray        ** array);
GST_API
gboolean              gst_structure_get_list             (GstStructure        * structure,
                                                          const gchar         * fieldname,
                                                          GValueArray        ** array);

GST_API
gboolean              gst_structure_get_flags            (const GstStructure  * structure,
                                                          const gchar         * fieldname,
                                                          GType                 flags_type,
                                                          guint               * value);

GST_API
gchar *               gst_structure_to_string            (const GstStructure * structure) G_GNUC_MALLOC;
GST_API
gchar *               gst_structure_serialize            (const GstStructure * structure,
                                                          GstSerializeFlags flags) G_GNUC_MALLOC;

GST_API
GstStructure *        gst_structure_from_string  (const gchar * string,
                                                  gchar      ** end) G_GNUC_MALLOC;
GST_API
gboolean              gst_structure_fixate_field_nearest_int      (GstStructure * structure,
                                                                   const char   * field_name,
                                                                   int            target);
GST_API
gboolean              gst_structure_fixate_field_nearest_double   (GstStructure * structure,
                                                                   const char   * field_name,
                                                                   double         target);
GST_API
gboolean              gst_structure_fixate_field_boolean          (GstStructure * structure,
                                                                   const char   * field_name,
                                                                   gboolean       target);
GST_API
gboolean              gst_structure_fixate_field_string           (GstStructure * structure,
                                                                   const char   * field_name,
                                                                   const gchar  * target);
GST_API
gboolean              gst_structure_fixate_field_nearest_fraction (GstStructure * structure,
                                                                   const char   * field_name,
                                                                   const gint     target_numerator,
                                                                   const gint     target_denominator);
GST_API
gboolean              gst_structure_fixate_field  (GstStructure * structure,
                                                   const char   * field_name);
GST_API
void                  gst_structure_fixate        (GstStructure * structure);

GST_API
gboolean              gst_structure_is_equal      (const GstStructure * structure1,
                                                   const GstStructure * structure2);
GST_API
gboolean              gst_structure_is_subset     (const GstStructure * subset,
                                                   const GstStructure * superset);
GST_API
gboolean              gst_structure_can_intersect (const GstStructure * struct1,
                                                   const GstStructure * struct2);
GST_API
GstStructure *        gst_structure_intersect     (const GstStructure * struct1,
                                                   const GstStructure * struct2) G_GNUC_MALLOC;

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstStructure, gst_structure_free)

G_END_DECLS

#endif

