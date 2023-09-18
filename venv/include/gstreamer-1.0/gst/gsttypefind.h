/* GStreamer
 * Copyright (C) 2003 Benjamin Otte <in7y118@public.uni-hamburg.de>
 *
 * gsttypefind.h: typefinding subsystem
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


#ifndef __GST_TYPE_FIND_H__
#define __GST_TYPE_FIND_H__

#include <gst/gstcaps.h>
#include <gst/gstplugin.h>
#include <gst/gstpluginfeature.h>

G_BEGIN_DECLS
/**
 * GST_TYPE_FIND_REGISTER_DEFINE_CUSTOM:
 * @type_find: The type find name in lower case, with words separated by '_'.
 * Used to generate `gst_type_find_register_*(GstPlugin* plugin)`.
 * @register_func: pointer to a method with the format: `gboolean register_func (GstPlugin* plugin);`
 *
 * A convenience macro to define the entry point of a
 * type find `gst_type_find_register_*(GstPlugin* plugin)` which uses
 * register_func as the main registration method for the type find.
 * As an example, you may define the type find named "custom-typefind"
 * as following using `type_find_register_custom`:
 *
 * ```
 * GST_TYPE_FIND_REGISTER_DEFINE_CUSTOM (plugin, type_find_register_custom)
 * ```
 *
 * Since: 1.20
 */
#define GST_TYPE_FIND_REGISTER_DEFINE_CUSTOM(type_find, register_func) \
G_BEGIN_DECLS \
gboolean G_PASTE (gst_type_find_register_, type_find) (GstPlugin * plugin) \
{ \
  return register_func (plugin); \
} \
G_END_DECLS

/**
 * GST_TYPE_FIND_REGISTER_DEFINE:
 * @t_f: The type find name in lower case, with words separated by '_'.
 * Used to generate `gst_type_find_register_*(GstPlugin* plugin)`.
 * @t_f_n: The public name of the type find
 * @r: The #GstRank of the type find (higher rank means more importance when autoplugging, see #GstRank)
 * @func: The #GstTypeFindFunction to use
 * @extensions: (nullable): Optional comma-separated list of extensions
 *     that could belong to this type
 * @possible_caps: (nullable): Optionally the caps that could be returned when typefinding
 *                 succeeds
 * @data: Optional user data. This user data must be available until the plugin
 *        is unloaded.
 * @data_notify: a #GDestroyNotify that will be called on @data when the plugin
 *        is unloaded.
 *
 * A convenience macro to define the entry point of a
 * type find `gst_type_find_register_*(GstPlugin* plugin)`.
 *
 * Since: 1.20
 */
#define GST_TYPE_FIND_REGISTER_DEFINE(t_f, t_f_n, r, func, extensions, possible_caps, data, data_notify) \
G_BEGIN_DECLS \
gboolean G_PASTE (gst_type_find_register_, t_f) (GstPlugin * plugin) \
{ \
  return gst_type_find_register (plugin, t_f_n, r, func, extensions, possible_caps, data, data_notify); \
} \
G_END_DECLS

/**
 * GST_TYPE_FIND_REGISTER_DECLARE:
 * @t_f: The type find name in lower case, with words separated by '_'.
 *
 * This macro can be used to declare a new type find.
 * It has to be used in combination with #GST_TYPE_FIND_REGISTER_DEFINE macro
 * and must be placed outside any block to declare the type find registration
 * function.
 *
 * Since: 1.20
 */
#define GST_TYPE_FIND_REGISTER_DECLARE(t_f) \
G_BEGIN_DECLS \
gboolean G_PASTE(gst_type_find_register_, t_f) (GstPlugin * plugin); \
G_END_DECLS

/**
 * GST_TYPE_FIND_REGISTER:
 * @t_f: The type find name in lower case, with words separated by '_'.
 * @plugin: The #GstPlugin where to register the type find.

 *
 * This macro can be used to register a type find into a #GstPlugin.
 * This method will be usually called in the plugin init function
 * but can also be called with a NULL plugin.
 *
 * Since: 1.20
 */
#define GST_TYPE_FIND_REGISTER(t_f, plugin) G_PASTE(gst_type_find_register_, t_f) (plugin)


#define GST_TYPE_TYPE_FIND  (gst_type_find_get_type())

typedef struct _GstTypeFind GstTypeFind;

/**
 * GstTypeFindFunction:
 * @find: A #GstTypeFind structure
 * @user_data: optional data to pass to the function
 *
 * A function that will be called by typefinding.
 */
typedef void (* GstTypeFindFunction) (GstTypeFind *find, gpointer user_data);

/**
 * GstTypeFindProbability:
 * @GST_TYPE_FIND_NONE: type undetected.
 * @GST_TYPE_FIND_MINIMUM: unlikely typefind.
 * @GST_TYPE_FIND_POSSIBLE: possible type detected.
 * @GST_TYPE_FIND_LIKELY: likely a type was detected.
 * @GST_TYPE_FIND_NEARLY_CERTAIN: nearly certain that a type was detected.
 * @GST_TYPE_FIND_MAXIMUM: very certain a type was detected.
 *
 * The probability of the typefind function. Higher values have more certainty
 * in doing a reliable typefind.
 */
typedef enum {
  GST_TYPE_FIND_NONE = 0,
  GST_TYPE_FIND_MINIMUM = 1,
  GST_TYPE_FIND_POSSIBLE = 50,
  GST_TYPE_FIND_LIKELY = 80,
  GST_TYPE_FIND_NEARLY_CERTAIN = 99,
  GST_TYPE_FIND_MAXIMUM = 100
} GstTypeFindProbability;

/**
 * GstTypeFind:
 * @peek: Method to peek data.
 * @suggest: Method to suggest #GstCaps with a given probability.
 * @data: The data used by the caller of the typefinding function.
 * @get_length: Returns the length of current data.
 *
 * Object that stores typefind callbacks. To use with #GstTypeFindFactory.
 */
struct _GstTypeFind {
  /* private to the caller of the typefind function */
  const guint8 *  (* peek)       (gpointer         data,
                                  gint64           offset,
                                  guint            size);

  void            (* suggest)    (gpointer         data,
                                  guint            probability,
                                  GstCaps         *caps);

  gpointer         data;

  /* optional */
  guint64         (* get_length) (gpointer data);

  /* <private> */
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * gst_type_find_get_type: (attributes doc.skip=true)
 */
GST_API
GType     gst_type_find_get_type   (void);

/* typefind function interface */

GST_API
const guint8 *  gst_type_find_peek       (GstTypeFind   * find,
                                          gint64          offset,
                                          guint           size);
GST_API
void            gst_type_find_suggest    (GstTypeFind   * find,
                                          guint           probability,
                                          GstCaps       * caps);
GST_API
void            gst_type_find_suggest_empty_simple (GstTypeFind * find,
                                                    guint         probability,
                                                    const char  * media_type);
GST_API
void            gst_type_find_suggest_simple (GstTypeFind * find,
                                              guint         probability,
                                              const char  * media_type,
                                              const char  * fieldname, ...) G_GNUC_NULL_TERMINATED;
GST_API
guint64   gst_type_find_get_length (GstTypeFind   * find);

/* registration interface */

GST_API
gboolean  gst_type_find_register   (GstPlugin            * plugin,
                                    const gchar          * name,
                                    guint                  rank,
                                    GstTypeFindFunction    func,
                                    const gchar          * extensions,
                                    GstCaps              * possible_caps,
                                    gpointer               data,
                                    GDestroyNotify         data_notify);

G_END_DECLS

#endif /* __GST_TYPE_FIND_H__ */
