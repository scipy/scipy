/* GStreamer
 * Copyright (C) <2003> David A. Schleef <ds@schleef.org>
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

#ifndef __GST_VALUE_H__
#define __GST_VALUE_H__

#include <gst/gstconfig.h>
#include <gst/gstcaps.h>
#include <gst/gststructure.h>
#include <gst/gstcapsfeatures.h>

G_BEGIN_DECLS

/**
 * GST_MAKE_FOURCC:
 * @a: the first character
 * @b: the second character
 * @c: the third character
 * @d: the fourth character
 *
 * Transform four characters into a #guint32 fourcc value with host
 * endianness.
 *
 * |[
 * guint32 fourcc = GST_MAKE_FOURCC ('M', 'J', 'P', 'G');
 * ]|
 *
 */
#define GST_MAKE_FOURCC(a,b,c,d) \
  ( (guint32)(a) | ((guint32) (b)) << 8  | ((guint32) (c)) << 16 | ((guint32) (d)) << 24 )

/**
 * GST_STR_FOURCC:
 * @f: a string with at least four characters
 *
 * Transform an input string into a #guint32 fourcc value with host
 * endianness.
 * Caller is responsible for ensuring the input string consists of at least
 * four characters.
 *
 * |[
 * guint32 fourcc = GST_STR_FOURCC ("MJPG");
 * ]|
 *
 */
#define GST_STR_FOURCC(f)               ((guint32)(((f)[0])|((f)[1]<<8)|((f)[2]<<16)|((f)[3]<<24)))

/**
 * GST_FOURCC_FORMAT: (skip):
 *
 * Can be used together with #GST_FOURCC_ARGS to properly output a
 * #guint32 fourcc value in a printf\()-style text message.
 *
 * |[
 * printf ("fourcc: %" GST_FOURCC_FORMAT "\n", GST_FOURCC_ARGS (fcc));
 * ]|
 *
 */
#define GST_FOURCC_FORMAT "c%c%c%c"

/**
 * GST_FOURCC_ARGS: (skip):
 * @fourcc: a #guint32 fourcc value to output
 *
 * Can be used together with #GST_FOURCC_FORMAT to properly output a
 * #guint32 fourcc value in a printf\()-style text message.
 */

#define __GST_PRINT_CHAR(c) \
  g_ascii_isprint(c) ? (c) : '.'
#define GST_FOURCC_ARGS(fourcc)               \
  __GST_PRINT_CHAR((fourcc) & 0xff),          \
  __GST_PRINT_CHAR(((fourcc) >> 8) & 0xff),   \
  __GST_PRINT_CHAR(((fourcc) >> 16) & 0xff),  \
  __GST_PRINT_CHAR(((fourcc) >> 24) & 0xff)
/**
 * GST_VALUE_HOLDS_INT_RANGE:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstIntRange value.
 */
#define GST_VALUE_HOLDS_INT_RANGE(x)      ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_int_range_type)

/**
 * GST_VALUE_HOLDS_INT64_RANGE:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstInt64Range value.
 */
#define GST_VALUE_HOLDS_INT64_RANGE(x)    ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_int64_range_type)

/**
 * GST_VALUE_HOLDS_DOUBLE_RANGE:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstDoubleRange value.
 */
#define GST_VALUE_HOLDS_DOUBLE_RANGE(x)   ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_double_range_type)

/**
 * GST_VALUE_HOLDS_FRACTION_RANGE:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstFractionRange value.
 */
#define GST_VALUE_HOLDS_FRACTION_RANGE(x) ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_fraction_range_type)

/**
 * GST_VALUE_HOLDS_LIST:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstValueList value.
 */
#define GST_VALUE_HOLDS_LIST(x)         ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_value_list_type)

/**
 * GST_VALUE_HOLDS_ARRAY:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstValueArray value.
 */
#define GST_VALUE_HOLDS_ARRAY(x)        ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_value_array_type)

/**
 * GST_VALUE_HOLDS_CAPS:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstCaps value.
 */
#define GST_VALUE_HOLDS_CAPS(x)         ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_caps_type)

/**
 * GST_VALUE_HOLDS_STRUCTURE:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstStructure value.
 */
#define GST_VALUE_HOLDS_STRUCTURE(x)      ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_structure_type)

/**
 * GST_VALUE_HOLDS_CAPS_FEATURES:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstCapsFeatures value.
 */
#define GST_VALUE_HOLDS_CAPS_FEATURES(x)  ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_caps_features_type)

/**
 * GST_VALUE_HOLDS_BUFFER:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstBuffer value.
 */
#define GST_VALUE_HOLDS_BUFFER(x)       ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_buffer_type)

/**
 * GST_VALUE_HOLDS_SAMPLE:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstSample value.
 */
#define GST_VALUE_HOLDS_SAMPLE(x)       ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_sample_type)

/**
 * GST_VALUE_HOLDS_FRACTION:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstFraction value.
 */
#define GST_VALUE_HOLDS_FRACTION(x)     ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_fraction_type)

/**
 * GST_VALUE_HOLDS_DATE_TIME:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstDateTime value.
 */
#define GST_VALUE_HOLDS_DATE_TIME(x)    ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_date_time_type)

/**
 * GST_VALUE_HOLDS_BITMASK:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstBitmask value.
 */
#define GST_VALUE_HOLDS_BITMASK(x)      ((gpointer)(x) != NULL && G_VALUE_TYPE(x) == _gst_bitmask_type)

/**
 * GST_VALUE_HOLDS_FLAG_SET:
 * @x: the #GValue to check
 *
 * Checks if the given #GValue contains a #GstFlagSet value.
 *
 * Since: 1.6
 */
#define GST_VALUE_HOLDS_FLAG_SET(x)     (G_TYPE_FUNDAMENTAL (G_VALUE_TYPE ((x))) == GST_TYPE_FLAG_SET)

/**
 * GST_FLAG_SET_MASK_EXACT: (value 4294967295) (type guint)
 *
 * A mask value with all bits set, for use as a
 * GstFlagSet mask where all flag bits must match
 * exactly
 *
 * Since: 1.6
 */
#define GST_FLAG_SET_MASK_EXACT ((guint)(-1))

GST_API GType _gst_int_range_type;

/**
 * GstIntRange:
 *
 * A fundamental type that describes a #gint range
 */

/**
 * GST_TYPE_INT_RANGE:
 *
 * a #GValue type that represents an integer range
 *
 * Returns: the #GType of GstIntRange
 */
#define GST_TYPE_INT_RANGE               (_gst_int_range_type)

GST_API GType _gst_int64_range_type;

/**
 * GstInt64Range:
 *
 * A fundamental type that describes a #gint64 range
 */

/**
 * GST_TYPE_INT64_RANGE:
 *
 * a #GValue type that represents an #gint64 range
 *
 * Returns: the #GType of GstInt64Range
 */
#define GST_TYPE_INT64_RANGE             (_gst_int64_range_type)

GST_API GType _gst_double_range_type;

/**
 * GstDoubleRange:
 *
 * A fundamental type that describes a #gdouble range
 */

/**
 * GST_TYPE_DOUBLE_RANGE:
 *
 * a #GValue type that represents a floating point range with double precision
 *
 * Returns: the #GType of GstIntRange
 */
#define GST_TYPE_DOUBLE_RANGE            (_gst_double_range_type)

GST_API GType _gst_fraction_range_type;

/**
 * GstFractionRange:
 *
 * A fundamental type that describes a #GstFractionRange range
 */

/**
 * GST_TYPE_FRACTION_RANGE:
 *
 * a #GValue type that represents a GstFraction range
 *
 * Returns: the #GType of GstFractionRange
 */
#define GST_TYPE_FRACTION_RANGE           (_gst_fraction_range_type)

GST_API GType _gst_value_list_type;

/**
 * GstValueList:
 *
 * A fundamental type that describes an unordered list of #GValue
 */

/**
 * GST_TYPE_LIST:
 *
 * a #GValue type that represents an unordered list of #GValue values. This
 * is used for example to express a list of possible values for a field in
 * a caps structure, like a list of possible sample rates, of which only one
 * will be chosen in the end. This means that all values in the list are
 * meaningful on their own.
 *
 * Returns: the #GType of GstValueList (which is not explicitly typed)
 */
#define GST_TYPE_LIST                    (_gst_value_list_type)

GST_API GType _gst_value_array_type;

/**
 * GstValueArray:
 *
 * A fundamental type that describes an ordered list of #GValue
 */

/**
 * GST_TYPE_ARRAY:
 *
 * a #GValue type that represents an ordered list of #GValue values. This is
 * used to express a set of values that is meaningful only in their specific
 * combination and order of values. Each value on its own is not particularly
 * meaningful, only the ordered array in its entirety is meaningful. This is
 * used for example to express channel layouts for multichannel audio where
 * each channel needs to be mapped to a position in the room.
 *
 * Returns: the #GType of GstArrayList (which is not explicitly typed)
 */
#define GST_TYPE_ARRAY                   (_gst_value_array_type)

GST_API GType _gst_fraction_type;

/**
 * GstFraction:
 *
 * A fundamental type that describes a fraction of an integer numerator
 * over an integer denominator
 */

/**
 * GST_TYPE_FRACTION:
 *
 * a #GValue type that represents a fraction of an integer numerator over
 * an integer denominator
 *
 * Returns: the #GType of GstFraction (which is not explicitly typed)
 */

#define GST_TYPE_FRACTION                (_gst_fraction_type)

GST_API GType _gst_bitmask_type;

/**
 * GstBitmask:
 *
 * A fundamental type that describes a 64-bit bitmask
 */

/**
 * GST_TYPE_BITMASK:
 *
 * a #GValue type that represents a 64-bit bitmask.
 *
 * Returns: the #GType of GstBitmask (which is not explicitly typed)
 */

#define GST_TYPE_BITMASK                 (_gst_bitmask_type)

GST_API GType _gst_flagset_type;

/**
 * GstFlagSet:
 *
 * A fundamental type that describes a 32-bit flag bitfield, with 32-bit
 * mask indicating which of the bits in the field are explicitly set.
 */

/**
 * GST_TYPE_FLAG_SET:
 *
 * a #GValue type that represents a 32-bit flag bitfield, with 32-bit
 * mask indicating which of the bits in the field are explicitly set.
 * Useful for negotiation.
 *
 * Returns: the #GType of GstFlags (which is not explicitly typed)
 *
 * Since: 1.6
 */
#define GST_TYPE_FLAG_SET                   (_gst_flagset_type)

/**
 * GST_TYPE_G_THREAD:
 *
 * a boxed #GValue type for #GThread that represents a thread.
 *
 * Returns: the #GType of GstGThread
 */

#define GST_TYPE_G_THREAD                gst_g_thread_get_type ()

/**
 * GST_VALUE_LESS_THAN:
 *
 * Indicates that the first value provided to a comparison function
 * (gst_value_compare()) is lesser than the second one.
 */
#define GST_VALUE_LESS_THAN              (-1)

/**
 * GST_VALUE_EQUAL:
 *
 * Indicates that the first value provided to a comparison function
 * (gst_value_compare()) is equal to the second one.
 */
#define GST_VALUE_EQUAL                   0

/**
 * GST_VALUE_GREATER_THAN:
 *
 * Indicates that the first value provided to a comparison function
 * (gst_value_compare()) is greater than the second one.
 */
#define GST_VALUE_GREATER_THAN            1

/**
 * GST_VALUE_UNORDERED:
 *
 * Indicates that the comparison function (gst_value_compare()) can not
 * determine a order for the two provided values.
 */
#define GST_VALUE_UNORDERED               2

/**
 * GstValueCompareFunc:
 * @value1: first value for comparison
 * @value2: second value for comparison
 *
 * Used together with gst_value_compare() to compare #GValue items.
 *
 * Returns: one of GST_VALUE_LESS_THAN, GST_VALUE_EQUAL, GST_VALUE_GREATER_THAN
 * or GST_VALUE_UNORDERED
 */
typedef gint     (* GstValueCompareFunc)     (const GValue *value1,
                                              const GValue *value2);

/**
 * GstValueSerializeFunc:
 * @value1: a #GValue
 *
 * Used by gst_value_serialize() to obtain a non-binary form of the #GValue.
 *
 * Free-function: g_free
 *
 * Returns: (transfer full): the string representation of the value
 */
typedef gchar *  (* GstValueSerializeFunc)   (const GValue *value1);

/**
 * GstValueDeserializeFunc:
 * @dest: a #GValue
 * @s: a string
 *
 * Used by gst_value_deserialize() to parse a non-binary form into the #GValue.
 *
 * Returns: %TRUE for success
 */
typedef gboolean (* GstValueDeserializeFunc) (GValue       *dest,
                                              const gchar  *s);

/**
 * GstValueDeserializeWithPSpecFunc:
 * @dest: a #GValue
 * @s: a string
 * @pspec: a #GParamSpec describing the expected value
 *
 * Used by gst_value_deserialize_with_pspec() to parse a non-binary form into the #GValue.
 *
 * Returns: %TRUE for success
 * Since: 1.20
 */
typedef gboolean (* GstValueDeserializeWithPSpecFunc) (GValue       *dest,
                                                       const gchar  *s,
                                                       GParamSpec   *pspec);


typedef struct _GstValueTable GstValueTable;
/**
 * GstValueTable:
 * @type: a #GType
 * @compare: a #GstValueCompareFunc
 * @serialize: a #GstValueSerializeFunc
 * @deserialize: a #GstValueDeserializeFunc
 * @deserialize_with_pspec: a #GstValueDeserializeWithPSpecFunc
 *
 * VTable for the #GValue @type.
 */
struct _GstValueTable {
  GType type;
  GstValueCompareFunc compare;
  GstValueSerializeFunc serialize;
  GstValueDeserializeFunc deserialize;

  /**
   * GstValueTable.deserialize_with_pspec:
   *
   * a #GstValueDeserializeWithPSpecFunc
   *
   * Since: 1.20
   */
  GstValueDeserializeWithPSpecFunc deserialize_with_pspec;

  /*< private >*/
  gpointer _gst_reserved [GST_PADDING - 1];
};

GST_API
GType gst_int_range_get_type (void);

GST_API
GType gst_int64_range_get_type (void);

GST_API
GType gst_double_range_get_type (void);

GST_API
GType gst_fraction_range_get_type (void);

GST_API
GType gst_fraction_get_type (void);

GST_API
GType gst_value_list_get_type (void);

GST_API
GType gst_value_array_get_type (void);

GST_API
GType gst_bitmask_get_type (void);

GST_API
GType gst_flagset_get_type (void);

/* Hide this compatibility type from introspection */
#ifndef __GI_SCANNER__
GST_API
GType gst_g_thread_get_type (void);
#endif

GST_API
void            gst_value_register              (const GstValueTable   *table);

GST_API
void            gst_value_init_and_copy         (GValue                *dest,
                                                 const GValue          *src);
GST_API
gchar *         gst_value_serialize             (const GValue          *value) G_GNUC_MALLOC;

GST_API
gboolean        gst_value_deserialize           (GValue                *dest,
                                                 const gchar           *src);

GST_API
gboolean        gst_value_deserialize_with_pspec (GValue               *dest,
                                                 const gchar           *src,
                                                 GParamSpec            *pspec);

/* list */

GST_API
void            gst_value_list_append_value     (GValue         *value,
                                                 const GValue   *append_value);
GST_API
void            gst_value_list_append_and_take_value (GValue         *value,
                                                 GValue   *append_value);
GST_API
void            gst_value_list_prepend_value    (GValue         *value,
                                                 const GValue   *prepend_value);
GST_API
void            gst_value_list_concat           (GValue         *dest,
                                                 const GValue   *value1,
                                                 const GValue   *value2);
GST_API
void            gst_value_list_merge            (GValue         *dest,
                                                 const GValue   *value1,
                                                 const GValue   *value2);
GST_API
guint           gst_value_list_get_size         (const GValue   *value);

GST_API
const GValue *  gst_value_list_get_value        (const GValue   *value,
                                                 guint          index);

GST_API
GValue *        gst_value_list_init             (GValue *value,
						 guint prealloc);
/* array */

GST_API
void            gst_value_array_append_value    (GValue         *value,
                                                 const GValue   *append_value);
GST_API
void            gst_value_array_append_and_take_value    (GValue         *value,
                                                 GValue   *append_value);
GST_API
void            gst_value_array_prepend_value   (GValue         *value,
                                                 const GValue   *prepend_value);
GST_API
guint           gst_value_array_get_size        (const GValue   *value);

GST_API
const GValue *  gst_value_array_get_value       (const GValue   *value,
                                                 guint          index);
GST_API
GValue *        gst_value_array_init            (GValue *value,
						 guint prealloc);

/* int range */

GST_API
void            gst_value_set_int_range         (GValue         *value,
                                                 gint           start,
                                                 gint           end);
GST_API
void            gst_value_set_int_range_step    (GValue         *value,
                                                 gint           start,
                                                 gint           end,
                                                 gint           step);
GST_API
gint            gst_value_get_int_range_min     (const GValue   *value);

GST_API
gint            gst_value_get_int_range_max     (const GValue   *value);

GST_API
gint            gst_value_get_int_range_step    (const GValue   *value);

/* int64 range */

GST_API
void            gst_value_set_int64_range       (GValue         *value,
                                                 gint64         start,
                                                 gint64         end);
GST_API
void            gst_value_set_int64_range_step  (GValue         *value,
                                                 gint64         start,
                                                 gint64         end,
                                                 gint64         step);
GST_API
gint64          gst_value_get_int64_range_min   (const GValue   *value);

GST_API
gint64          gst_value_get_int64_range_max   (const GValue   *value);

GST_API
gint64          gst_value_get_int64_range_step  (const GValue   *value);

/* double range */

GST_API
void            gst_value_set_double_range      (GValue         *value,
                                                 gdouble        start,
                                                 gdouble        end);
GST_API
gdouble         gst_value_get_double_range_min  (const GValue   *value);

GST_API
gdouble         gst_value_get_double_range_max  (const GValue   *value);

/* caps */

GST_API
const GstCaps * gst_value_get_caps              (const GValue   *value);

GST_API
void            gst_value_set_caps              (GValue         *value,
                                                 const GstCaps  *caps);

/* structure */

GST_API
const GstStructure *
                gst_value_get_structure         (const GValue   *value);

GST_API
void            gst_value_set_structure         (GValue         *value,
                                                 const GstStructure  *structure);

/* caps features */

GST_API
const GstCapsFeatures *
                gst_value_get_caps_features     (const GValue   *value);

GST_API
void            gst_value_set_caps_features     (GValue         *value,
                                                 const GstCapsFeatures  *features);

/* fraction */

GST_API
void            gst_value_set_fraction          (GValue         *value,
                                                 gint           numerator,
                                                 gint           denominator);
GST_API
gint            gst_value_get_fraction_numerator   (const GValue  *value);

GST_API
gint            gst_value_get_fraction_denominator (const GValue *value);

GST_API
gboolean        gst_value_fraction_multiply        (GValue         *product,
                                                    const GValue   *factor1,
                                                    const GValue   *factor2);
GST_API
gboolean        gst_value_fraction_subtract     (GValue * dest,
                                                 const GValue * minuend,
                                                 const GValue * subtrahend);

/* fraction range */

GST_API
void            gst_value_set_fraction_range    (GValue         *value,
                                                 const GValue   *start,
                                                 const GValue   *end);
GST_API
void            gst_value_set_fraction_range_full (GValue       *value,
                                                 gint numerator_start,
                                                 gint denominator_start,
                                                 gint numerator_end,
                                                 gint denominator_end);
GST_API
const GValue    *gst_value_get_fraction_range_min (const GValue *value);

GST_API
const GValue    *gst_value_get_fraction_range_max (const GValue *value);

/* bitmask */

GST_API
guint64         gst_value_get_bitmask           (const GValue   *value);

GST_API
void            gst_value_set_bitmask           (GValue         *value,
                                                 guint64         bitmask);
/* flagset */

GST_API
void            gst_value_set_flagset (GValue * value, guint flags, guint mask);

GST_API
guint           gst_value_get_flagset_flags (const GValue * value);

GST_API
guint           gst_value_get_flagset_mask (const GValue * value);

/* compare */

GST_API
gint            gst_value_compare               (const GValue   *value1,
                                                 const GValue   *value2);
GST_API
gboolean        gst_value_can_compare           (const GValue   *value1,
                                                 const GValue   *value2);
GST_API
gboolean        gst_value_is_subset             (const GValue   *value1,
                                                 const GValue   *value2);

/* union */

GST_API
gboolean        gst_value_union                 (GValue         *dest,
                                                 const GValue   *value1,
                                                 const GValue   *value2);
GST_API
gboolean        gst_value_can_union             (const GValue   *value1,
                                                 const GValue   *value2);

/* intersection */

GST_API
gboolean        gst_value_intersect             (GValue         *dest,
                                                 const GValue   *value1,
                                                 const GValue   *value2);
GST_API
gboolean        gst_value_can_intersect         (const GValue   *value1,
                                                 const GValue   *value2);

/* subtraction */

GST_API
gboolean        gst_value_subtract              (GValue         *dest,
                                                 const GValue   *minuend,
                                                 const GValue   *subtrahend);
GST_API
gboolean        gst_value_can_subtract          (const GValue   *minuend,
                                                 const GValue   *subtrahend);

/* fixation */

GST_API
gboolean        gst_value_is_fixed              (const GValue   *value);

GST_API
gboolean        gst_value_fixate                (GValue         *dest,
                                                 const GValue   *src);

/* Flagset registration wrapper */

GST_API
GType		gst_flagset_register (GType flags_type);

G_END_DECLS

#endif


