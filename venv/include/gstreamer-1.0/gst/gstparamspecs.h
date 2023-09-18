/* GStreamer - GParamSpecs for some of our types
 * Copyright (C) 2007 Tim-Philipp Müller  <tim centricular net>
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

#ifndef __GST_PARAMSPECS_H__
#define __GST_PARAMSPECS_H__

#include <gst/gstvalue.h>

G_BEGIN_DECLS

/* --- paramspec flags */

/**
 * GST_PARAM_CONTROLLABLE: (value 512)
 *
 * Use this flag on GObject properties to signal they can make sense to be.
 * controlled over time. This hint is used by the GstController.
 */
#define	GST_PARAM_CONTROLLABLE	(1 << (G_PARAM_USER_SHIFT + 1))

/**
 * GST_PARAM_MUTABLE_READY: (value 1024)
 *
 * Use this flag on GObject properties of GstElements to indicate that
 * they can be changed when the element is in the READY or lower state.
 */
#define GST_PARAM_MUTABLE_READY  (1 << (G_PARAM_USER_SHIFT + 2))

/**
 * GST_PARAM_MUTABLE_PAUSED: (value 2048)
 *
 * Use this flag on GObject properties of GstElements to indicate that
 * they can be changed when the element is in the PAUSED or lower state.
 * This flag implies GST_PARAM_MUTABLE_READY.
 */
#define GST_PARAM_MUTABLE_PAUSED  (1 << (G_PARAM_USER_SHIFT + 3))

/**
 * GST_PARAM_MUTABLE_PLAYING: (value 4096)
 *
 * Use this flag on GObject properties of GstElements to indicate that
 * they can be changed when the element is in the PLAYING or lower state.
 * This flag implies GST_PARAM_MUTABLE_PAUSED.
 */
#define GST_PARAM_MUTABLE_PLAYING  (1 << (G_PARAM_USER_SHIFT + 4))

/**
 * GST_PARAM_DOC_SHOW_DEFAULT: (value 8192)
 *
 * Use this flag on GObject properties of GstObject to indicate that
 * during `gst-inspect` and friends, the default value should be used
 * as default instead of the current value.
 *
 * Since: 1.18
 */
#define GST_PARAM_DOC_SHOW_DEFAULT  (1 << (G_PARAM_USER_SHIFT + 5))

/**
 * GST_PARAM_CONDITIONALLY_AVAILABLE: (value 16384)
 *
 * Use this flag on GObject properties of GstObject to indicate that
 * they might not be available depending on environment such as OS, device, etc,
 * so such properties will be installed conditionally only if the GstObject is
 * able to support it.
 *
 * Since: 1.18
 */
#define GST_PARAM_CONDITIONALLY_AVAILABLE  (1 << (G_PARAM_USER_SHIFT + 6))

/**
 * GST_PARAM_USER_SHIFT: (value 65536)
 *
 * Bits based on GST_PARAM_USER_SHIFT can be used by 3rd party applications.
 */
#define	GST_PARAM_USER_SHIFT	(1 << (G_PARAM_USER_SHIFT + 8))


/* --- type macros --- */

/**
 * GstParamFraction:
 *
 * A fundamental type that describes a #GParamSpec for fractional
 * properties
 */

#define GST_TYPE_PARAM_FRACTION           (gst_param_spec_fraction_get_type ())
#define GST_IS_PARAM_SPEC_FRACTION(pspec) (G_TYPE_CHECK_INSTANCE_TYPE ((pspec), GST_TYPE_PARAM_FRACTION))
#define GST_PARAM_SPEC_FRACTION(pspec)    (G_TYPE_CHECK_INSTANCE_CAST ((pspec), GST_TYPE_PARAM_FRACTION, GstParamSpecFraction))

/**
 * GstParamArray:
 *
 * A fundamental type that describes a #GParamSpec for arrays of
 * values
 *
 * Since: 1.12
 */

#define GST_TYPE_PARAM_ARRAY_LIST           (gst_param_spec_array_get_type ())
#define GST_IS_PARAM_SPEC_ARRAY_LIST(pspec) (G_TYPE_CHECK_INSTANCE_TYPE ((pspec), GST_TYPE_PARAM_ARRAY_LIST))
#define GST_PARAM_SPEC_ARRAY_LIST(pspec)    (G_TYPE_CHECK_INSTANCE_CAST ((pspec), GST_TYPE_PARAM_ARRAY_LIST, GstParamSpecArray))


/* --- get_type functions --- */

GST_API
GType  gst_param_spec_fraction_get_type (void);

GST_API
GType  gst_param_spec_array_get_type (void);


/* --- typedefs & structures --- */

typedef struct _GstParamSpecFraction GstParamSpecFraction;
typedef struct _GstParamSpecArray GstParamSpecArray;

/**
 * GstParamSpecFraction:
 * @parent_instance: super class
 * @min_num: minimal numerator
 * @min_den: minimal denominator
 * @max_num: maximal numerator
 * @max_den: maximal denominator
 * @def_num: default numerator
 * @def_den: default denominator
 *
 * A GParamSpec derived structure that contains the meta data for fractional
 * properties.
 */
struct _GstParamSpecFraction {
  GParamSpec    parent_instance;

  gint          min_num, min_den;
  gint          max_num, max_den;
  gint          def_num, def_den;
};

/**
 * GstParamSpecArray:
 * @parent_instance: super class
 * @element_spec: the #GParamSpec of the type of values in the array
 *
 * A GParamSpec derived structure for arrays of values.
 */
struct _GstParamSpecArray {
  GParamSpec    parent_instance;

  GParamSpec * element_spec;
};


/* --- GParamSpec prototypes --- */

GST_API
GParamSpec  * gst_param_spec_fraction (const gchar * name,
                                       const gchar * nick,
                                       const gchar * blurb,
                                       gint min_num, gint min_denom,
                                       gint max_num, gint max_denom,
                                       gint default_num, gint default_denom,
                                       GParamFlags flags) G_GNUC_MALLOC;
GST_API
GParamSpec  * gst_param_spec_array    (const gchar * name,
                                       const gchar * nick,
                                       const gchar * blurb,
                                       GParamSpec * element_spec,
                                       GParamFlags flags) G_GNUC_MALLOC;

G_END_DECLS

#endif /* __GST_PARAMSPECS_H__ */

