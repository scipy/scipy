/* GStreamer
 *
 * Copyright (C) 2007 Sebastian Dröge <slomo@circular-chaos.org>
 *
 * gstlfocontrolsource.h: Control source that provides some periodic waveforms
 *                        as control values.
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

#ifndef __GST_LFO_CONTROL_SOURCE_H__
#define __GST_LFO_CONTROL_SOURCE_H__

#include <glib-object.h>
#include <gst/gst.h>
#include <gst/controller/controller-enumtypes.h>

G_BEGIN_DECLS

#define GST_TYPE_LFO_CONTROL_SOURCE \
  (gst_lfo_control_source_get_type ())
#define GST_LFO_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_LFO_CONTROL_SOURCE, GstLFOControlSource))
#define GST_LFO_CONTROL_SOURCE_CLASS(vtable) \
  (G_TYPE_CHECK_CLASS_CAST ((vtable), GST_TYPE_LFO_CONTROL_SOURCE, GstLFOControlSourceClass))
#define GST_IS_LFO_CONTROL_SOURCE(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_LFO_CONTROL_SOURCE))
#define GST_IS_LFO_CONTROL_SOURCE_CLASS(vtable) \
  (G_TYPE_CHECK_CLASS_TYPE ((vtable), GST_TYPE_LFO_CONTROL_SOURCE))
#define GST_LFO_CONTROL_SOURCE_GET_CLASS(inst) \
  (G_TYPE_INSTANCE_GET_CLASS ((inst), GST_TYPE_LFO_CONTROL_SOURCE, GstLFOControlSourceClass))

typedef struct _GstLFOControlSource GstLFOControlSource;
typedef struct _GstLFOControlSourceClass GstLFOControlSourceClass;
typedef struct _GstLFOControlSourcePrivate GstLFOControlSourcePrivate;

/**
 * GstLFOWaveform:
 * @GST_LFO_WAVEFORM_SINE: sine waveform
 * @GST_LFO_WAVEFORM_SQUARE: square waveform
 * @GST_LFO_WAVEFORM_SAW: saw waveform
 * @GST_LFO_WAVEFORM_REVERSE_SAW: reverse saw waveform
 * @GST_LFO_WAVEFORM_TRIANGLE: triangle waveform
 *
 * The various waveform modes available.
 */
typedef enum
{
  GST_LFO_WAVEFORM_SINE,
  GST_LFO_WAVEFORM_SQUARE,
  GST_LFO_WAVEFORM_SAW,
  GST_LFO_WAVEFORM_REVERSE_SAW,
  GST_LFO_WAVEFORM_TRIANGLE
} GstLFOWaveform;

/**
 * GstLFOControlSource:
 *
 * The instance structure of #GstControlSource.
 */
struct _GstLFOControlSource {
  GstControlSource parent;

  /* <private> */
  GstLFOControlSourcePrivate *priv;
  GMutex lock;
  gpointer _gst_reserved[GST_PADDING];
};

struct _GstLFOControlSourceClass {
  GstControlSourceClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_CONTROLLER_API
GType gst_lfo_control_source_get_type (void);

/* Functions */

GST_CONTROLLER_API
GstControlSource *gst_lfo_control_source_new (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstLFOControlSource, gst_object_unref)

G_END_DECLS

#endif /* __GST_LFO_CONTROL_SOURCE_H__ */
