/* GStreamer
 * Copyright (C) 2004 Ronald Bultje <rbultje@ronald.bitfreak.net>
 *           (C) 2015 Wim Taymans <wim.taymans@gmail.com>
 *
 * audioconverter.h: audio format conversion library
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

#ifndef __GST_AUDIO_CONVERTER_H__
#define __GST_AUDIO_CONVERTER_H__

#include <gst/gst.h>
#include <gst/audio/audio.h>

G_BEGIN_DECLS

/**
 * GstAudioConverter:
 *
 * Opaque #GstAudioConverter struct.
 *
 * Since: 1.8
 */
typedef struct _GstAudioConverter GstAudioConverter;

/**
 * GST_AUDIO_CONVERTER_OPT_RESAMPLER_METHOD:
 *
 * #GstAudioResamplerMethod, The resampler method to use when
 * changing sample rates.
 * Default is #GST_AUDIO_RESAMPLER_METHOD_BLACKMAN_NUTTALL.
 */
#define GST_AUDIO_CONVERTER_OPT_RESAMPLER_METHOD   "GstAudioConverter.resampler-method"

/**
 * GST_AUDIO_CONVERTER_OPT_DITHER_METHOD:
 *
 * #GstAudioDitherMethod, The dither method to use when
 * changing bit depth.
 * Default is #GST_AUDIO_DITHER_NONE.
 */
#define GST_AUDIO_CONVERTER_OPT_DITHER_METHOD   "GstAudioConverter.dither-method"

/**
 * GST_AUDIO_CONVERTER_OPT_NOISE_SHAPING_METHOD:
 *
 * #GstAudioNoiseShapingMethod, The noise shaping method to use
 * to mask noise from quantization errors.
 * Default is #GST_AUDIO_NOISE_SHAPING_NONE.
 */
#define GST_AUDIO_CONVERTER_OPT_NOISE_SHAPING_METHOD   "GstAudioConverter.noise-shaping-method"

/**
 * GST_AUDIO_CONVERTER_OPT_QUANTIZATION:
 *
 * #G_TYPE_UINT, The quantization amount. Components will be
 * quantized to multiples of this value.
 * Default is 1
 */
#define GST_AUDIO_CONVERTER_OPT_QUANTIZATION   "GstAudioConverter.quantization"

/**
 * GST_AUDIO_CONVERTER_OPT_MIX_MATRIX:
 *
 * #GST_TYPE_LIST, The channel mapping matrix.
 *
 * The matrix coefficients must be between -1 and 1: the number of rows is equal
 * to the number of output channels and the number of columns is equal to the
 * number of input channels.
 *
 * ## Example matrix generation code
 * To generate the matrix using code:
 *
 * |[
 * GValue v = G_VALUE_INIT;
 * GValue v2 = G_VALUE_INIT;
 * GValue v3 = G_VALUE_INIT;
 *
 * g_value_init (&v2, GST_TYPE_ARRAY);
 * g_value_init (&v3, G_TYPE_DOUBLE);
 * g_value_set_double (&v3, 1);
 * gst_value_array_append_value (&v2, &v3);
 * g_value_unset (&v3);
 * [ Repeat for as many double as your input channels - unset and reinit v3 ]
 * g_value_init (&v, GST_TYPE_ARRAY);
 * gst_value_array_append_value (&v, &v2);
 * g_value_unset (&v2);
 * [ Repeat for as many v2's as your output channels - unset and reinit v2]
 * g_object_set_property (G_OBJECT (audiomixmatrix), "matrix", &v);
 * g_value_unset (&v);
 * ]|
 */
#define GST_AUDIO_CONVERTER_OPT_MIX_MATRIX   "GstAudioConverter.mix-matrix"

/**
 * GST_AUDIO_CONVERTER_OPT_DITHER_THRESHOLD:
 *
 * Threshold for the output bit depth at/below which to apply dithering.
 *
 * Default is 20 bit.
 *
 * Since: 1.22
 */
#define GST_AUDIO_CONVERTER_OPT_DITHER_THRESHOLD   "GstAudioConverter.dither-threshold"

/**
 * GstAudioConverterFlags:
 * @GST_AUDIO_CONVERTER_FLAG_NONE: no flag
 * @GST_AUDIO_CONVERTER_FLAG_IN_WRITABLE: the input sample arrays are writable and can be
 *    used as temporary storage during conversion.
 * @GST_AUDIO_CONVERTER_FLAG_VARIABLE_RATE: allow arbitrary rate updates with
 *    gst_audio_converter_update_config().
 *
 * Extra flags passed to gst_audio_converter_new() and gst_audio_converter_samples().
 */
typedef enum {
  GST_AUDIO_CONVERTER_FLAG_NONE            = 0,
  GST_AUDIO_CONVERTER_FLAG_IN_WRITABLE     = (1 << 0),
  GST_AUDIO_CONVERTER_FLAG_VARIABLE_RATE   = (1 << 1)
} GstAudioConverterFlags;

GST_AUDIO_API
GstAudioConverter *  gst_audio_converter_new             (GstAudioConverterFlags flags,
                                                          GstAudioInfo *in_info,
                                                          GstAudioInfo *out_info,
                                                          GstStructure *config);

GST_AUDIO_API
GType                gst_audio_converter_get_type        (void);

GST_AUDIO_API
void                 gst_audio_converter_free            (GstAudioConverter * convert);

GST_AUDIO_API
void                 gst_audio_converter_reset           (GstAudioConverter * convert);

GST_AUDIO_API
gboolean             gst_audio_converter_update_config   (GstAudioConverter * convert,
                                                          gint in_rate, gint out_rate,
                                                          GstStructure *config);

GST_AUDIO_API
const GstStructure * gst_audio_converter_get_config      (GstAudioConverter * convert,
                                                          gint *in_rate, gint *out_rate);

GST_AUDIO_API
gsize                gst_audio_converter_get_out_frames  (GstAudioConverter *convert,
                                                          gsize in_frames);

GST_AUDIO_API
gsize                gst_audio_converter_get_in_frames   (GstAudioConverter *convert,
                                                          gsize out_frames);

GST_AUDIO_API
gsize                gst_audio_converter_get_max_latency (GstAudioConverter *convert);

GST_AUDIO_API
gboolean             gst_audio_converter_samples         (GstAudioConverter * convert,
                                                          GstAudioConverterFlags flags,
                                                          gpointer in[], gsize in_frames,
                                                          gpointer out[], gsize out_frames);

GST_AUDIO_API
gboolean             gst_audio_converter_supports_inplace (GstAudioConverter *convert);

GST_AUDIO_API
gboolean             gst_audio_converter_is_passthrough   (GstAudioConverter *convert);

GST_AUDIO_API
gboolean             gst_audio_converter_convert          (GstAudioConverter * convert,
                                                           GstAudioConverterFlags flags,
                                                           gpointer in, gsize in_size,
                                                           gpointer *out, gsize *out_size);

G_END_DECLS

#endif /* __GST_AUDIO_CONVERTER_H__ */
