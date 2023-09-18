/* GStreamer
 * Copyright (C) 2007 Sebastian Dröge <slomo@circular-chaos.org>
 *           (C) 2015 Wim Taymans <wim.taymans@gmail.com>
 *
 * gstaudioquantize.h: quantizes audio to the target format and optionally
 *                     applies dithering and noise shaping.
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

#include <gst/gst.h>

#include <gst/audio/audio.h>


#ifndef __GST_AUDIO_QUANTIZE_H__
#define __GST_AUDIO_QUANTIZE_H__

/**
 * GstAudioDitherMethod:
 * @GST_AUDIO_DITHER_NONE: No dithering
 * @GST_AUDIO_DITHER_RPDF: Rectangular dithering
 * @GST_AUDIO_DITHER_TPDF: Triangular dithering (default)
 * @GST_AUDIO_DITHER_TPDF_HF: High frequency triangular dithering
 *
 * Set of available dithering methods.
 */
typedef enum
{
  GST_AUDIO_DITHER_NONE = 0,
  GST_AUDIO_DITHER_RPDF,
  GST_AUDIO_DITHER_TPDF,
  GST_AUDIO_DITHER_TPDF_HF
} GstAudioDitherMethod;

/**
 * GstAudioNoiseShapingMethod:
 * @GST_AUDIO_NOISE_SHAPING_NONE: No noise shaping (default)
 * @GST_AUDIO_NOISE_SHAPING_ERROR_FEEDBACK: Error feedback
 * @GST_AUDIO_NOISE_SHAPING_SIMPLE: Simple 2-pole noise shaping
 * @GST_AUDIO_NOISE_SHAPING_MEDIUM: Medium 5-pole noise shaping
 * @GST_AUDIO_NOISE_SHAPING_HIGH: High 8-pole noise shaping
 *
 * Set of available noise shaping methods
 */
typedef enum
{
  GST_AUDIO_NOISE_SHAPING_NONE = 0,
  GST_AUDIO_NOISE_SHAPING_ERROR_FEEDBACK,
  GST_AUDIO_NOISE_SHAPING_SIMPLE,
  GST_AUDIO_NOISE_SHAPING_MEDIUM,
  GST_AUDIO_NOISE_SHAPING_HIGH
} GstAudioNoiseShapingMethod;

/**
 * GstAudioQuantizeFlags:
 * @GST_AUDIO_QUANTIZE_FLAG_NONE: no flags
 * @GST_AUDIO_QUANTIZE_FLAG_NON_INTERLEAVED: samples are non-interleaved
 *
 * Extra flags that can be passed to gst_audio_quantize_new()
 */
typedef enum
{
  GST_AUDIO_QUANTIZE_FLAG_NONE            = 0,
  GST_AUDIO_QUANTIZE_FLAG_NON_INTERLEAVED = (1 << 0)
} GstAudioQuantizeFlags;


typedef struct _GstAudioQuantize GstAudioQuantize;

GST_AUDIO_API
GstAudioQuantize *  gst_audio_quantize_new      (GstAudioDitherMethod dither,
                                                 GstAudioNoiseShapingMethod ns,
                                                 GstAudioQuantizeFlags flags,
                                                 GstAudioFormat format,
                                                 guint channels,
                                                 guint quantizer);

GST_AUDIO_API
void                gst_audio_quantize_free     (GstAudioQuantize * quant);

GST_AUDIO_API
void                gst_audio_quantize_reset    (GstAudioQuantize * quant);

GST_AUDIO_API
void                gst_audio_quantize_samples  (GstAudioQuantize * quant,
                                                 const gpointer in[],
                                                 gpointer out[], guint samples);

#endif /* __GST_AUDIO_QUANTIZE_H__ */
