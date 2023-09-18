/* GStreamer
 * Copyright (C) <2015> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_AUDIO_RESAMPLER_H__
#define __GST_AUDIO_RESAMPLER_H__

#include <gst/gst.h>
#include <gst/audio/audio.h>

G_BEGIN_DECLS

/**
 * GstAudioResampler:
 *
 * Opaque #GstAudioResampler struct.
 *
 * Since: 1.10
 */
typedef struct _GstAudioResampler GstAudioResampler;

/**
 * GST_AUDIO_RESAMPLER_OPT_CUTOFF:
 *
 * G_TYPE_DOUBLE, Cutoff parameter for the filter. 0.940 is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_CUTOFF      "GstAudioResampler.cutoff"
/**
 * GST_AUDIO_RESAMPLER_OPT_STOP_ATTENUATION:
 *
 * G_TYPE_DOUBLE, stopband attenuation in decibels. The attenuation
 * after the stopband for the kaiser window. 85 dB is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_STOP_ATTENUATION "GstAudioResampler.stop-attenutation"
/**
 * GST_AUDIO_RESAMPLER_OPT_TRANSITION_BANDWIDTH:
 *
 * G_TYPE_DOUBLE, transition bandwidth. The width of the
 * transition band for the kaiser window. 0.087 is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_TRANSITION_BANDWIDTH "GstAudioResampler.transition-bandwidth"

/**
 * GST_AUDIO_RESAMPLER_OPT_CUBIC_B:
 *
 * G_TYPE_DOUBLE, B parameter of the cubic filter.
 * Values between 0.0 and 2.0 are accepted. 1.0 is the default.
 *
 * Below are some values of popular filters:
 *                    B       C
 * Hermite           0.0     0.0
 * Spline            1.0     0.0
 * Catmull-Rom       0.0     1/2
 */
#define GST_AUDIO_RESAMPLER_OPT_CUBIC_B      "GstAudioResampler.cubic-b"
/**
 * GST_AUDIO_RESAMPLER_OPT_CUBIC_C:
 *
 * G_TYPE_DOUBLE, C parameter of the cubic filter.
 * Values between 0.0 and 2.0 are accepted. 0.0 is the default.
 *
 * See #GST_AUDIO_RESAMPLER_OPT_CUBIC_B for some more common values
 */
#define GST_AUDIO_RESAMPLER_OPT_CUBIC_C      "GstAudioResampler.cubic-c"

/**
 * GST_AUDIO_RESAMPLER_OPT_N_TAPS:
 *
 * G_TYPE_INT: the number of taps to use for the filter.
 * 0 is the default and selects the taps automatically.
 */
#define GST_AUDIO_RESAMPLER_OPT_N_TAPS      "GstAudioResampler.n-taps"

/**
 * GstAudioResamplerFilterMode:
 * @GST_AUDIO_RESAMPLER_FILTER_MODE_INTERPOLATED: Use interpolated filter tables. This
 *     uses less memory but more CPU and is slightly less accurate but it allows for more
 *     efficient variable rate resampling with gst_audio_resampler_update().
 * @GST_AUDIO_RESAMPLER_FILTER_MODE_FULL: Use full filter table. This uses more memory
 *     but less CPU.
 * @GST_AUDIO_RESAMPLER_FILTER_MODE_AUTO: Automatically choose between interpolated
 *     and full filter tables.
 *
 * Select for the filter tables should be set up.
 *
 * Since: 1.10
 */
typedef enum {
  GST_AUDIO_RESAMPLER_FILTER_MODE_INTERPOLATED = (0),
  GST_AUDIO_RESAMPLER_FILTER_MODE_FULL,
  GST_AUDIO_RESAMPLER_FILTER_MODE_AUTO,
} GstAudioResamplerFilterMode;
/**
 * GST_AUDIO_RESAMPLER_OPT_FILTER_MODE:
 *
 * GST_TYPE_AUDIO_RESAMPLER_FILTER_MODE: how the filter tables should be
 * constructed.
 * GST_AUDIO_RESAMPLER_FILTER_MODE_AUTO is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_FILTER_MODE      "GstAudioResampler.filter-mode"
/**
 * GST_AUDIO_RESAMPLER_OPT_FILTER_MODE_THRESHOLD:
 *
 * G_TYPE_UINT: the amount of memory to use for full filter tables before
 * switching to interpolated filter tables.
 * 1048576 is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_FILTER_MODE_THRESHOLD "GstAudioResampler.filter-mode-threshold"

/**
 * GstAudioResamplerFilterInterpolation:
 * @GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_NONE: no interpolation
 * @GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_LINEAR: linear interpolation of the
 *   filter coefficients.
 * @GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_CUBIC: cubic interpolation of the
 *   filter coefficients.
 *
 * The different filter interpolation methods.
 *
 * Since: 1.10
 */
typedef enum {
  GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_NONE = (0),
  GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_LINEAR,
  GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_CUBIC,
} GstAudioResamplerFilterInterpolation;
/**
 * GST_AUDIO_RESAMPLER_OPT_FILTER_INTERPOLATION:
 *
 * GST_TYPE_AUDIO_RESAMPLER_INTERPOLATION: how the filter coefficients should be
 *    interpolated.
 * GST_AUDIO_RESAMPLER_FILTER_INTERPOLATION_CUBIC is default.
 */
#define GST_AUDIO_RESAMPLER_OPT_FILTER_INTERPOLATION "GstAudioResampler.filter-interpolation"
/**
 * GST_AUDIO_RESAMPLER_OPT_FILTER_OVERSAMPLE:
 *
 * G_TYPE_UINT, oversampling to use when interpolating filters
 * 8 is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_FILTER_OVERSAMPLE "GstAudioResampler.filter-oversample"

/**
 * GST_AUDIO_RESAMPLER_OPT_MAX_PHASE_ERROR:
 *
 * G_TYPE_DOUBLE: The maximum allowed phase error when switching sample
 * rates.
 * 0.1 is the default.
 */
#define GST_AUDIO_RESAMPLER_OPT_MAX_PHASE_ERROR "GstAudioResampler.max-phase-error"

/**
 * GstAudioResamplerMethod:
 * @GST_AUDIO_RESAMPLER_METHOD_NEAREST: Duplicates the samples when
 *    upsampling and drops when downsampling
 * @GST_AUDIO_RESAMPLER_METHOD_LINEAR: Uses linear interpolation to reconstruct
 *    missing samples and averaging to downsample
 * @GST_AUDIO_RESAMPLER_METHOD_CUBIC: Uses cubic interpolation
 * @GST_AUDIO_RESAMPLER_METHOD_BLACKMAN_NUTTALL: Uses Blackman-Nuttall windowed sinc interpolation
 * @GST_AUDIO_RESAMPLER_METHOD_KAISER: Uses Kaiser windowed sinc interpolation
 *
 * Different subsampling and upsampling methods
 *
 * Since: 1.10
 */
typedef enum {
  GST_AUDIO_RESAMPLER_METHOD_NEAREST,
  GST_AUDIO_RESAMPLER_METHOD_LINEAR,
  GST_AUDIO_RESAMPLER_METHOD_CUBIC,
  GST_AUDIO_RESAMPLER_METHOD_BLACKMAN_NUTTALL,
  GST_AUDIO_RESAMPLER_METHOD_KAISER
} GstAudioResamplerMethod;

/**
 * GstAudioResamplerFlags:
 * @GST_AUDIO_RESAMPLER_FLAG_NONE: no flags
 * @GST_AUDIO_RESAMPLER_FLAG_NON_INTERLEAVED_IN: input samples are non-interleaved.
 *    an array of blocks of samples, one for each channel, should be passed to the
 *    resample function.
 * @GST_AUDIO_RESAMPLER_FLAG_NON_INTERLEAVED_OUT: output samples are non-interleaved.
 *    an array of blocks of samples, one for each channel, should be passed to the
 *    resample function.
 * @GST_AUDIO_RESAMPLER_FLAG_VARIABLE_RATE: optimize for dynamic updates of the sample
 *    rates with gst_audio_resampler_update(). This will select an interpolating filter
 *    when #GST_AUDIO_RESAMPLER_FILTER_MODE_AUTO is configured.
 *
 * Different resampler flags.
 *
 * Since: 1.10
 */
typedef enum {
  GST_AUDIO_RESAMPLER_FLAG_NONE                 = (0),
  GST_AUDIO_RESAMPLER_FLAG_NON_INTERLEAVED_IN   = (1 << 0),
  GST_AUDIO_RESAMPLER_FLAG_NON_INTERLEAVED_OUT  = (1 << 1),
  GST_AUDIO_RESAMPLER_FLAG_VARIABLE_RATE        = (1 << 2),
} GstAudioResamplerFlags;

#define GST_AUDIO_RESAMPLER_QUALITY_MIN 0
#define GST_AUDIO_RESAMPLER_QUALITY_MAX 10
#define GST_AUDIO_RESAMPLER_QUALITY_DEFAULT 4

GST_AUDIO_API
void           gst_audio_resampler_options_set_quality   (GstAudioResamplerMethod method,
                                                          guint quality,
                                                          gint in_rate, gint out_rate,
                                                          GstStructure *options);

GST_AUDIO_API
GstAudioResampler * gst_audio_resampler_new              (GstAudioResamplerMethod method,
                                                          GstAudioResamplerFlags flags,
                                                          GstAudioFormat format, gint channels,
                                                          gint in_rate, gint out_rate,
                                                          GstStructure *options);

GST_AUDIO_API
void                gst_audio_resampler_free             (GstAudioResampler *resampler);

GST_AUDIO_API
void                gst_audio_resampler_reset            (GstAudioResampler *resampler);

GST_AUDIO_API
gboolean            gst_audio_resampler_update           (GstAudioResampler *resampler,
                                                          gint in_rate, gint out_rate,
                                                          GstStructure *options);

GST_AUDIO_API
gsize               gst_audio_resampler_get_out_frames   (GstAudioResampler *resampler,
                                                          gsize in_frames);

GST_AUDIO_API
gsize               gst_audio_resampler_get_in_frames    (GstAudioResampler *resampler,
                                                          gsize out_frames);

GST_AUDIO_API
gsize               gst_audio_resampler_get_max_latency  (GstAudioResampler *resampler);

GST_AUDIO_API
void                gst_audio_resampler_resample         (GstAudioResampler * resampler,
                                                          gpointer in[], gsize in_frames,
                                                          gpointer out[], gsize out_frames);

G_END_DECLS

#endif /* __GST_AUDIO_RESAMPLER_H__ */
