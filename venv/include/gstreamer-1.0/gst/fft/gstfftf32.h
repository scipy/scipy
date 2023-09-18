/* GStreamer
 * Copyright (C) <2007> Sebastian Dröge <slomo@circular-chaos.org>
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

#ifndef __GST_FFT_F32_H__
#define __GST_FFT_F32_H__

#include <glib.h>
#include <gst/gst.h>

#include "gstfft.h"

G_BEGIN_DECLS

typedef struct _GstFFTF32 GstFFTF32;
typedef struct _GstFFTF32Complex GstFFTF32Complex;

/* Copy of kiss_fft_f32_cpx for documentation reasons,
 * do NOT change! */

/**
 * GstFFTF32Complex:
 * @r: Real part
 * @i: Imaginary part
 *
 * Data type for complex numbers composed of
 * 32 bit float.
 */
struct _GstFFTF32Complex
{
  gfloat r;
  gfloat i;
};

/* Functions */

GST_FFT_API
GstFFTF32 *   gst_fft_f32_new           (gint len, gboolean inverse);

GST_FFT_API
void          gst_fft_f32_free          (GstFFTF32 *self);

GST_FFT_API
void          gst_fft_f32_fft           (GstFFTF32 *self, const gfloat *timedata,
                                         GstFFTF32Complex *freqdata);

GST_FFT_API
void          gst_fft_f32_inverse_fft   (GstFFTF32 *self, const GstFFTF32Complex *freqdata,
                                         gfloat *timedata);

GST_FFT_API
void          gst_fft_f32_window        (GstFFTF32 *self, gfloat *timedata, GstFFTWindow window);

G_END_DECLS

#endif /* __GST_FFT_F32_H__ */
