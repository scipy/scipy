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

#ifndef __GST_FFT_S16_H__
#define __GST_FFT_S16_H__

#include <glib.h>
#include <gst/gst.h>

#include "gstfft.h"

G_BEGIN_DECLS

typedef struct _GstFFTS16 GstFFTS16;
typedef struct _GstFFTS16Complex GstFFTS16Complex;

/* Copy of kiss_fft_s16_cpx for documentation reasons,
 * do NOT change! */

/**
 * GstFFTS16Complex:
 * @r: Real part
 * @i: Imaginary part
 *
 * Data type for complex numbers composed of
 * signed 16 bit integers.
 */
struct _GstFFTS16Complex
{
  gint16 r;
  gint16 i;
};

/* Functions */

GST_FFT_API
GstFFTS16 *     gst_fft_s16_new         (gint len, gboolean inverse);

GST_FFT_API
void            gst_fft_s16_free        (GstFFTS16 *self);

GST_FFT_API
void            gst_fft_s16_fft         (GstFFTS16 *self, const gint16 *timedata,
                                         GstFFTS16Complex *freqdata);

GST_FFT_API
void            gst_fft_s16_inverse_fft (GstFFTS16 *self, const GstFFTS16Complex *freqdata,
                                         gint16 *timedata);

GST_FFT_API
void            gst_fft_s16_window      (GstFFTS16 *self, gint16 *timedata, GstFFTWindow window);

G_END_DECLS

#endif /* __GST_FFT_S16_H__ */
