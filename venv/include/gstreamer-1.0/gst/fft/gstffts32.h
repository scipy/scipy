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

#ifndef __GST_FFT_S32_H__
#define __GST_FFT_S32_H__

#include <glib.h>
#include <gst/gst.h>

#include "gstfft.h"

G_BEGIN_DECLS

typedef struct _GstFFTS32 GstFFTS32;
typedef struct _GstFFTS32Complex GstFFTS32Complex;

/* Copy of kiss_fft_s32_cpx for documentation reasons,
 * do NOT change! */

/**
 * GstFFTS32Complex:
 * @r: Real part
 * @i: Imaginary part
 *
 * Data type for complex numbers composed of
 * signed 32 bit integers.
 */
struct _GstFFTS32Complex
{
  gint32 r;
  gint32 i;
};

/* Functions */

GST_FFT_API
GstFFTS32 *     gst_fft_s32_new         (gint len, gboolean inverse);

GST_FFT_API
void            gst_fft_s32_free        (GstFFTS32 *self);

GST_FFT_API
void            gst_fft_s32_fft         (GstFFTS32 *self, const gint32 *timedata,
                                         GstFFTS32Complex *freqdata);

GST_FFT_API
void            gst_fft_s32_inverse_fft (GstFFTS32 *self, const GstFFTS32Complex *freqdata,
                                         gint32 *timedata);

GST_FFT_API
void            gst_fft_s32_window      (GstFFTS32 *self, gint32 *timedata, GstFFTWindow window);

G_END_DECLS

#endif /* __GST_FFT_S32_H__ */
