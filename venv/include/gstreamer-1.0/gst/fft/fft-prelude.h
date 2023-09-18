/* GStreamer FFT Library
 * Copyright (C) 2018 GStreamer developers
 *
 * fft-prelude.h: prelude include header for gst-fft library
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

#ifndef __GST_FFT_PRELUDE_H__
#define __GST_FFT_PRELUDE_H__

#include <gst/gst.h>

#ifdef BUILDING_GST_FFT
#define GST_FFT_API GST_API_EXPORT         /* from config.h */
#else
#define GST_FFT_API GST_API_IMPORT
#endif

#endif /* __GST_FFT_PRELUDE_H__ */
