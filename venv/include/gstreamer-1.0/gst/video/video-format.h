/* GStreamer
 * Copyright (C) <2011> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_VIDEO_FORMAT_H__
#define __GST_VIDEO_FORMAT_H__

#include <gst/gst.h>

G_BEGIN_DECLS

#include <gst/video/video-enumtypes.h>
#include <gst/video/video-tile.h>

/**
 * GstVideoFormat:
 * @GST_VIDEO_FORMAT_UNKNOWN: Unknown or unset video format id
 * @GST_VIDEO_FORMAT_ENCODED: Encoded video format. Only ever use that in caps for
 *                            special video formats in combination with non-system
 *                            memory GstCapsFeatures where it does not make sense
 *                            to specify a real video format.
 * @GST_VIDEO_FORMAT_I420: planar 4:2:0 YUV
 * @GST_VIDEO_FORMAT_YV12: planar 4:2:0 YVU (like I420 but UV planes swapped)
 * @GST_VIDEO_FORMAT_YUY2: packed 4:2:2 YUV (Y0-U0-Y1-V0 Y2-U2-Y3-V2 Y4 ...)
 * @GST_VIDEO_FORMAT_UYVY: packed 4:2:2 YUV (U0-Y0-V0-Y1 U2-Y2-V2-Y3 U4 ...)
 * @GST_VIDEO_FORMAT_VYUY: packed 4:2:2 YUV (V0-Y0-U0-Y1 V2-Y2-U2-Y3 V4 ...)
 * @GST_VIDEO_FORMAT_AYUV: packed 4:4:4 YUV with alpha channel (A0-Y0-U0-V0 ...)
 * @GST_VIDEO_FORMAT_RGBx: sparse rgb packed into 32 bit, space last
 * @GST_VIDEO_FORMAT_BGRx: sparse reverse rgb packed into 32 bit, space last
 * @GST_VIDEO_FORMAT_xRGB: sparse rgb packed into 32 bit, space first
 * @GST_VIDEO_FORMAT_xBGR: sparse reverse rgb packed into 32 bit, space first
 * @GST_VIDEO_FORMAT_RGBA: rgb with alpha channel last
 * @GST_VIDEO_FORMAT_BGRA: reverse rgb with alpha channel last
 * @GST_VIDEO_FORMAT_ARGB: rgb with alpha channel first
 * @GST_VIDEO_FORMAT_ABGR: reverse rgb with alpha channel first
 * @GST_VIDEO_FORMAT_RGB: RGB packed into 24 bits without padding (`R-G-B-R-G-B`)
 * @GST_VIDEO_FORMAT_BGR: reverse RGB packed into 24 bits without padding (`B-G-R-B-G-R`)
 * @GST_VIDEO_FORMAT_Y41B: planar 4:1:1 YUV
 * @GST_VIDEO_FORMAT_Y42B: planar 4:2:2 YUV
 * @GST_VIDEO_FORMAT_YVYU: packed 4:2:2 YUV (Y0-V0-Y1-U0 Y2-V2-Y3-U2 Y4 ...)
 * @GST_VIDEO_FORMAT_Y444: planar 4:4:4 YUV
 * @GST_VIDEO_FORMAT_v210: packed 4:2:2 10-bit YUV, complex format
 * @GST_VIDEO_FORMAT_v216: packed 4:2:2 16-bit YUV, Y0-U0-Y1-V1 order
 * @GST_VIDEO_FORMAT_NV12: planar 4:2:0 YUV with interleaved UV plane
 * @GST_VIDEO_FORMAT_NV21: planar 4:2:0 YUV with interleaved VU plane
 * @GST_VIDEO_FORMAT_NV12_10LE32: 10-bit variant of @GST_VIDEO_FORMAT_NV12, packed into 32bit words (MSB 2 bits padding) (Since: 1.14)
 * @GST_VIDEO_FORMAT_GRAY8: 8-bit grayscale
 * @GST_VIDEO_FORMAT_GRAY10_LE32: 10-bit grayscale, packed into 32bit words (2 bits padding) (Since: 1.14)
 * @GST_VIDEO_FORMAT_GRAY16_BE: 16-bit grayscale, most significant byte first
 * @GST_VIDEO_FORMAT_GRAY16_LE: 16-bit grayscale, least significant byte first
 * @GST_VIDEO_FORMAT_v308: packed 4:4:4 YUV (Y-U-V ...)
 * @GST_VIDEO_FORMAT_IYU2: packed 4:4:4 YUV (U-Y-V ...) (Since: 1.10)
 * @GST_VIDEO_FORMAT_RGB16: rgb 5-6-5 bits per component
 * @GST_VIDEO_FORMAT_BGR16: reverse rgb 5-6-5 bits per component
 * @GST_VIDEO_FORMAT_RGB15: rgb 5-5-5 bits per component
 * @GST_VIDEO_FORMAT_BGR15: reverse rgb 5-5-5 bits per component
 * @GST_VIDEO_FORMAT_UYVP: packed 10-bit 4:2:2 YUV (U0-Y0-V0-Y1 U2-Y2-V2-Y3 U4 ...)
 * @GST_VIDEO_FORMAT_A420: planar 4:4:2:0 AYUV
 * @GST_VIDEO_FORMAT_RGB8P: 8-bit paletted RGB
 * @GST_VIDEO_FORMAT_YUV9: planar 4:1:0 YUV
 * @GST_VIDEO_FORMAT_YVU9: planar 4:1:0 YUV (like YUV9 but UV planes swapped)
 * @GST_VIDEO_FORMAT_IYU1: packed 4:1:1 YUV (Cb-Y0-Y1-Cr-Y2-Y3 ...)
 * @GST_VIDEO_FORMAT_ARGB64: rgb with alpha channel first, 16 bits (native endianness) per channel
 * @GST_VIDEO_FORMAT_AYUV64: packed 4:4:4 YUV with alpha channel, 16 bits (native endianness) per channel (A0-Y0-U0-V0 ...)
 * @GST_VIDEO_FORMAT_r210: packed 4:4:4 RGB, 10 bits per channel
 * @GST_VIDEO_FORMAT_I420_10BE: planar 4:2:0 YUV, 10 bits per channel
 * @GST_VIDEO_FORMAT_I420_10LE: planar 4:2:0 YUV, 10 bits per channel
 * @GST_VIDEO_FORMAT_I422_10BE: planar 4:2:2 YUV, 10 bits per channel
 * @GST_VIDEO_FORMAT_I422_10LE: planar 4:2:2 YUV, 10 bits per channel
 * @GST_VIDEO_FORMAT_Y444_10BE: planar 4:4:4 YUV, 10 bits per channel (Since: 1.2)
 * @GST_VIDEO_FORMAT_Y444_10LE: planar 4:4:4 YUV, 10 bits per channel (Since: 1.2)
 * @GST_VIDEO_FORMAT_GBR: planar 4:4:4 RGB, 8 bits per channel (Since: 1.2)
 * @GST_VIDEO_FORMAT_GBR_10BE: planar 4:4:4 RGB, 10 bits per channel (Since: 1.2)
 * @GST_VIDEO_FORMAT_GBR_10LE: planar 4:4:4 RGB, 10 bits per channel (Since: 1.2)
 * @GST_VIDEO_FORMAT_NV16: planar 4:2:2 YUV with interleaved UV plane (Since: 1.2)
 * @GST_VIDEO_FORMAT_NV16_10LE32: 10-bit variant of @GST_VIDEO_FORMAT_NV16, packed into 32bit words (MSB 2 bits padding) (Since: 1.14)
 * @GST_VIDEO_FORMAT_NV24: planar 4:4:4 YUV with interleaved UV plane (Since: 1.2)
 * @GST_VIDEO_FORMAT_NV12_64Z32: NV12 with 64x32 tiling in zigzag pattern (Since: 1.4)
 * @GST_VIDEO_FORMAT_A420_10BE: planar 4:4:2:0 YUV, 10 bits per channel (Since: 1.6)
 * @GST_VIDEO_FORMAT_A420_10LE: planar 4:4:2:0 YUV, 10 bits per channel (Since: 1.6)
 * @GST_VIDEO_FORMAT_A422_10BE: planar 4:4:2:2 YUV, 10 bits per channel (Since: 1.6)
 * @GST_VIDEO_FORMAT_A422_10LE: planar 4:4:2:2 YUV, 10 bits per channel (Since: 1.6)
 * @GST_VIDEO_FORMAT_A444_10BE: planar 4:4:4:4 YUV, 10 bits per channel (Since: 1.6)
 * @GST_VIDEO_FORMAT_A444_10LE: planar 4:4:4:4 YUV, 10 bits per channel (Since: 1.6)
 * @GST_VIDEO_FORMAT_NV61: planar 4:2:2 YUV with interleaved VU plane (Since: 1.6)
 * @GST_VIDEO_FORMAT_P010_10BE: planar 4:2:0 YUV with interleaved UV plane, 10 bits per channel (Since: 1.10)
 * @GST_VIDEO_FORMAT_P010_10LE: planar 4:2:0 YUV with interleaved UV plane, 10 bits per channel (Since: 1.10)
 * @GST_VIDEO_FORMAT_GBRA: planar 4:4:4:4 ARGB, 8 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_GBRA_10BE: planar 4:4:4:4 ARGB, 10 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_GBRA_10LE: planar 4:4:4:4 ARGB, 10 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_GBR_12BE: planar 4:4:4 RGB, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_GBR_12LE: planar 4:4:4 RGB, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_GBRA_12BE: planar 4:4:4:4 ARGB, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_GBRA_12LE: planar 4:4:4:4 ARGB, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_I420_12BE: planar 4:2:0 YUV, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_I420_12LE: planar 4:2:0 YUV, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_I422_12BE: planar 4:2:2 YUV, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_I422_12LE: planar 4:2:2 YUV, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_Y444_12BE: planar 4:4:4 YUV, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_Y444_12LE: planar 4:4:4 YUV, 12 bits per channel (Since: 1.12)
 * @GST_VIDEO_FORMAT_NV12_10LE40: Fully packed variant of NV12_10LE32 (Since: 1.16)
 * @GST_VIDEO_FORMAT_Y210: packed 4:2:2 YUV, 10 bits per channel (Since: 1.16)
 * @GST_VIDEO_FORMAT_Y410: packed 4:4:4 YUV, 10 bits per channel(A-V-Y-U...) (Since: 1.16)
 * @GST_VIDEO_FORMAT_VUYA: packed 4:4:4 YUV with alpha channel (V0-U0-Y0-A0...) (Since: 1.16)
 * @GST_VIDEO_FORMAT_BGR10A2_LE: packed 4:4:4 RGB with alpha channel(B-G-R-A), 10 bits for R/G/B channel and MSB 2 bits for alpha channel (Since: 1.16)
 * @GST_VIDEO_FORMAT_RGB10A2_LE: packed 4:4:4 RGB with alpha channel(R-G-B-A), 10 bits for R/G/B channel and MSB 2 bits for alpha channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_Y444_16BE: planar 4:4:4 YUV, 16 bits per channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_Y444_16LE: planar 4:4:4 YUV, 16 bits per channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_P016_BE: planar 4:2:0 YUV with interleaved UV plane, 16 bits per channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_P016_LE: planar 4:2:0 YUV with interleaved UV plane, 16 bits per channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_P012_BE: planar 4:2:0 YUV with interleaved UV plane, 12 bits per channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_P012_LE: planar 4:2:0 YUV with interleaved UV plane, 12 bits per channel (Since: 1.18)
 * @GST_VIDEO_FORMAT_Y212_BE: packed 4:2:2 YUV, 12 bits per channel (Y-U-Y-V) (Since: 1.18)
 * @GST_VIDEO_FORMAT_Y212_LE: packed 4:2:2 YUV, 12 bits per channel (Y-U-Y-V) (Since: 1.18)
 * @GST_VIDEO_FORMAT_Y412_BE: packed 4:4:4:4 YUV, 12 bits per channel(U-Y-V-A...) (Since: 1.18)
 * @GST_VIDEO_FORMAT_Y412_LE: packed 4:4:4:4 YUV, 12 bits per channel(U-Y-V-A...) (Since: 1.18)
 * @GST_VIDEO_FORMAT_NV12_4L4: NV12 with 4x4 tiles in linear order (Since: 1.18)
 * @GST_VIDEO_FORMAT_NV12_32L32: NV12 with 32x32 tiles in linear order (Since: 1.18)
 * @GST_VIDEO_FORMAT_RGBP: planar 4:4:4 RGB, 8 bits per channel (Since: 1.20)
 * @GST_VIDEO_FORMAT_BGRP: planar 4:4:4 RGB, 8 bits per channel (Since: 1.20)
 * @GST_VIDEO_FORMAT_AV12: Planar 4:2:0 YUV with interleaved UV plane with alpha as 3rd plane (Since: 1.20)
 * @GST_VIDEO_FORMAT_ARGB64_LE: RGB with alpha channel first, 16 bits per channel
 * @GST_VIDEO_FORMAT_ARGB64_BE: RGB with alpha channel first, 16 bits per channel
 * @GST_VIDEO_FORMAT_RGBA64_LE: RGB with alpha channel last, 16 bits per channel
 * @GST_VIDEO_FORMAT_RGBA64_BE: RGB with alpha channel last, 16 bits per channel
 * @GST_VIDEO_FORMAT_BGRA64_LE: reverse RGB with alpha channel last, 16 bits per channel
 * @GST_VIDEO_FORMAT_BGRA64_BE: reverse RGB with alpha channel last, 16 bits per channel
 * @GST_VIDEO_FORMAT_ABGR64_LE: reverse RGB with alpha channel first, 16 bits per channel
 * @GST_VIDEO_FORMAT_ABGR64_BE: reverse RGB with alpha channel first, 16 bits per channel
 * @GST_VIDEO_FORMAT_NV12_16L32S: NV12 with 16x32 Y tiles and 16x16 UV tiles. (Since: 1.22)
 * @GST_VIDEO_FORMAT_NV12_8L128 : NV12 with 8x128 tiles in linear order (Since: 1.22)
 * @GST_VIDEO_FORMAT_NV12_10BE_8L128 : NV12 10bit big endian with 8x128 tiles in linear order (Since: 1.22)
 *
 * Enum value describing the most common video formats.
 *
 * See the [GStreamer raw video format design document](https://gstreamer.freedesktop.org/documentation/additional/design/mediatype-video-raw.html#formats)
 * for details about the layout and packing of these formats in memory.
 */
typedef enum {
  GST_VIDEO_FORMAT_UNKNOWN,
  GST_VIDEO_FORMAT_ENCODED,
  GST_VIDEO_FORMAT_I420,
  GST_VIDEO_FORMAT_YV12,
  GST_VIDEO_FORMAT_YUY2,
  GST_VIDEO_FORMAT_UYVY,
  GST_VIDEO_FORMAT_AYUV,
  GST_VIDEO_FORMAT_RGBx,
  GST_VIDEO_FORMAT_BGRx,
  GST_VIDEO_FORMAT_xRGB,
  GST_VIDEO_FORMAT_xBGR,
  GST_VIDEO_FORMAT_RGBA,
  GST_VIDEO_FORMAT_BGRA,
  GST_VIDEO_FORMAT_ARGB,
  GST_VIDEO_FORMAT_ABGR,
  GST_VIDEO_FORMAT_RGB,
  GST_VIDEO_FORMAT_BGR,
  GST_VIDEO_FORMAT_Y41B,
  GST_VIDEO_FORMAT_Y42B,
  GST_VIDEO_FORMAT_YVYU,
  GST_VIDEO_FORMAT_Y444,
  GST_VIDEO_FORMAT_v210,
  GST_VIDEO_FORMAT_v216,
  GST_VIDEO_FORMAT_NV12,
  GST_VIDEO_FORMAT_NV21,
  GST_VIDEO_FORMAT_GRAY8,
  GST_VIDEO_FORMAT_GRAY16_BE,
  GST_VIDEO_FORMAT_GRAY16_LE,
  GST_VIDEO_FORMAT_v308,
  GST_VIDEO_FORMAT_RGB16,
  GST_VIDEO_FORMAT_BGR16,
  GST_VIDEO_FORMAT_RGB15,
  GST_VIDEO_FORMAT_BGR15,
  GST_VIDEO_FORMAT_UYVP,
  GST_VIDEO_FORMAT_A420,
  GST_VIDEO_FORMAT_RGB8P,
  GST_VIDEO_FORMAT_YUV9,
  GST_VIDEO_FORMAT_YVU9,
  GST_VIDEO_FORMAT_IYU1,
  GST_VIDEO_FORMAT_ARGB64,
  GST_VIDEO_FORMAT_AYUV64,
  GST_VIDEO_FORMAT_r210,
  GST_VIDEO_FORMAT_I420_10BE,
  GST_VIDEO_FORMAT_I420_10LE,
  GST_VIDEO_FORMAT_I422_10BE,
  GST_VIDEO_FORMAT_I422_10LE,
  GST_VIDEO_FORMAT_Y444_10BE,
  GST_VIDEO_FORMAT_Y444_10LE,
  GST_VIDEO_FORMAT_GBR,
  GST_VIDEO_FORMAT_GBR_10BE,
  GST_VIDEO_FORMAT_GBR_10LE,
  GST_VIDEO_FORMAT_NV16,
  GST_VIDEO_FORMAT_NV24,
  GST_VIDEO_FORMAT_NV12_64Z32,
  GST_VIDEO_FORMAT_A420_10BE,
  GST_VIDEO_FORMAT_A420_10LE,
  GST_VIDEO_FORMAT_A422_10BE,
  GST_VIDEO_FORMAT_A422_10LE,
  GST_VIDEO_FORMAT_A444_10BE,
  GST_VIDEO_FORMAT_A444_10LE,
  GST_VIDEO_FORMAT_NV61,
  GST_VIDEO_FORMAT_P010_10BE,
  GST_VIDEO_FORMAT_P010_10LE,
  GST_VIDEO_FORMAT_IYU2,
  GST_VIDEO_FORMAT_VYUY,
  GST_VIDEO_FORMAT_GBRA,
  GST_VIDEO_FORMAT_GBRA_10BE,
  GST_VIDEO_FORMAT_GBRA_10LE,
  GST_VIDEO_FORMAT_GBR_12BE,
  GST_VIDEO_FORMAT_GBR_12LE,
  GST_VIDEO_FORMAT_GBRA_12BE,
  GST_VIDEO_FORMAT_GBRA_12LE,
  GST_VIDEO_FORMAT_I420_12BE,
  GST_VIDEO_FORMAT_I420_12LE,
  GST_VIDEO_FORMAT_I422_12BE,
  GST_VIDEO_FORMAT_I422_12LE,
  GST_VIDEO_FORMAT_Y444_12BE,
  GST_VIDEO_FORMAT_Y444_12LE,
  GST_VIDEO_FORMAT_GRAY10_LE32,
  GST_VIDEO_FORMAT_NV12_10LE32,
  GST_VIDEO_FORMAT_NV16_10LE32,
  GST_VIDEO_FORMAT_NV12_10LE40,
  GST_VIDEO_FORMAT_Y210,
  GST_VIDEO_FORMAT_Y410,
  GST_VIDEO_FORMAT_VUYA,
  GST_VIDEO_FORMAT_BGR10A2_LE,
  GST_VIDEO_FORMAT_RGB10A2_LE,
  GST_VIDEO_FORMAT_Y444_16BE,
  GST_VIDEO_FORMAT_Y444_16LE,
  GST_VIDEO_FORMAT_P016_BE,
  GST_VIDEO_FORMAT_P016_LE,
  GST_VIDEO_FORMAT_P012_BE,
  GST_VIDEO_FORMAT_P012_LE,
  GST_VIDEO_FORMAT_Y212_BE,
  GST_VIDEO_FORMAT_Y212_LE,
  GST_VIDEO_FORMAT_Y412_BE,
  GST_VIDEO_FORMAT_Y412_LE,
  /**
   * GST_VIDEO_FORMAT_NV12_4L4:
   *
   * NV12 with 4x4 tiles in linear order.
   *
   * Since: 1.18
   */
  GST_VIDEO_FORMAT_NV12_4L4,
  /**
   * GST_VIDEO_FORMAT_NV12_32L32:
   *
   * NV12 with 32x32 tiles in linear order.
   *
   * Since: 1.18
   */
  GST_VIDEO_FORMAT_NV12_32L32,

  /**
   * GST_VIDEO_FORMAT_RGBP:
   *
   * Planar 4:4:4 RGB, R-G-B order
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_RGBP,

  /**
   * GST_VIDEO_FORMAT_BGRP:
   *
   * Planar 4:4:4 RGB, B-G-R order
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_BGRP,

  /**
   * GST_VIDEO_FORMAT_AV12:
   *
   * Planar 4:2:0 YUV with interleaved UV plane with alpha as
   * 3rd plane.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_AV12,

  /**
   * GST_VIDEO_FORMAT_ARGB64_LE:
   *
   * RGB with alpha channel first, 16 bits (little endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_ARGB64_LE,

  /**
   * GST_VIDEO_FORMAT_ARGB64_BE:
   *
   * RGB with alpha channel first, 16 bits (big endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_ARGB64_BE,

  /**
   * GST_VIDEO_FORMAT_RGBA64_LE:
   *
   * RGB with alpha channel last, 16 bits (little endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_RGBA64_LE,

  /**
   * GST_VIDEO_FORMAT_RGBA64_BE:
   *
   * RGB with alpha channel last, 16 bits (big endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_RGBA64_BE,

  /**
   * GST_VIDEO_FORMAT_BGRA64_LE:
   *
   * Reverse RGB with alpha channel last, 16 bits (little endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_BGRA64_LE,

  /**
   * GST_VIDEO_FORMAT_BGRA64_BE:
   *
   * Reverse RGB with alpha channel last, 16 bits (big endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_BGRA64_BE,

  /**
   * GST_VIDEO_FORMAT_ABGR64_LE:
   *
   * Reverse RGB with alpha channel first, 16 bits (little endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_ABGR64_LE,

  /**
   * GST_VIDEO_FORMAT_ABGR64_BE:
   *
   * Reverse RGB with alpha channel first, 16 bits (big endian)
   * per channel.
   *
   * Since: 1.20
   */
  GST_VIDEO_FORMAT_ABGR64_BE,

  /**
   * GST_VIDEO_FORMAT_NV12_16L32S:
   *
   * NV12 with 16x32 Y tiles and 16x16 UV tiles.
   *
   * Since: 1.22
   */
  GST_VIDEO_FORMAT_NV12_16L32S,

  /**
   * GST_VIDEO_FORMAT_NV12_8L128:
   *
   * NV12 with 8x128 tiles in linear order.
   *
   * Since: 1.22
   */
  GST_VIDEO_FORMAT_NV12_8L128,

  /**
   * GST_VIDEO_FORMAT_NV12_10BE_8L128:
   *
   * NV12 10bit big endian with 8x128 tiles in linear order.
   *
   * Since: 1.22
   */
  GST_VIDEO_FORMAT_NV12_10BE_8L128,
} GstVideoFormat;

#define GST_VIDEO_MAX_PLANES 4
#define GST_VIDEO_MAX_COMPONENTS 4

typedef struct _GstVideoFormatInfo GstVideoFormatInfo;

/**
 * GstVideoFormatFlags:
 * @GST_VIDEO_FORMAT_FLAG_YUV: The video format is YUV, components are numbered
 *   0=Y, 1=U, 2=V.
 * @GST_VIDEO_FORMAT_FLAG_RGB: The video format is RGB, components are numbered
 *   0=R, 1=G, 2=B.
 * @GST_VIDEO_FORMAT_FLAG_GRAY: The video is gray, there is one gray component
 *   with index 0.
 * @GST_VIDEO_FORMAT_FLAG_ALPHA: The video format has an alpha components with
 *   the number 3.
 * @GST_VIDEO_FORMAT_FLAG_LE: The video format has data stored in little
 *   endianness.
 * @GST_VIDEO_FORMAT_FLAG_PALETTE: The video format has a palette. The palette
 *   is stored in the second plane and indexes are stored in the first plane.
 * @GST_VIDEO_FORMAT_FLAG_COMPLEX: The video format has a complex layout that
 *   can't be described with the usual information in the #GstVideoFormatInfo.
 * @GST_VIDEO_FORMAT_FLAG_UNPACK: This format can be used in a
 *   #GstVideoFormatUnpack and #GstVideoFormatPack function.
 * @GST_VIDEO_FORMAT_FLAG_TILED: The format is tiled, there is tiling information
 *   in the last plane.
 * @GST_VIDEO_FORMAT_FLAG_SUBTILES: The tile size varies per plane
 *   according to the subsampling. (Since: 1.22)
 *
 * The different video flags that a format info can have.
 */
typedef enum
{
  GST_VIDEO_FORMAT_FLAG_YUV      = (1 << 0),
  GST_VIDEO_FORMAT_FLAG_RGB      = (1 << 1),
  GST_VIDEO_FORMAT_FLAG_GRAY     = (1 << 2),
  GST_VIDEO_FORMAT_FLAG_ALPHA    = (1 << 3),
  GST_VIDEO_FORMAT_FLAG_LE       = (1 << 4),
  GST_VIDEO_FORMAT_FLAG_PALETTE  = (1 << 5),
  GST_VIDEO_FORMAT_FLAG_COMPLEX  = (1 << 6),
  GST_VIDEO_FORMAT_FLAG_UNPACK   = (1 << 7),
  GST_VIDEO_FORMAT_FLAG_TILED    = (1 << 8),
  /**
   * GST_VIDEO_FORMAT_FLAG_SUBTILES:
   *
   * The tile size varies per plane according to the subsampling.
   *
   * Since: 1.22
   */
  GST_VIDEO_FORMAT_FLAG_SUBTILES = (1 << 9)
} GstVideoFormatFlags;

/* YUV components */
#define GST_VIDEO_COMP_Y  0
#define GST_VIDEO_COMP_U  1
#define GST_VIDEO_COMP_V  2

/* RGB components */
#define GST_VIDEO_COMP_R  0
#define GST_VIDEO_COMP_G  1
#define GST_VIDEO_COMP_B  2

/* alpha component */
#define GST_VIDEO_COMP_A  3

/* palette components */
#define GST_VIDEO_COMP_INDEX    0
#define GST_VIDEO_COMP_PALETTE  1

#include <gst/video/video-chroma.h>

/**
 * GstVideoPackFlags:
 * @GST_VIDEO_PACK_FLAG_NONE: No flag
 * @GST_VIDEO_PACK_FLAG_TRUNCATE_RANGE: When the source has a smaller depth
 *   than the target format, set the least significant bits of the target
 *   to 0. This is likely slightly faster but less accurate. When this flag
 *   is not specified, the most significant bits of the source are duplicated
 *   in the least significant bits of the destination.
 * @GST_VIDEO_PACK_FLAG_INTERLACED: The source is interlaced. The unpacked
 *   format will be interlaced as well with each line containing
 *   information from alternating fields. (Since: 1.2)
 *
 * The different flags that can be used when packing and unpacking.
 */
typedef enum
{
  GST_VIDEO_PACK_FLAG_NONE           = 0,
  GST_VIDEO_PACK_FLAG_TRUNCATE_RANGE = (1 << 0),
  GST_VIDEO_PACK_FLAG_INTERLACED     = (1 << 1)
} GstVideoPackFlags;

/**
 * GstVideoFormatUnpack:
 * @info: a #GstVideoFormatInfo
 * @flags: flags to control the unpacking
 * @dest: a destination array
 * @data: pointers to the data planes
 * @stride: strides of the planes
 * @x: the x position in the image to start from
 * @y: the y position in the image to start from
 * @width: the amount of pixels to unpack.
 *
 * Unpacks @width pixels from the given planes and strides containing data of
 * format @info. The pixels will be unpacked into @dest with each component
 * interleaved as per @info's unpack_format, which will usually be one of
 * #GST_VIDEO_FORMAT_ARGB, #GST_VIDEO_FORMAT_AYUV, #GST_VIDEO_FORMAT_ARGB64 or
 * #GST_VIDEO_FORMAT_AYUV64 depending on the format to unpack.
 * @dest should at least be big enough to hold @width * bytes_per_pixel bytes
 * where bytes_per_pixel relates to the unpack format and will usually be
 * either 4 or 8 depending on the unpack format. bytes_per_pixel will be
 * the same as the pixel stride for plane 0 for the above formats.
 *
 * For subsampled formats, the components will be duplicated in the destination
 * array. Reconstruction of the missing components can be performed in a
 * separate step after unpacking.
 */
typedef void (*GstVideoFormatUnpack)         (const GstVideoFormatInfo *info,
                                              GstVideoPackFlags flags, gpointer dest,
                                              const gpointer data[GST_VIDEO_MAX_PLANES],
                                              const gint stride[GST_VIDEO_MAX_PLANES],
                                              gint x, gint y, gint width);
/**
 * GstVideoFormatPack:
 * @info: a #GstVideoFormatInfo
 * @flags: flags to control the packing
 * @src: a source array
 * @sstride: the source array stride
 * @data: pointers to the destination data planes
 * @stride: strides of the destination planes
 * @chroma_site: the chroma siting of the target when subsampled (not used)
 * @y: the y position in the image to pack to
 * @width: the amount of pixels to pack.
 *
 * Packs @width pixels from @src to the given planes and strides in the
 * format @info. The pixels from source have each component interleaved
 * and will be packed into the planes in @data.
 *
 * This function operates on pack_lines lines, meaning that @src should
 * contain at least pack_lines lines with a stride of @sstride and @y
 * should be a multiple of pack_lines.
 *
 * Subsampled formats will use the horizontally and vertically cosited
 * component from the source. Subsampling should be performed before
 * packing.
 *
 * Because this function does not have a x coordinate, it is not possible to
 * pack pixels starting from an unaligned position. For tiled images this
 * means that packing should start from a tile coordinate. For subsampled
 * formats this means that a complete pixel needs to be packed.
 */
/* FIXME(2.0): remove the chroma_site, it is unused and is not relevant for
 * packing, chroma subsampling based on chroma-site should be done in a separate
 * step before packing*/
typedef void (*GstVideoFormatPack)           (const GstVideoFormatInfo *info,
                                              GstVideoPackFlags flags,
                                              const gpointer src, gint sstride,
                                              gpointer data[GST_VIDEO_MAX_PLANES],
                                              const gint stride[GST_VIDEO_MAX_PLANES],
                                              GstVideoChromaSite chroma_site,
                                              gint y, gint width);

/**
 * GstVideoFormatInfo:
 * @format: #GstVideoFormat
 * @name: string representation of the format
 * @description: use readable description of the format
 * @flags: #GstVideoFormatFlags
 * @bits: The number of bits used to pack data items. This can be less than 8
 *    when multiple pixels are stored in a byte. for values > 8 multiple bytes
 *    should be read according to the endianness flag before applying the shift
 *    and mask.
 * @n_components: the number of components in the video format.
 * @shift: the number of bits to shift away to get the component data
 * @depth: the depth in bits for each component
 * @pixel_stride: the pixel stride of each component. This is the amount of
 *    bytes to the pixel immediately to the right. When bits < 8, the stride is
 *    expressed in bits. For 24-bit RGB, this would be 3 bytes, for example,
 *    while it would be 4 bytes for RGBx or ARGB.
 * @n_planes: the number of planes for this format. The number of planes can be
 *    less than the amount of components when multiple components are packed into
 *    one plane.
 * @plane: the plane number where a component can be found
 * @poffset: the offset in the plane where the first pixel of the components
 *    can be found.
 * @w_sub: subsampling factor of the width for the component. Use
 *     GST_VIDEO_SUB_SCALE to scale a width.
 * @h_sub: subsampling factor of the height for the component. Use
 *     GST_VIDEO_SUB_SCALE to scale a height.
 * @unpack_format: the format of the unpacked pixels. This format must have the
 *     #GST_VIDEO_FORMAT_FLAG_UNPACK flag set.
 * @unpack_func: an unpack function for this format
 * @pack_lines: the amount of lines that will be packed
 * @pack_func: an pack function for this format
 * @tile_mode: The tiling mode
 * @tile_ws: The width of a tile, in bytes, represented as a shift. DEPRECATED,
 * use tile_info[] array instead.
 * @tile_hs: The height of a tile, in bytes, represented as a shift. DEPREACTED,
 * use tile_info[] array instead.
 * @tile_info: Per-plane tile information
 *
 * Information for a video format.
 */
struct _GstVideoFormatInfo {
  GstVideoFormat format;
  const gchar *name;
  const gchar *description;
  GstVideoFormatFlags flags;
  guint bits;
  guint n_components;
  guint shift[GST_VIDEO_MAX_COMPONENTS];
  guint depth[GST_VIDEO_MAX_COMPONENTS];
  gint  pixel_stride[GST_VIDEO_MAX_COMPONENTS];
  guint n_planes;
  guint plane[GST_VIDEO_MAX_COMPONENTS];
  guint poffset[GST_VIDEO_MAX_COMPONENTS];
  guint w_sub[GST_VIDEO_MAX_COMPONENTS];
  guint h_sub[GST_VIDEO_MAX_COMPONENTS];

  GstVideoFormat unpack_format;
  GstVideoFormatUnpack unpack_func;
  gint pack_lines;
  GstVideoFormatPack pack_func;

  GstVideoTileMode tile_mode;
  G_DEPRECATED_FOR(tile_info) guint tile_ws;
  G_DEPRECATED_FOR(tile_info) guint tile_hs;

  /**
   * GstVideoFormatInfo.tile_info:
   *
   * Information about the tiles for each of the planes.
   *
   * Since: 1.22
   */
  GstVideoTileInfo tile_info[GST_VIDEO_MAX_PLANES];
};

/**
 * GST_VIDEO_FORMAT_INFO_IS_VALID_RAW:
 *
 * Tests that the given #GstVideoFormatInfo represents a valid un-encoded
 * format.
 *
 * Since: 1.22
 */
#define GST_VIDEO_FORMAT_INFO_IS_VALID_RAW(info)              \
  (info != NULL && (info)->format > GST_VIDEO_FORMAT_ENCODED)

#define GST_VIDEO_FORMAT_INFO_FORMAT(info)       ((info)->format)
#define GST_VIDEO_FORMAT_INFO_NAME(info)         ((info)->name)
#define GST_VIDEO_FORMAT_INFO_FLAGS(info)        ((info)->flags)

#define GST_VIDEO_FORMAT_INFO_IS_YUV(info)       (((info)->flags & GST_VIDEO_FORMAT_FLAG_YUV) != 0)
#define GST_VIDEO_FORMAT_INFO_IS_RGB(info)       (((info)->flags & GST_VIDEO_FORMAT_FLAG_RGB) != 0)
#define GST_VIDEO_FORMAT_INFO_IS_GRAY(info)      (((info)->flags & GST_VIDEO_FORMAT_FLAG_GRAY) != 0)
#define GST_VIDEO_FORMAT_INFO_HAS_ALPHA(info)    (((info)->flags & GST_VIDEO_FORMAT_FLAG_ALPHA) != 0)
#define GST_VIDEO_FORMAT_INFO_IS_LE(info)        (((info)->flags & GST_VIDEO_FORMAT_FLAG_LE) != 0)
#define GST_VIDEO_FORMAT_INFO_HAS_PALETTE(info)  (((info)->flags & GST_VIDEO_FORMAT_FLAG_PALETTE) != 0)
#define GST_VIDEO_FORMAT_INFO_IS_COMPLEX(info)   (((info)->flags & GST_VIDEO_FORMAT_FLAG_COMPLEX) != 0)
#define GST_VIDEO_FORMAT_INFO_IS_TILED(info)     (((info)->flags & GST_VIDEO_FORMAT_FLAG_TILED) != 0)
/**
 * GST_VIDEO_FORMAT_INFO_HAS_SUBTILES:
 * @info: a #GstVideoFormatInfo
 *
 * This macro checks if %GST_VIDEO_FORMAT_FLAG_SUBTILES is set. When this
 * flag is set, it means that the tile sizes must be scaled as per the
 * subsampling.
 *
 * Returns: %TRUE if the format uses subsampled tile sizes.
 *
 * Since: 1.22
 */
#define GST_VIDEO_FORMAT_INFO_HAS_SUBTILES(info) (((info)->flags & GST_VIDEO_FORMAT_FLAG_SUBTILES) != 0)

#define GST_VIDEO_FORMAT_INFO_BITS(info)         ((info)->bits)
#define GST_VIDEO_FORMAT_INFO_N_COMPONENTS(info) ((info)->n_components)
#define GST_VIDEO_FORMAT_INFO_SHIFT(info,c)      ((info)->shift[c])
#define GST_VIDEO_FORMAT_INFO_DEPTH(info,c)      ((info)->depth[c])
/**
 * GST_VIDEO_FORMAT_INFO_PSTRIDE:
 * @info: a #GstVideoFormatInfo
 * @c: the component index
 *
 * pixel stride for the given component. This is the amount of bytes to the
 * pixel immediately to the right, so basically bytes from one pixel to the
 * next. When bits < 8, the stride is expressed in bits.
 *
 * Examples: for 24-bit RGB, the pixel stride would be 3 bytes, while it
 * would be 4 bytes for RGBx or ARGB, and 8 bytes for ARGB64 or AYUV64.
 * For planar formats such as I420 the pixel stride is usually 1. For
 * YUY2 it would be 2 bytes.
 */
#define GST_VIDEO_FORMAT_INFO_PSTRIDE(info,c)    ((info)->pixel_stride[c])
/**
 * GST_VIDEO_FORMAT_INFO_N_PLANES:
 * @info: a #GstVideoFormatInfo
 *
 * Number of planes. This is the number of planes the pixel layout is
 * organized in in memory. The number of planes can be less than the
 * number of components (e.g. Y,U,V,A or R, G, B, A) when multiple
 * components are packed into one plane.
 *
 * Examples: RGB/RGBx/RGBA: 1 plane, 3/3/4 components;
 * I420: 3 planes, 3 components; NV21/NV12: 2 planes, 3 components.
 */
#define GST_VIDEO_FORMAT_INFO_N_PLANES(info)     ((info)->n_planes)
/**
 * GST_VIDEO_FORMAT_INFO_PLANE:
 * @info: a #GstVideoFormatInfo
 * @c: the component index
 *
 * Plane number where the given component can be found. A plane may
 * contain data for multiple components.
 */
#define GST_VIDEO_FORMAT_INFO_PLANE(info,c)      ((info)->plane[c])
#define GST_VIDEO_FORMAT_INFO_POFFSET(info,c)    ((info)->poffset[c])
#define GST_VIDEO_FORMAT_INFO_W_SUB(info,c)      ((info)->w_sub[c])
#define GST_VIDEO_FORMAT_INFO_H_SUB(info,c)      ((info)->h_sub[c])

/* rounds up */
#define GST_VIDEO_SUB_SCALE(scale,val)   (-((-((gint)(val)))>>(scale)))

#define GST_VIDEO_FORMAT_INFO_SCALE_WIDTH(info,c,w)  GST_VIDEO_SUB_SCALE ((info)->w_sub[c],(w))
#define GST_VIDEO_FORMAT_INFO_SCALE_HEIGHT(info,c,h) GST_VIDEO_SUB_SCALE ((info)->h_sub[c],(h))

#define GST_VIDEO_FORMAT_INFO_DATA(info,planes,comp) \
  (((guint8*)(planes)[(info)->plane[comp]]) + (info)->poffset[comp])
/**
 * GST_VIDEO_FORMAT_INFO_STRIDE:
 * @info: a #GstVideoFormatInfo
 * @strides: an array of strides
 * @comp: the component index
 *
 * Row stride in bytes, that is number of bytes from the first pixel component
 * of a row to the first pixel component in the next row. This might include
 * some row padding (memory not actually used for anything, to make sure the
 * beginning of the next row is aligned in a particular way).
 */
#define GST_VIDEO_FORMAT_INFO_STRIDE(info,strides,comp) ((strides)[(info)->plane[comp]])
#define GST_VIDEO_FORMAT_INFO_OFFSET(info,offsets,comp) \
  (((offsets)[(info)->plane[comp]]) + (info)->poffset[comp])

#define GST_VIDEO_FORMAT_INFO_TILE_MODE(info) ((info)->tile_mode)
#define GST_VIDEO_FORMAT_INFO_TILE_WS(info) ((info)->tile_ws)
#define GST_VIDEO_FORMAT_INFO_TILE_HS(info) ((info)->tile_hs)

/**
 * GST_VIDEO_FORMAT_INFO_TILE_SIZE:
 * @info: a #GstVideoFormatInfo
 * @plane: the plane index
 *
 * Provides the size in bytes of a tile in the specified @plane. This replaces
 * the width and height shift, which was limited to power of two dimensions.
 *
 * Since: 1.22
 */
#define GST_VIDEO_FORMAT_INFO_TILE_SIZE(info,plane) ((info)->tile_info[plane].size)

/**
 * GST_VIDEO_FORMAT_INFO_TILE_WIDTH:
 * @info: a #GstVideoFormatInfo
 * @plane: the plane index
 *
 * See #GstVideoTileInfo.width.
 *
 * Return the width of one tile in pixels, zero if its not an integer.
 *
 * Since: 1.22
 */
#define GST_VIDEO_FORMAT_INFO_TILE_WIDTH(info,plane) ((info)->tile_info[plane].width)

/**
 * GST_VIDEO_FORMAT_INFO_TILE_HEIGHT:
 * @info: a #GstVideoFormatInfo
 * @plane: the plane index
 *
 * See #GstVideoTileInfo.height.
 *
 * Returns the tile height.
 *
 * Since: 1.22
 */
#define GST_VIDEO_FORMAT_INFO_TILE_HEIGHT(info,plane) ((info)->tile_info[plane].height)

/**
 * GST_VIDEO_FORMAT_INFO_TILE_STRIDE:
 * @info: a #GstVideoFormatInfo
 * @plane: the plane index
 *
 * See #GstVideoTileInfo.stride.
 *
 * Returns the stride of one tile, regardless of the internal details of the
 * tile (could be a complex system with subtile) the tiles size should alway
 * match the tile width multiplied by the tile stride.
 *
 * Since: 1.22
 */
#define GST_VIDEO_FORMAT_INFO_TILE_STRIDE(info,plane) ((info)->tile_info[plane].stride)


GST_VIDEO_API
void gst_video_format_info_component                  (const GstVideoFormatInfo *info, guint plane, gint components[GST_VIDEO_MAX_COMPONENTS]);

GST_VIDEO_API
gint gst_video_format_info_extrapolate_stride        (const GstVideoFormatInfo * finfo,
                                                      gint plane, gint stride);

/* format properties */

GST_VIDEO_API
GstVideoFormat gst_video_format_from_masks           (gint depth, gint bpp, gint endianness,
                                                      guint red_mask, guint green_mask,
                                                      guint blue_mask, guint alpha_mask) G_GNUC_CONST;

GST_VIDEO_API
GstVideoFormat gst_video_format_from_fourcc          (guint32 fourcc) G_GNUC_CONST;

GST_VIDEO_API
GstVideoFormat gst_video_format_from_string          (const gchar *format) G_GNUC_CONST;

GST_VIDEO_API
guint32        gst_video_format_to_fourcc            (GstVideoFormat format) G_GNUC_CONST;

GST_VIDEO_API
const gchar *  gst_video_format_to_string            (GstVideoFormat format) G_GNUC_CONST;

GST_VIDEO_API
const GstVideoFormatInfo *
               gst_video_format_get_info             (GstVideoFormat format) G_GNUC_CONST;

GST_VIDEO_API
gconstpointer  gst_video_format_get_palette          (GstVideoFormat format, gsize *size);

#define GST_VIDEO_SIZE_RANGE "(int) [ 1, max ]"
#define GST_VIDEO_FPS_RANGE "(fraction) [ 0, max ]"

#if G_BYTE_ORDER == G_LITTLE_ENDIAN
# define GST_VIDEO_NE(s) G_STRINGIFY(s)"_LE"
# define GST_VIDEO_OE(s) G_STRINGIFY(s)"_BE"
#else
# define GST_VIDEO_NE(s) G_STRINGIFY(s)"_BE"
# define GST_VIDEO_OE(s) G_STRINGIFY(s)"_LE"
#endif

/**
 * GST_VIDEO_FORMATS_ALL:
 *
 * List of all video formats, for use in template caps strings.
 *
 * Formats are sorted by decreasing "quality", using these criteria by priority:
 *   - number of components
 *   - depth
 *   - subsampling factor of the width
 *   - subsampling factor of the height
 *   - number of planes
 *   - native endianness preferred
 *   - pixel stride
 *   - poffset
 *   - prefer non-complex formats
 *   - prefer YUV formats over RGB ones
 *   - prefer I420 over YV12
 *   - format name
 */
#if G_BYTE_ORDER == G_BIG_ENDIAN
#define GST_VIDEO_FORMATS_ALL "{ ABGR64_BE, BGRA64_BE, AYUV64, ARGB64_BE, ARGB64, " \
    "RGBA64_BE, ABGR64_LE, BGRA64_LE, ARGB64_LE, RGBA64_LE, GBRA_12BE, GBRA_12LE, Y412_BE, " \
    "Y412_LE, A444_10BE, GBRA_10BE, A444_10LE, GBRA_10LE, A422_10BE, A422_10LE, " \
    "A420_10BE, A420_10LE, Y410, RGB10A2_LE, BGR10A2_LE, GBRA, ABGR, VUYA, BGRA, " \
    "AYUV, ARGB, RGBA, A420, AV12, Y444_16BE, Y444_16LE, v216, P016_BE, P016_LE, Y444_12BE, " \
    "GBR_12BE, Y444_12LE, GBR_12LE, I422_12BE, I422_12LE, Y212_BE, Y212_LE, I420_12BE, " \
    "I420_12LE, P012_BE, P012_LE, Y444_10BE, GBR_10BE, Y444_10LE, GBR_10LE, r210, " \
    "I422_10BE, I422_10LE, NV16_10LE32, Y210, v210, UYVP, I420_10BE, I420_10LE, " \
    "P010_10BE, P010_10LE, NV12_10LE32, NV12_10LE40, NV12_10BE_8L128, Y444, RGBP, GBR, BGRP, NV24, xBGR, BGRx, " \
    "xRGB, RGBx, BGR, IYU2, v308, RGB, Y42B, NV61, NV16, VYUY, UYVY, YVYU, YUY2, I420, " \
    "YV12, NV21, NV12, NV12_8L128, NV12_64Z32, NV12_4L4, NV12_32L32, NV12_16L32S, Y41B, IYU1, YVU9, YUV9, RGB16, " \
    "BGR16, RGB15, BGR15, RGB8P, GRAY16_BE, GRAY16_LE, GRAY10_LE32, GRAY8 }"
#elif G_BYTE_ORDER == G_LITTLE_ENDIAN
#define GST_VIDEO_FORMATS_ALL "{ ABGR64_LE, BGRA64_LE, AYUV64, ARGB64_LE, ARGB64, " \
    "RGBA64_LE, ABGR64_BE, BGRA64_BE, ARGB64_BE, RGBA64_BE, GBRA_12LE, GBRA_12BE, Y412_LE, " \
    "Y412_BE, A444_10LE, GBRA_10LE, A444_10BE, GBRA_10BE, A422_10LE, A422_10BE, " \
    "A420_10LE, A420_10BE, RGB10A2_LE, BGR10A2_LE, Y410, GBRA, ABGR, VUYA, BGRA, " \
    "AYUV, ARGB, RGBA, A420, AV12, Y444_16LE, Y444_16BE, v216, P016_LE, P016_BE, Y444_12LE, " \
    "GBR_12LE, Y444_12BE, GBR_12BE, I422_12LE, I422_12BE, Y212_LE, Y212_BE, I420_12LE, " \
    "I420_12BE, P012_LE, P012_BE, Y444_10LE, GBR_10LE, Y444_10BE, GBR_10BE, r210, " \
    "I422_10LE, I422_10BE, NV16_10LE32, Y210, v210, UYVP, I420_10LE, I420_10BE, " \
    "P010_10LE, NV12_10LE32, NV12_10LE40, P010_10BE, NV12_10BE_8L128, Y444, RGBP, GBR, BGRP, NV24, xBGR, BGRx, " \
    "xRGB, RGBx, BGR, IYU2, v308, RGB, Y42B, NV61, NV16, VYUY, UYVY, YVYU, YUY2, I420, " \
    "YV12, NV21, NV12, NV12_8L128, NV12_64Z32, NV12_4L4, NV12_32L32, NV12_16L32S, Y41B, IYU1, YVU9, YUV9, RGB16, " \
    "BGR16, RGB15, BGR15, RGB8P, GRAY16_LE, GRAY16_BE, GRAY10_LE32, GRAY8 }"
#endif

GST_VIDEO_API
const GstVideoFormat * gst_video_formats_raw (guint * len);

/**
 * GST_VIDEO_CAPS_MAKE:
 * @format: string format that describes the pixel layout, as string
 *     (e.g. "I420", "RGB", "YV12", "YUY2", "AYUV", etc.)
 *
 * Generic caps string for video, for use in pad templates.
 */
#define GST_VIDEO_CAPS_MAKE(format)                                     \
    "video/x-raw, "                                                     \
    "format = (string) " format ", "                                    \
    "width = " GST_VIDEO_SIZE_RANGE ", "                                \
    "height = " GST_VIDEO_SIZE_RANGE ", "                               \
    "framerate = " GST_VIDEO_FPS_RANGE

/**
 * GST_VIDEO_CAPS_MAKE_WITH_FEATURES:
 * @format: string format that describes the pixel layout, as string
 *     (e.g. "I420", "RGB", "YV12", "YUY2", "AYUV", etc.)
 * @features: Requires caps features as a string, e.g.
 *     "memory:SystemMemory".
 *
 * Generic caps string for video, for use in pad templates.
 *
 * Since: 1.2
 */
#define GST_VIDEO_CAPS_MAKE_WITH_FEATURES(features,format)              \
    "video/x-raw(" features "), "                                       \
    "format = (string) " format ", "                                    \
    "width = " GST_VIDEO_SIZE_RANGE ", "                                \
    "height = " GST_VIDEO_SIZE_RANGE ", "                               \
    "framerate = " GST_VIDEO_FPS_RANGE

GST_VIDEO_API
GstCaps * gst_video_make_raw_caps (const GstVideoFormat formats[], guint len);

GST_VIDEO_API
GstCaps * gst_video_make_raw_caps_with_features (const GstVideoFormat formats[], guint len,
                                                 GstCapsFeatures * features);


G_END_DECLS

#endif /* __GST_VIDEO_FORMAT_H__ */
