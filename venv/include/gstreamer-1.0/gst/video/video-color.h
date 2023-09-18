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

#ifndef __GST_VIDEO_COLOR_H__
#define __GST_VIDEO_COLOR_H__

#include <gst/gst.h>

#include <gst/video/video-format.h>

G_BEGIN_DECLS

/**
 * GstVideoColorRange:
 * @GST_VIDEO_COLOR_RANGE_UNKNOWN: unknown range
 * @GST_VIDEO_COLOR_RANGE_0_255: [0..255] for 8 bit components
 * @GST_VIDEO_COLOR_RANGE_16_235: [16..235] for 8 bit components. Chroma has
 *                 [16..240] range.
 *
 * Possible color range values. These constants are defined for 8 bit color
 * values and can be scaled for other bit depths.
 */
typedef enum {
  GST_VIDEO_COLOR_RANGE_UNKNOWN = 0,
  GST_VIDEO_COLOR_RANGE_0_255,
  GST_VIDEO_COLOR_RANGE_16_235
} GstVideoColorRange;

/**
 * GstVideoColorMatrix:
 * @GST_VIDEO_COLOR_MATRIX_UNKNOWN: unknown matrix
 * @GST_VIDEO_COLOR_MATRIX_RGB: identity matrix. Order of coefficients is
 * actually GBR, also IEC 61966-2-1 (sRGB)
 * @GST_VIDEO_COLOR_MATRIX_FCC: FCC Title 47 Code of Federal Regulations 73.682 (a)(20)
 * @GST_VIDEO_COLOR_MATRIX_BT709: ITU-R BT.709 color matrix, also ITU-R BT1361
 * / IEC 61966-2-4 xvYCC709 / SMPTE RP177 Annex B
 * @GST_VIDEO_COLOR_MATRIX_BT601: ITU-R BT.601 color matrix, also SMPTE170M / ITU-R BT1358 525 / ITU-R BT1700 NTSC
 * @GST_VIDEO_COLOR_MATRIX_SMPTE240M: SMPTE 240M color matrix
 * @GST_VIDEO_COLOR_MATRIX_BT2020: ITU-R BT.2020 color matrix. Since: 1.6
 *
 * The color matrix is used to convert between Y'PbPr and
 * non-linear RGB (R'G'B')
 */
typedef enum {
  GST_VIDEO_COLOR_MATRIX_UNKNOWN = 0,
  GST_VIDEO_COLOR_MATRIX_RGB,
  GST_VIDEO_COLOR_MATRIX_FCC,
  GST_VIDEO_COLOR_MATRIX_BT709,
  GST_VIDEO_COLOR_MATRIX_BT601,
  GST_VIDEO_COLOR_MATRIX_SMPTE240M,
  GST_VIDEO_COLOR_MATRIX_BT2020
} GstVideoColorMatrix;

GST_VIDEO_API
gboolean gst_video_color_matrix_get_Kr_Kb (GstVideoColorMatrix matrix, gdouble * Kr, gdouble * Kb);

/**
 * GstVideoTransferFunction:
 * @GST_VIDEO_TRANSFER_UNKNOWN: unknown transfer function
 * @GST_VIDEO_TRANSFER_GAMMA10: linear RGB, gamma 1.0 curve
 * @GST_VIDEO_TRANSFER_GAMMA18: Gamma 1.8 curve
 * @GST_VIDEO_TRANSFER_GAMMA20: Gamma 2.0 curve
 * @GST_VIDEO_TRANSFER_GAMMA22: Gamma 2.2 curve
 * @GST_VIDEO_TRANSFER_BT709: Gamma 2.2 curve with a linear segment in the lower
 *                           range, also ITU-R BT470M / ITU-R BT1700 625 PAL &
 *                           SECAM / ITU-R BT1361
 * @GST_VIDEO_TRANSFER_SMPTE240M: Gamma 2.2 curve with a linear segment in the
 *                               lower range
 * @GST_VIDEO_TRANSFER_SRGB: Gamma 2.4 curve with a linear segment in the lower
 *                          range. IEC 61966-2-1 (sRGB or sYCC)
 * @GST_VIDEO_TRANSFER_GAMMA28: Gamma 2.8 curve, also ITU-R BT470BG
 * @GST_VIDEO_TRANSFER_LOG100: Logarithmic transfer characteristic
 *                             100:1 range
 * @GST_VIDEO_TRANSFER_LOG316: Logarithmic transfer characteristic
 *                             316.22777:1 range (100 * sqrt(10) : 1)
 * @GST_VIDEO_TRANSFER_BT2020_12: Gamma 2.2 curve with a linear segment in the lower
 *                                range. Used for BT.2020 with 12 bits per
 *                                component. Since: 1.6
 * @GST_VIDEO_TRANSFER_ADOBERGB: Gamma 2.19921875. Since: 1.8
 * @GST_VIDEO_TRANSFER_BT2020_10: Rec. ITU-R BT.2020-2 with 10 bits per component.
 *                                (functionally the same as the values
 *                                GST_VIDEO_TRANSFER_BT709 and GST_VIDEO_TRANSFER_BT601).
 *                                Since: 1.18
 * @GST_VIDEO_TRANSFER_SMPTE2084: SMPTE ST 2084 for 10, 12, 14, and 16-bit systems.
 *                                Known as perceptual quantization (PQ)
 *                                Since: 1.18
 * @GST_VIDEO_TRANSFER_ARIB_STD_B67: Association of Radio Industries and Businesses (ARIB)
 *                                   STD-B67 and Rec. ITU-R BT.2100-1 hybrid loggamma (HLG) system
 *                                   Since: 1.18
 * @GST_VIDEO_TRANSFER_BT601: also known as SMPTE170M / ITU-R BT1358 525 or 625 / ITU-R BT1700 NTSC
 *                            Functionally the same as the values
 *                            GST_VIDEO_TRANSFER_BT709, and GST_VIDEO_TRANSFER_BT2020_10.
 *                            Since: 1.18
 *
 * The video transfer function defines the formula for converting between
 * non-linear RGB (R'G'B') and linear RGB
 */
typedef enum {
  GST_VIDEO_TRANSFER_UNKNOWN = 0,
  GST_VIDEO_TRANSFER_GAMMA10,
  GST_VIDEO_TRANSFER_GAMMA18,
  GST_VIDEO_TRANSFER_GAMMA20,
  GST_VIDEO_TRANSFER_GAMMA22,
  GST_VIDEO_TRANSFER_BT709,
  GST_VIDEO_TRANSFER_SMPTE240M,
  GST_VIDEO_TRANSFER_SRGB,
  GST_VIDEO_TRANSFER_GAMMA28,
  GST_VIDEO_TRANSFER_LOG100,
  GST_VIDEO_TRANSFER_LOG316,
  GST_VIDEO_TRANSFER_BT2020_12,
  GST_VIDEO_TRANSFER_ADOBERGB,
  GST_VIDEO_TRANSFER_BT2020_10,
  GST_VIDEO_TRANSFER_SMPTE2084,
  GST_VIDEO_TRANSFER_ARIB_STD_B67,
  /**
   * GST_VIDEO_TRANSFER_BT601:
   *
   * also known as SMPTE170M / ITU-R BT1358 525 or 625 / ITU-R BT1700 NTSC
   *
   * Since: 1.18
   */
  GST_VIDEO_TRANSFER_BT601
} GstVideoTransferFunction;

GST_VIDEO_DEPRECATED_FOR(gst_video_transfer_function_encode)
gdouble      gst_video_color_transfer_encode    (GstVideoTransferFunction func, gdouble val);
GST_VIDEO_API
gdouble      gst_video_transfer_function_encode (GstVideoTransferFunction func, gdouble val);

GST_VIDEO_DEPRECATED_FOR(gst_video_transfer_function_decode)
gdouble      gst_video_color_transfer_decode    (GstVideoTransferFunction func, gdouble val);
GST_VIDEO_API
gdouble      gst_video_transfer_function_decode (GstVideoTransferFunction func, gdouble val);

/**
 * GstVideoColorPrimaries:
 * @GST_VIDEO_COLOR_PRIMARIES_UNKNOWN: unknown color primaries
 * @GST_VIDEO_COLOR_PRIMARIES_BT709: BT709 primaries, also ITU-R BT1361 / IEC
 * 61966-2-4 / SMPTE RP177 Annex B
 * @GST_VIDEO_COLOR_PRIMARIES_BT470M: BT470M primaries, also FCC Title 47 Code
 * of Federal Regulations 73.682 (a)(20)
 * @GST_VIDEO_COLOR_PRIMARIES_BT470BG: BT470BG primaries, also ITU-R BT601-6
 * 625 / ITU-R BT1358 625 / ITU-R BT1700 625 PAL & SECAM
 * @GST_VIDEO_COLOR_PRIMARIES_SMPTE170M: SMPTE170M primaries, also ITU-R
 * BT601-6 525 / ITU-R BT1358 525 / ITU-R BT1700 NTSC
 * @GST_VIDEO_COLOR_PRIMARIES_SMPTE240M: SMPTE240M primaries
 * @GST_VIDEO_COLOR_PRIMARIES_FILM: Generic film (colour filters using
 * Illuminant C)
 * @GST_VIDEO_COLOR_PRIMARIES_BT2020: ITU-R BT2020 primaries. Since: 1.6
 * @GST_VIDEO_COLOR_PRIMARIES_ADOBERGB: Adobe RGB primaries. Since: 1.8
 * @GST_VIDEO_COLOR_PRIMARIES_SMPTEST428: SMPTE ST 428 primaries (CIE 1931
 * XYZ). Since: 1.16
 * @GST_VIDEO_COLOR_PRIMARIES_SMPTERP431: SMPTE RP 431 primaries (ST 431-2
 * (2011) / DCI P3). Since: 1.16
 * @GST_VIDEO_COLOR_PRIMARIES_SMPTEEG432: SMPTE EG 432 primaries (ST 432-1
 * (2010) / P3 D65). Since: 1.16
 * @GST_VIDEO_COLOR_PRIMARIES_EBU3213: EBU 3213 primaries (JEDEC P22
 * phosphors). Since: 1.16
 *
 * The color primaries define the how to transform linear RGB values to and from
 * the CIE XYZ colorspace.
 */
typedef enum {
  GST_VIDEO_COLOR_PRIMARIES_UNKNOWN = 0,
  GST_VIDEO_COLOR_PRIMARIES_BT709,
  GST_VIDEO_COLOR_PRIMARIES_BT470M,
  GST_VIDEO_COLOR_PRIMARIES_BT470BG,
  GST_VIDEO_COLOR_PRIMARIES_SMPTE170M,
  GST_VIDEO_COLOR_PRIMARIES_SMPTE240M,
  GST_VIDEO_COLOR_PRIMARIES_FILM,
  GST_VIDEO_COLOR_PRIMARIES_BT2020,
  GST_VIDEO_COLOR_PRIMARIES_ADOBERGB,
  GST_VIDEO_COLOR_PRIMARIES_SMPTEST428,
  GST_VIDEO_COLOR_PRIMARIES_SMPTERP431,
  GST_VIDEO_COLOR_PRIMARIES_SMPTEEG432,
  GST_VIDEO_COLOR_PRIMARIES_EBU3213,
} GstVideoColorPrimaries;

/**
 * GstVideoColorPrimariesInfo:
 * @primaries: a #GstVideoColorPrimaries
 * @Wx: reference white x coordinate
 * @Wy: reference white y coordinate
 * @Rx: red x coordinate
 * @Ry: red y coordinate
 * @Gx: green x coordinate
 * @Gy: green y coordinate
 * @Bx: blue x coordinate
 * @By: blue y coordinate
 *
 * Structure describing the chromaticity coordinates of an RGB system. These
 * values can be used to construct a matrix to transform RGB to and from the
 * XYZ colorspace.
 *
 * Since: 1.6
 */
typedef struct {
  GstVideoColorPrimaries primaries;
  gdouble Wx, Wy;
  gdouble Rx, Ry;
  gdouble Gx, Gy;
  gdouble Bx, By;
} GstVideoColorPrimariesInfo;

GST_VIDEO_API
const GstVideoColorPrimariesInfo *
                gst_video_color_primaries_get_info     (GstVideoColorPrimaries primaries);

GST_VIDEO_API
gboolean gst_video_color_primaries_is_equivalent       (GstVideoColorPrimaries primaries,
                                                        GstVideoColorPrimaries other);

/**
 * GstVideoColorimetry:
 * @range: the color range. This is the valid range for the samples.
 *         It is used to convert the samples to Y'PbPr values.
 * @matrix: the color matrix. Used to convert between Y'PbPr and
 *          non-linear RGB (R'G'B')
 * @transfer: the transfer function. used to convert between R'G'B' and RGB
 * @primaries: color primaries. used to convert between R'G'B' and CIE XYZ
 *
 * Structure describing the color info.
 */
typedef struct {
  GstVideoColorRange        range;
  GstVideoColorMatrix       matrix;
  GstVideoTransferFunction  transfer;
  GstVideoColorPrimaries    primaries;
} GstVideoColorimetry;

/* predefined colorimetry */
#define GST_VIDEO_COLORIMETRY_BT601       "bt601"
#define GST_VIDEO_COLORIMETRY_BT709       "bt709"
#define GST_VIDEO_COLORIMETRY_SMPTE240M   "smpte240m"
#define GST_VIDEO_COLORIMETRY_SRGB        "sRGB"
#define GST_VIDEO_COLORIMETRY_BT2020      "bt2020"
#define GST_VIDEO_COLORIMETRY_BT2020_10   "bt2020-10"
#define GST_VIDEO_COLORIMETRY_BT2100_PQ   "bt2100-pq"
#define GST_VIDEO_COLORIMETRY_BT2100_HLG  "bt2100-hlg"

GST_VIDEO_API
gboolean     gst_video_colorimetry_matches     (const GstVideoColorimetry *cinfo, const gchar *color);

GST_VIDEO_API
gboolean     gst_video_colorimetry_from_string (GstVideoColorimetry *cinfo, const gchar *color);

GST_VIDEO_API
gchar *      gst_video_colorimetry_to_string   (const GstVideoColorimetry *cinfo);

GST_VIDEO_API
gboolean     gst_video_colorimetry_is_equal    (const GstVideoColorimetry *cinfo, const GstVideoColorimetry *other);

GST_VIDEO_API
gboolean     gst_video_colorimetry_is_equivalent (const GstVideoColorimetry *cinfo,
                                                  guint bitdepth,
                                                  const GstVideoColorimetry *other,
                                                  guint other_bitdepth);

/* compute offset and scale */

GST_VIDEO_API
void         gst_video_color_range_offsets     (GstVideoColorRange range,
                                                const GstVideoFormatInfo *info,
                                                gint offset[GST_VIDEO_MAX_COMPONENTS],
                                                gint scale[GST_VIDEO_MAX_COMPONENTS]);

/* conversion between GStreamer color{matrix,transfer,primaries} enum and
 * values defined by ISO/IEC 23001-8 and ITU-T H.273 specification.
 * Also H264 and H265 specifications follow the color{matrix,transfer,primaries}
 * values */

GST_VIDEO_API
guint                     gst_video_color_matrix_to_iso      (GstVideoColorMatrix matrix);

GST_VIDEO_API
guint                     gst_video_transfer_function_to_iso    (GstVideoTransferFunction func);

GST_VIDEO_API
guint                     gst_video_color_primaries_to_iso   (GstVideoColorPrimaries primaries);

GST_VIDEO_API
GstVideoColorMatrix       gst_video_color_matrix_from_iso    (guint value);

GST_VIDEO_API
GstVideoTransferFunction  gst_video_transfer_function_from_iso  (guint value);

GST_VIDEO_API
GstVideoColorPrimaries    gst_video_color_primaries_from_iso (guint value);

GST_VIDEO_API
gboolean                  gst_video_transfer_function_is_equivalent (GstVideoTransferFunction from_func,
                                                                    guint from_bpp,
                                                                    GstVideoTransferFunction to_func,
                                                                    guint to_bpp);

G_END_DECLS

#endif /* __GST_VIDEO_COLOR_H__ */
