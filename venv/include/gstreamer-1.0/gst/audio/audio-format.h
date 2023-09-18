/* GStreamer
 * Copyright (C) <1999> Erik Walthinsen <omega@cse.ogi.edu>
 * Library       <2001> Thomas Vander Stichele <thomas@apestaart.org>
 *               <2011> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef __GST_AUDIO_FORMAT_H__
#define __GST_AUDIO_FORMAT_H__

G_BEGIN_DECLS

#if G_BYTE_ORDER == G_BIG_ENDIAN
#define _GST_AUDIO_FORMAT_NE(fmt) GST_AUDIO_FORMAT_ ## fmt ## BE
#define _GST_AUDIO_FORMAT_OE(fmt) GST_AUDIO_FORMAT_ ## fmt ## LE
#elif G_BYTE_ORDER == G_LITTLE_ENDIAN
#define _GST_AUDIO_FORMAT_NE(fmt) GST_AUDIO_FORMAT_ ## fmt ## LE
#define _GST_AUDIO_FORMAT_OE(fmt) GST_AUDIO_FORMAT_ ## fmt ## BE
#endif

/**
 * GstAudioFormat:
 * @GST_AUDIO_FORMAT_UNKNOWN: unknown or unset audio format
 * @GST_AUDIO_FORMAT_ENCODED: encoded audio format
 * @GST_AUDIO_FORMAT_S8: 8 bits in 8 bits, signed
 * @GST_AUDIO_FORMAT_U8: 8 bits in 8 bits, unsigned
 * @GST_AUDIO_FORMAT_S16LE: 16 bits in 16 bits, signed, little endian
 * @GST_AUDIO_FORMAT_S16BE: 16 bits in 16 bits, signed, big endian
 * @GST_AUDIO_FORMAT_U16LE: 16 bits in 16 bits, unsigned, little endian
 * @GST_AUDIO_FORMAT_U16BE: 16 bits in 16 bits, unsigned, big endian
 * @GST_AUDIO_FORMAT_S24_32LE: 24 bits in 32 bits, signed, little endian
 * @GST_AUDIO_FORMAT_S24_32BE: 24 bits in 32 bits, signed, big endian
 * @GST_AUDIO_FORMAT_U24_32LE: 24 bits in 32 bits, unsigned, little endian
 * @GST_AUDIO_FORMAT_U24_32BE: 24 bits in 32 bits, unsigned, big endian
 * @GST_AUDIO_FORMAT_S32LE: 32 bits in 32 bits, signed, little endian
 * @GST_AUDIO_FORMAT_S32BE: 32 bits in 32 bits, signed, big endian
 * @GST_AUDIO_FORMAT_U32LE: 32 bits in 32 bits, unsigned, little endian
 * @GST_AUDIO_FORMAT_U32BE: 32 bits in 32 bits, unsigned, big endian
 * @GST_AUDIO_FORMAT_S24LE: 24 bits in 24 bits, signed, little endian
 * @GST_AUDIO_FORMAT_S24BE: 24 bits in 24 bits, signed, big endian
 * @GST_AUDIO_FORMAT_U24LE: 24 bits in 24 bits, unsigned, little endian
 * @GST_AUDIO_FORMAT_U24BE: 24 bits in 24 bits, unsigned, big endian
 * @GST_AUDIO_FORMAT_S20LE: 20 bits in 24 bits, signed, little endian
 * @GST_AUDIO_FORMAT_S20BE: 20 bits in 24 bits, signed, big endian
 * @GST_AUDIO_FORMAT_U20LE: 20 bits in 24 bits, unsigned, little endian
 * @GST_AUDIO_FORMAT_U20BE: 20 bits in 24 bits, unsigned, big endian
 * @GST_AUDIO_FORMAT_S18LE: 18 bits in 24 bits, signed, little endian
 * @GST_AUDIO_FORMAT_S18BE: 18 bits in 24 bits, signed, big endian
 * @GST_AUDIO_FORMAT_U18LE: 18 bits in 24 bits, unsigned, little endian
 * @GST_AUDIO_FORMAT_U18BE: 18 bits in 24 bits, unsigned, big endian
 * @GST_AUDIO_FORMAT_F32LE: 32-bit floating point samples, little endian
 * @GST_AUDIO_FORMAT_F32BE: 32-bit floating point samples, big endian
 * @GST_AUDIO_FORMAT_F64LE: 64-bit floating point samples, little endian
 * @GST_AUDIO_FORMAT_F64BE: 64-bit floating point samples, big endian
 * @GST_AUDIO_FORMAT_S16: 16 bits in 16 bits, signed, native endianness
 * @GST_AUDIO_FORMAT_U16: 16 bits in 16 bits, unsigned, native endianness
 * @GST_AUDIO_FORMAT_S24_32: 24 bits in 32 bits, signed, native endianness
 * @GST_AUDIO_FORMAT_U24_32: 24 bits in 32 bits, unsigned, native endianness
 * @GST_AUDIO_FORMAT_S32: 32 bits in 32 bits, signed, native endianness
 * @GST_AUDIO_FORMAT_U32: 32 bits in 32 bits, unsigned, native endianness
 * @GST_AUDIO_FORMAT_S24: 24 bits in 24 bits, signed, native endianness
 * @GST_AUDIO_FORMAT_U24: 24 bits in 24 bits, unsigned, native endianness
 * @GST_AUDIO_FORMAT_S20: 20 bits in 24 bits, signed, native endianness
 * @GST_AUDIO_FORMAT_U20: 20 bits in 24 bits, unsigned, native endianness
 * @GST_AUDIO_FORMAT_S18: 18 bits in 24 bits, signed, native endianness
 * @GST_AUDIO_FORMAT_U18: 18 bits in 24 bits, unsigned, native endianness
 * @GST_AUDIO_FORMAT_F32: 32-bit floating point samples, native endianness
 * @GST_AUDIO_FORMAT_F64: 64-bit floating point samples, native endianness
 *
 * Enum value describing the most common audio formats.
 */
typedef enum {
  GST_AUDIO_FORMAT_UNKNOWN,
  GST_AUDIO_FORMAT_ENCODED,
  /* 8 bit */
  GST_AUDIO_FORMAT_S8,
  GST_AUDIO_FORMAT_U8,
  /* 16 bit */
  GST_AUDIO_FORMAT_S16LE,
  GST_AUDIO_FORMAT_S16BE,
  GST_AUDIO_FORMAT_U16LE,
  GST_AUDIO_FORMAT_U16BE,
  /* 24 bit in low 3 bytes of 32 bits*/
  GST_AUDIO_FORMAT_S24_32LE,
  GST_AUDIO_FORMAT_S24_32BE,
  GST_AUDIO_FORMAT_U24_32LE,
  GST_AUDIO_FORMAT_U24_32BE,
  /* 32 bit */
  GST_AUDIO_FORMAT_S32LE,
  GST_AUDIO_FORMAT_S32BE,
  GST_AUDIO_FORMAT_U32LE,
  GST_AUDIO_FORMAT_U32BE,
  /* 24 bit in 3 bytes*/
  GST_AUDIO_FORMAT_S24LE,
  GST_AUDIO_FORMAT_S24BE,
  GST_AUDIO_FORMAT_U24LE,
  GST_AUDIO_FORMAT_U24BE,
  /* 20 bit in 3 bytes*/
  GST_AUDIO_FORMAT_S20LE,
  GST_AUDIO_FORMAT_S20BE,
  GST_AUDIO_FORMAT_U20LE,
  GST_AUDIO_FORMAT_U20BE,
  /* 18 bit in 3 bytes*/
  GST_AUDIO_FORMAT_S18LE,
  GST_AUDIO_FORMAT_S18BE,
  GST_AUDIO_FORMAT_U18LE,
  GST_AUDIO_FORMAT_U18BE,
  /* float */
  GST_AUDIO_FORMAT_F32LE,
  GST_AUDIO_FORMAT_F32BE,
  GST_AUDIO_FORMAT_F64LE,
  GST_AUDIO_FORMAT_F64BE,
  /* native endianness equivalents */
  GST_AUDIO_FORMAT_S16 = _GST_AUDIO_FORMAT_NE(S16),
  GST_AUDIO_FORMAT_U16 = _GST_AUDIO_FORMAT_NE(U16),
  GST_AUDIO_FORMAT_S24_32 = _GST_AUDIO_FORMAT_NE(S24_32),
  GST_AUDIO_FORMAT_U24_32 = _GST_AUDIO_FORMAT_NE(U24_32),
  GST_AUDIO_FORMAT_S32 = _GST_AUDIO_FORMAT_NE(S32),
  GST_AUDIO_FORMAT_U32 = _GST_AUDIO_FORMAT_NE(U32),
  GST_AUDIO_FORMAT_S24 = _GST_AUDIO_FORMAT_NE(S24),
  GST_AUDIO_FORMAT_U24 = _GST_AUDIO_FORMAT_NE(U24),
  GST_AUDIO_FORMAT_S20 = _GST_AUDIO_FORMAT_NE(S20),
  GST_AUDIO_FORMAT_U20 = _GST_AUDIO_FORMAT_NE(U20),
  GST_AUDIO_FORMAT_S18 = _GST_AUDIO_FORMAT_NE(S18),
  GST_AUDIO_FORMAT_U18 = _GST_AUDIO_FORMAT_NE(U18),
  GST_AUDIO_FORMAT_F32 = _GST_AUDIO_FORMAT_NE(F32),
  GST_AUDIO_FORMAT_F64 = _GST_AUDIO_FORMAT_NE(F64)
} GstAudioFormat;


typedef struct _GstAudioFormatInfo GstAudioFormatInfo;

/**
 * GstAudioFormatFlags:
 * @GST_AUDIO_FORMAT_FLAG_INTEGER: integer samples
 * @GST_AUDIO_FORMAT_FLAG_FLOAT: float samples
 * @GST_AUDIO_FORMAT_FLAG_SIGNED: signed samples
 * @GST_AUDIO_FORMAT_FLAG_COMPLEX: complex layout
 * @GST_AUDIO_FORMAT_FLAG_UNPACK: the format can be used in
 * #GstAudioFormatUnpack and #GstAudioFormatPack functions
 *
 * The different audio flags that a format info can have.
 */
typedef enum
{
  GST_AUDIO_FORMAT_FLAG_INTEGER  = (1 << 0),
  GST_AUDIO_FORMAT_FLAG_FLOAT    = (1 << 1),
  GST_AUDIO_FORMAT_FLAG_SIGNED   = (1 << 2),
  GST_AUDIO_FORMAT_FLAG_COMPLEX  = (1 << 4),
  GST_AUDIO_FORMAT_FLAG_UNPACK   = (1 << 5)
} GstAudioFormatFlags;

/**
 * GstAudioPackFlags:
 * @GST_AUDIO_PACK_FLAG_NONE: No flag
 * @GST_AUDIO_PACK_FLAG_TRUNCATE_RANGE: When the source has a smaller depth
 *   than the target format, set the least significant bits of the target
 *   to 0. This is likely slightly faster but less accurate. When this flag
 *   is not specified, the most significant bits of the source are duplicated
 *   in the least significant bits of the destination.
 *
 * The different flags that can be used when packing and unpacking.
 */
typedef enum
{
  GST_AUDIO_PACK_FLAG_NONE             = 0,
  GST_AUDIO_PACK_FLAG_TRUNCATE_RANGE   = (1 << 0)
} GstAudioPackFlags;

/**
 * GstAudioFormatUnpack:
 * @info: a #GstAudioFormatInfo
 * @flags: #GstAudioPackFlags
 * @dest: (array) (element-type guint8): a destination array
 * @data: (array) (element-type guint8): pointer to the audio data
 * @length: the amount of samples to unpack.
 *
 * Unpacks @length samples from the given data of format @info.
 * The samples will be unpacked into @dest which each channel
 * interleaved. @dest should at least be big enough to hold @length *
 * channels * size(unpack_format) bytes.
 */
typedef void (*GstAudioFormatUnpack)         (const GstAudioFormatInfo *info,
                                              GstAudioPackFlags flags, gpointer dest,
                                              gconstpointer data, gint length);
/**
 * GstAudioFormatPack:
 * @info: a #GstAudioFormatInfo
 * @flags: #GstAudioPackFlags
 * @src: (array) (element-type guint8): a source array
 * @data: (array) (element-type guint8): pointer to the destination
 *   data
 * @length: the amount of samples to pack.
 *
 * Packs @length samples from @src to the data array in format @info.
 * The samples from source have each channel interleaved
 * and will be packed into @data.
 */
typedef void (*GstAudioFormatPack)           (const GstAudioFormatInfo *info,
                                              GstAudioPackFlags flags, gconstpointer src,
                                              gpointer data, gint length);

/**
 * GstAudioFormatInfo:
 * @format: #GstAudioFormat
 * @name: string representation of the format
 * @description: user readable description of the format
 * @flags: #GstAudioFormatFlags
 * @endianness: the endianness
 * @width: amount of bits used for one sample
 * @depth: amount of valid bits in @width
 * @silence: @width/8 bytes with 1 silent sample
 * @unpack_format: the format of the unpacked samples
 * @unpack_func: function to unpack samples
 * @pack_func: function to pack samples
 *
 * Information for an audio format.
 */
struct _GstAudioFormatInfo {
  /*< public >*/
  GstAudioFormat format;
  const gchar *name;
  const gchar *description;
  GstAudioFormatFlags flags;
  gint endianness;
  gint width;
  gint depth;
  guint8 silence[8];

  GstAudioFormat unpack_format;
  GstAudioFormatUnpack unpack_func;
  GstAudioFormatPack pack_func;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType gst_audio_format_info_get_type (void);

/**
 * GST_AUDIO_FORMAT_INFO_IS_VALID_RAW:
 *
 * Tests that the given #GstAudioFormatInfo represents a valid un-encoded
 * format.
 *
 * Since: 1.22
 */
#define GST_AUDIO_FORMAT_INFO_IS_VALID_RAW(info)                        \
  (info != NULL && (info)->format > GST_AUDIO_FORMAT_ENCODED &&         \
   (info)->width > 0 && (info)->depth > 0)

#define GST_AUDIO_FORMAT_INFO_FORMAT(info)           ((info)->format)
#define GST_AUDIO_FORMAT_INFO_NAME(info)             ((info)->name)
#define GST_AUDIO_FORMAT_INFO_FLAGS(info)            ((info)->flags)

#define GST_AUDIO_FORMAT_INFO_IS_INTEGER(info)       !!((info)->flags & GST_AUDIO_FORMAT_FLAG_INTEGER)
#define GST_AUDIO_FORMAT_INFO_IS_FLOAT(info)         !!((info)->flags & GST_AUDIO_FORMAT_FLAG_FLOAT)
#define GST_AUDIO_FORMAT_INFO_IS_SIGNED(info)        !!((info)->flags & GST_AUDIO_FORMAT_FLAG_SIGNED)

#define GST_AUDIO_FORMAT_INFO_ENDIANNESS(info)       ((info)->endianness)
#define GST_AUDIO_FORMAT_INFO_IS_LITTLE_ENDIAN(info) ((info)->endianness == G_LITTLE_ENDIAN)
#define GST_AUDIO_FORMAT_INFO_IS_BIG_ENDIAN(info)    ((info)->endianness == G_BIG_ENDIAN)
#define GST_AUDIO_FORMAT_INFO_WIDTH(info)            ((info)->width)
#define GST_AUDIO_FORMAT_INFO_DEPTH(info)            ((info)->depth)


GST_AUDIO_API
GstAudioFormat gst_audio_format_build_integer    (gboolean sign, gint endianness,
                                                  gint width, gint depth) G_GNUC_CONST;

GST_AUDIO_API
GstAudioFormat gst_audio_format_from_string      (const gchar *format) G_GNUC_CONST;

GST_AUDIO_API
const gchar *  gst_audio_format_to_string        (GstAudioFormat format) G_GNUC_CONST;

GST_AUDIO_API
const GstAudioFormatInfo *
               gst_audio_format_get_info         (GstAudioFormat format) G_GNUC_CONST;

GST_AUDIO_API
void           gst_audio_format_info_fill_silence (const GstAudioFormatInfo *info,
                                                   gpointer dest, gsize length);
GST_AUDIO_API G_DEPRECATED_FOR(gst_audio_format_info_fill_silence)
void           gst_audio_format_fill_silence      (const GstAudioFormatInfo *info,
                                                   gpointer dest, gsize length);

/**
 * GST_AUDIO_RATE_RANGE:
 *
 * Maximum range of allowed sample rates, for use in template caps strings.
 */
#define GST_AUDIO_RATE_RANGE "(int) [ 1, max ]"
/**
 * GST_AUDIO_CHANNELS_RANGE:
 *
 * Maximum range of allowed channels, for use in template caps strings.
 */
#define GST_AUDIO_CHANNELS_RANGE "(int) [ 1, max ]"

/**
 * GST_AUDIO_NE:
 * @s: format string without endianness marker
 *
 * Turns audio format string @s into the format string for native endianness.
 */
/**
 * GST_AUDIO_OE:
 * @s: format string without endianness marker
 *
 * Turns audio format string @s into the format string for other endianness.
 */
#if G_BYTE_ORDER == G_LITTLE_ENDIAN
# define GST_AUDIO_NE(s) G_STRINGIFY(s)"LE"
# define GST_AUDIO_OE(s) G_STRINGIFY(s)"BE"
#else
# define GST_AUDIO_NE(s) G_STRINGIFY(s)"BE"
# define GST_AUDIO_OE(s) G_STRINGIFY(s)"LE"
#endif

/**
 * GST_AUDIO_FORMATS_ALL:
 *
 * List of all audio formats, for use in template caps strings.
 *
 * Formats are sorted by decreasing "quality", using these criteria by priority:
 *   - depth
 *   - width
 *   - Float > Signed > Unsigned
 *   - native endianness preferred
 */
#if G_BYTE_ORDER == G_BIG_ENDIAN
#define GST_AUDIO_FORMATS_ALL "{ F64BE, F64LE, " \
    "F32BE, F32LE, S32BE, S32LE, U32BE, U32LE, " \
    "S24_32BE, S24_32LE, U24_32BE, U24_32LE, " \
    "S24BE, S24LE, U24BE, U24LE, " \
    "S20BE, S20LE, U20BE, U20LE, " \
    "S18BE, S18LE, U18BE, U18LE, " \
    "S16BE, S16LE, U16BE, U16LE, " \
    "S8, U8 }"
#elif G_BYTE_ORDER == G_LITTLE_ENDIAN
#define GST_AUDIO_FORMATS_ALL "{ F64LE, F64BE, " \
    "F32LE, F32BE, S32LE, S32BE, U32LE, U32BE, " \
    "S24_32LE, S24_32BE, U24_32LE, U24_32BE, " \
    "S24LE, S24BE, U24LE, U24BE, " \
    "S20LE, S20BE, U20LE, U20BE, " \
    "S18LE, S18BE, U18LE, U18BE, " \
    "S16LE, S16BE, U16LE, U16BE, " \
    "S8, U8 }"
#endif

GST_AUDIO_API
const GstAudioFormat * gst_audio_formats_raw (guint * len);

/**
 * GST_AUDIO_CAPS_MAKE:
 * @format: string format that describes the sample layout, as string
 *     (e.g. "S16LE", "S8", etc.)
 *
 * Generic caps string for audio, for use in pad templates.
 */
#define GST_AUDIO_CAPS_MAKE(format)                                    \
    "audio/x-raw, "                                                    \
    "format = (string) " format ", "                                   \
    "rate = " GST_AUDIO_RATE_RANGE ", "                                \
    "channels = " GST_AUDIO_CHANNELS_RANGE

/**
 * GST_AUDIO_DEF_RATE:
 *
 * Standard sampling rate used in consumer audio.
 */
#define GST_AUDIO_DEF_RATE 44100
/**
 * GST_AUDIO_DEF_CHANNELS:
 *
 * Standard number of channels used in consumer audio.
 */
#define GST_AUDIO_DEF_CHANNELS 2
/**
 * GST_AUDIO_DEF_FORMAT:
 *
 * Standard format used in consumer audio.
 */
#define GST_AUDIO_DEF_FORMAT "S16LE"

/**
 * GstAudioLayout:
 * @GST_AUDIO_LAYOUT_INTERLEAVED: interleaved audio
 * @GST_AUDIO_LAYOUT_NON_INTERLEAVED: non-interleaved audio
 *
 * Layout of the audio samples for the different channels.
 */
typedef enum {
  GST_AUDIO_LAYOUT_INTERLEAVED = 0,
  GST_AUDIO_LAYOUT_NON_INTERLEAVED
} GstAudioLayout;

GST_AUDIO_API
GstCaps * gst_audio_make_raw_caps (const GstAudioFormat formats[], guint len,
                                   GstAudioLayout layout);

G_END_DECLS

#endif /* __GST_AUDIO_FORMAT_H__ */
