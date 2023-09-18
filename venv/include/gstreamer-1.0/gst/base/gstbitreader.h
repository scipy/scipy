/* GStreamer
 *
 * Copyright (C) 2008 Sebastian Dröge <sebastian.droege@collabora.co.uk>.
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

#ifndef __GST_BIT_READER_H__
#define __GST_BIT_READER_H__

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

/* FIXME: inline functions */

G_BEGIN_DECLS

#define GST_BIT_READER(reader) ((GstBitReader *) (reader))

/**
 * GstBitReader:
 * @data: (array length=size): Data from which the bit reader will
 *   read
 * @size: Size of @data in bytes
 * @byte: Current byte position
 * @bit: Bit position in the current byte
 *
 * A bit reader instance.
 */
typedef struct {
  const guint8 *data;
  guint size;

  guint byte;  /* Byte position */
  guint bit;   /* Bit position in the current byte */

  /* < private > */
  gpointer _gst_reserved[GST_PADDING];
} GstBitReader;

GST_BASE_API
GstBitReader *  gst_bit_reader_new              (const guint8 *data, guint size) G_GNUC_MALLOC;

GST_BASE_API
void            gst_bit_reader_free             (GstBitReader *reader);

GST_BASE_API
void            gst_bit_reader_init             (GstBitReader *reader, const guint8 *data, guint size);

GST_BASE_API
gboolean        gst_bit_reader_set_pos          (GstBitReader *reader, guint pos);

GST_BASE_API
guint           gst_bit_reader_get_pos          (const GstBitReader *reader);

GST_BASE_API
guint           gst_bit_reader_get_remaining    (const GstBitReader *reader);

GST_BASE_API
guint           gst_bit_reader_get_size         (const GstBitReader *reader);

GST_BASE_API
gboolean        gst_bit_reader_skip             (GstBitReader *reader, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_skip_to_byte     (GstBitReader *reader);

GST_BASE_API
gboolean        gst_bit_reader_get_bits_uint8   (GstBitReader *reader, guint8 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_get_bits_uint16  (GstBitReader *reader, guint16 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_get_bits_uint32  (GstBitReader *reader, guint32 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_get_bits_uint64  (GstBitReader *reader, guint64 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_peek_bits_uint8  (const GstBitReader *reader, guint8 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_peek_bits_uint16 (const GstBitReader *reader, guint16 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_peek_bits_uint32 (const GstBitReader *reader, guint32 *val, guint nbits);

GST_BASE_API
gboolean        gst_bit_reader_peek_bits_uint64 (const GstBitReader *reader, guint64 *val, guint nbits);

/**
 * GST_BIT_READER_INIT:
 * @data: Data from which the #GstBitReader should read
 * @size: Size of @data in bytes
 *
 * A #GstBitReader must be initialized with this macro, before it can be
 * used. This macro can used be to initialize a variable, but it cannot
 * be assigned to a variable. In that case you have to use
 * gst_bit_reader_init().
 */
#define GST_BIT_READER_INIT(data, size) {data, size, 0, 0}

/* Unchecked variants */

static inline void
gst_bit_reader_skip_unchecked (GstBitReader * reader, guint nbits)
{
  reader->bit += nbits;
  reader->byte += reader->bit / 8;
  reader->bit = reader->bit % 8;
}

static inline void
gst_bit_reader_skip_to_byte_unchecked (GstBitReader * reader)
{
  if (reader->bit) {
    reader->bit = 0;
    reader->byte++;
  }
}

#define __GST_BIT_READER_READ_BITS_UNCHECKED(bits) \
static inline guint##bits \
gst_bit_reader_peek_bits_uint##bits##_unchecked (const GstBitReader *reader, guint nbits) \
{ \
  guint##bits ret = 0; \
  const guint8 *data; \
  guint byte, bit; \
  \
  data = reader->data; \
  byte = reader->byte; \
  bit = reader->bit; \
  \
  while (nbits > 0) { \
    guint toread = MIN (nbits, 8 - bit); \
    \
    ret <<= toread; \
    ret |= (data[byte] & (0xff >> bit)) >> (8 - toread - bit); \
    \
    bit += toread; \
    if (bit >= 8) { \
      byte++; \
      bit = 0; \
    } \
    nbits -= toread; \
  } \
  \
  return ret; \
} \
\
static inline guint##bits \
gst_bit_reader_get_bits_uint##bits##_unchecked (GstBitReader *reader, guint nbits) \
{ \
  guint##bits ret; \
  \
  ret = gst_bit_reader_peek_bits_uint##bits##_unchecked (reader, nbits); \
  \
  gst_bit_reader_skip_unchecked (reader, nbits); \
  \
  return ret; \
}

__GST_BIT_READER_READ_BITS_UNCHECKED (8)
__GST_BIT_READER_READ_BITS_UNCHECKED (16)
__GST_BIT_READER_READ_BITS_UNCHECKED (32)
__GST_BIT_READER_READ_BITS_UNCHECKED (64)

#undef __GST_BIT_READER_READ_BITS_UNCHECKED

/* unchecked variants -- do not use */

static inline guint
_gst_bit_reader_get_size_unchecked (const GstBitReader * reader)
{
  return reader->size * 8;
}

static inline guint
_gst_bit_reader_get_pos_unchecked (const GstBitReader * reader)
{
  return reader->byte * 8 + reader->bit;
}

static inline guint
_gst_bit_reader_get_remaining_unchecked (const GstBitReader * reader)
{
  return reader->size * 8 - (reader->byte * 8 + reader->bit);
}

/* inlined variants -- do not use directly */
static inline guint
_gst_bit_reader_get_size_inline (const GstBitReader * reader)
{
  g_return_val_if_fail (reader != NULL, 0);

  return _gst_bit_reader_get_size_unchecked (reader);
}

static inline guint
_gst_bit_reader_get_pos_inline (const GstBitReader * reader)
{
  g_return_val_if_fail (reader != NULL, 0);

  return _gst_bit_reader_get_pos_unchecked (reader);
}

static inline guint
_gst_bit_reader_get_remaining_inline (const GstBitReader * reader)
{
  g_return_val_if_fail (reader != NULL, 0);

  return _gst_bit_reader_get_remaining_unchecked (reader);
}

static inline gboolean
_gst_bit_reader_skip_inline (GstBitReader * reader, guint nbits)
{
  g_return_val_if_fail (reader != NULL, FALSE);

  if (_gst_bit_reader_get_remaining_unchecked (reader) < nbits)
    return FALSE;

  gst_bit_reader_skip_unchecked (reader, nbits);

  return TRUE;
}

static inline gboolean
_gst_bit_reader_skip_to_byte_inline (GstBitReader * reader)
{
  g_return_val_if_fail (reader != NULL, FALSE);

  if (reader->byte > reader->size)
    return FALSE;

  gst_bit_reader_skip_to_byte_unchecked (reader);

  return TRUE;
}

#define __GST_BIT_READER_READ_BITS_INLINE(bits) \
static inline gboolean \
_gst_bit_reader_get_bits_uint##bits##_inline (GstBitReader *reader, guint##bits *val, guint nbits) \
{ \
  g_return_val_if_fail (reader != NULL, FALSE); \
  g_return_val_if_fail (val != NULL, FALSE); \
  g_return_val_if_fail (nbits <= bits, FALSE); \
  \
  if (_gst_bit_reader_get_remaining_unchecked (reader) < nbits) \
    return FALSE; \
\
  *val = gst_bit_reader_get_bits_uint##bits##_unchecked (reader, nbits); \
  return TRUE; \
} \
\
static inline gboolean \
_gst_bit_reader_peek_bits_uint##bits##_inline (const GstBitReader *reader, guint##bits *val, guint nbits) \
{ \
  g_return_val_if_fail (reader != NULL, FALSE); \
  g_return_val_if_fail (val != NULL, FALSE); \
  g_return_val_if_fail (nbits <= bits, FALSE); \
  \
  if (_gst_bit_reader_get_remaining_unchecked (reader) < nbits) \
    return FALSE; \
\
  *val = gst_bit_reader_peek_bits_uint##bits##_unchecked (reader, nbits); \
  return TRUE; \
}

__GST_BIT_READER_READ_BITS_INLINE (8)
__GST_BIT_READER_READ_BITS_INLINE (16)
__GST_BIT_READER_READ_BITS_INLINE (32)
__GST_BIT_READER_READ_BITS_INLINE (64)

#undef __GST_BIT_READER_READ_BITS_INLINE

#ifndef GST_BIT_READER_DISABLE_INLINES

#define gst_bit_reader_get_size(reader) \
    _gst_bit_reader_get_size_inline (reader)
#define gst_bit_reader_get_pos(reader) \
    _gst_bit_reader_get_pos_inline (reader)
#define gst_bit_reader_get_remaining(reader) \
    _gst_bit_reader_get_remaining_inline (reader)

/* we use defines here so we can add the G_LIKELY() */

#define gst_bit_reader_skip(reader, nbits)\
    G_LIKELY (_gst_bit_reader_skip_inline(reader, nbits))
#define gst_bit_reader_skip_to_byte(reader)\
    G_LIKELY (_gst_bit_reader_skip_to_byte_inline(reader))

#define gst_bit_reader_get_bits_uint8(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_get_bits_uint8_inline (reader, val, nbits))
#define gst_bit_reader_get_bits_uint16(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_get_bits_uint16_inline (reader, val, nbits))
#define gst_bit_reader_get_bits_uint32(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_get_bits_uint32_inline (reader, val, nbits))
#define gst_bit_reader_get_bits_uint64(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_get_bits_uint64_inline (reader, val, nbits))

#define gst_bit_reader_peek_bits_uint8(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_peek_bits_uint8_inline (reader, val, nbits))
#define gst_bit_reader_peek_bits_uint16(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_peek_bits_uint16_inline (reader, val, nbits))
#define gst_bit_reader_peek_bits_uint32(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_peek_bits_uint32_inline (reader, val, nbits))
#define gst_bit_reader_peek_bits_uint64(reader, val, nbits) \
    G_LIKELY (_gst_bit_reader_peek_bits_uint64_inline (reader, val, nbits))
#endif

G_END_DECLS

#endif /* __GST_BIT_READER_H__ */
