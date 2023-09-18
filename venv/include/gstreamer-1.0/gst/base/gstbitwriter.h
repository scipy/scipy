/*
 *  gstbitwriter.h - bitstream writer
 *
 *  Copyright (C) 2013 Intel Corporation
 *  Copyright (C) 2018 Igalia, S. L.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301 USA
 */

#ifndef GST_BIT_WRITER_H
#define GST_BIT_WRITER_H

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

#include <string.h>

G_BEGIN_DECLS

#define GST_BIT_WRITER_DATA(writer)     ((writer)->data)
#define GST_BIT_WRITER_BIT_SIZE(writer) ((writer)->bit_size)
#define GST_BIT_WRITER(writer)          ((GstBitWriter *) (writer))

typedef struct _GstBitWriter GstBitWriter;

/**
 * GstBitWriter:
 * @data: Allocated @data for bit writer to write
 * @bit_size: Size of written @data in bits
 *
 * A bit writer instance.
 *
 * Since: 1.16
 */
struct _GstBitWriter
{
  guint8 *data;
  guint bit_size;

  /*< private >*/
  guint bit_capacity; /* Capacity of the allocated data */
  gboolean auto_grow; /* Whether space can auto grow */
  gboolean owned;
  gpointer _gst_reserved[GST_PADDING];
};

GST_BASE_API
GstBitWriter *  gst_bit_writer_new              (void) G_GNUC_MALLOC;

GST_BASE_API
GstBitWriter *  gst_bit_writer_new_with_size    (guint32 size, gboolean fixed) G_GNUC_MALLOC;

GST_BASE_API
GstBitWriter *  gst_bit_writer_new_with_data    (guint8 *data, guint size,
						 gboolean initialized) G_GNUC_MALLOC;

GST_BASE_API
void            gst_bit_writer_free             (GstBitWriter *bitwriter);

GST_BASE_API
guint8 *        gst_bit_writer_free_and_get_data (GstBitWriter *bitwriter);

GST_BASE_API
GstBuffer *     gst_bit_writer_free_and_get_buffer (GstBitWriter *bitwriter);

GST_BASE_API
void            gst_bit_writer_init             (GstBitWriter *bitwriter);

GST_BASE_API
void            gst_bit_writer_init_with_size   (GstBitWriter *bitwriter,
						 guint32 size, gboolean fixed);

GST_BASE_API
void            gst_bit_writer_init_with_data   (GstBitWriter *bitwriter,  guint8 *data,
						 guint size, gboolean initialized);

GST_BASE_API
void            gst_bit_writer_reset            (GstBitWriter *bitwriter);

GST_BASE_API
guint8 *        gst_bit_writer_reset_and_get_data (GstBitWriter *bitwriter);

GST_BASE_API
GstBuffer *     gst_bit_writer_reset_and_get_buffer (GstBitWriter *bitwriter);

GST_BASE_API
guint           gst_bit_writer_get_size         (const GstBitWriter *bitwriter);

GST_BASE_API
guint8 *        gst_bit_writer_get_data         (const GstBitWriter *bitwriter);

GST_BASE_API
gboolean        gst_bit_writer_set_pos          (GstBitWriter *bitwriter, guint pos);

GST_BASE_API
guint           gst_bit_writer_get_remaining    (const GstBitWriter *bitwriter);

GST_BASE_API
gboolean        gst_bit_writer_put_bits_uint8   (GstBitWriter *bitwriter, guint8 value,
						 guint nbits);

GST_BASE_API
gboolean        gst_bit_writer_put_bits_uint16  (GstBitWriter *bitwriter, guint16 value,
						 guint nbits);

GST_BASE_API
gboolean        gst_bit_writer_put_bits_uint32  (GstBitWriter *bitwriter, guint32 value,
						 guint nbits);

GST_BASE_API
gboolean        gst_bit_writer_put_bits_uint64  (GstBitWriter *bitwriter, guint64 value,
						 guint nbits);

GST_BASE_API
gboolean        gst_bit_writer_put_bytes        (GstBitWriter *bitwriter, const guint8 *data,
						 guint nbytes);

GST_BASE_API
gboolean        gst_bit_writer_align_bytes      (GstBitWriter *bitwriter, guint8 trailing_bit);

static const guint8 _gst_bit_writer_bit_filling_mask[9] = {
    0x00, 0x01, 0x03, 0x07,
    0x0F, 0x1F, 0x3F, 0x7F,
    0xFF
};

/* Aligned to 256 bytes */
#define __GST_BITS_WRITER_ALIGNMENT_MASK 2047
#define __GST_BITS_WRITER_ALIGNED(bitsize)                   \
    (((bitsize) + __GST_BITS_WRITER_ALIGNMENT_MASK)&(~__GST_BITS_WRITER_ALIGNMENT_MASK))

static inline gboolean
_gst_bit_writer_check_remaining (GstBitWriter * bitwriter, guint32 bits)
{
  guint32 new_bit_size = bits + bitwriter->bit_size;
  guint32 clear_pos;

  g_assert (bitwriter->bit_size <= bitwriter->bit_capacity);
  if (new_bit_size <= bitwriter->bit_capacity)
    return TRUE;

  if (!bitwriter->auto_grow)
    return FALSE;

  /* auto grow space */
  new_bit_size = __GST_BITS_WRITER_ALIGNED (new_bit_size);
  g_assert (new_bit_size
      && ((new_bit_size & __GST_BITS_WRITER_ALIGNMENT_MASK) == 0));
  clear_pos = ((bitwriter->bit_size + 7) >> 3);
  bitwriter->data = (guint8 *) g_realloc (bitwriter->data, (new_bit_size >> 3));
  memset (bitwriter->data + clear_pos, 0, (new_bit_size >> 3) - clear_pos);
  bitwriter->bit_capacity = new_bit_size;
  return TRUE;
}

#undef __GST_BITS_WRITER_ALIGNMENT_MASK
#undef __GST_BITS_WRITER_ALIGNED

#define __GST_BIT_WRITER_WRITE_BITS_UNCHECKED(bits) \
static inline void \
gst_bit_writer_put_bits_uint##bits##_unchecked( \
    GstBitWriter *bitwriter, \
    guint##bits value, \
    guint nbits \
) \
{ \
    guint byte_pos, bit_offset; \
    guint8  *cur_byte; \
    guint fill_bits; \
    \
    byte_pos = (bitwriter->bit_size >> 3); \
    bit_offset = (bitwriter->bit_size & 0x07); \
    cur_byte = bitwriter->data + byte_pos; \
    g_assert (nbits <= bits); \
    g_assert( bit_offset < 8 && \
            bitwriter->bit_size <= bitwriter->bit_capacity); \
    \
    while (nbits) { \
        fill_bits = ((8 - bit_offset) < nbits ? (8 - bit_offset) : nbits); \
        nbits -= fill_bits; \
        bitwriter->bit_size += fill_bits; \
        \
        *cur_byte |= (((value >> nbits) & _gst_bit_writer_bit_filling_mask[fill_bits]) \
                      << (8 - bit_offset - fill_bits)); \
        ++cur_byte; \
        bit_offset = 0; \
    } \
    g_assert(cur_byte <= \
           (bitwriter->data + (bitwriter->bit_capacity >> 3))); \
}

__GST_BIT_WRITER_WRITE_BITS_UNCHECKED (8)
__GST_BIT_WRITER_WRITE_BITS_UNCHECKED (16)
__GST_BIT_WRITER_WRITE_BITS_UNCHECKED (32)
__GST_BIT_WRITER_WRITE_BITS_UNCHECKED (64)
#undef __GST_BIT_WRITER_WRITE_BITS_UNCHECKED

static inline guint
gst_bit_writer_get_size_unchecked (const GstBitWriter * bitwriter)
{
  return GST_BIT_WRITER_BIT_SIZE (bitwriter);
}

static inline guint8 *
gst_bit_writer_get_data_unchecked (const GstBitWriter * bitwriter)
{
  return GST_BIT_WRITER_DATA (bitwriter);
}

static inline gboolean
gst_bit_writer_set_pos_unchecked (GstBitWriter * bitwriter, guint pos)
{
  GST_BIT_WRITER_BIT_SIZE (bitwriter) = pos;
  return TRUE;
}

static inline guint
gst_bit_writer_get_remaining_unchecked (const GstBitWriter * bitwriter)
{
  return bitwriter->bit_capacity - bitwriter->bit_size;
}

static inline void
gst_bit_writer_put_bytes_unchecked (GstBitWriter * bitwriter,
    const guint8 * data, guint nbytes)
{
  if ((bitwriter->bit_size & 0x07) == 0) {
    memcpy (&bitwriter->data[bitwriter->bit_size >> 3], data, nbytes);
    bitwriter->bit_size += (nbytes << 3);
  } else {
    g_assert (0);
    while (nbytes) {
      gst_bit_writer_put_bits_uint8_unchecked (bitwriter, *data, 8);
      --nbytes;
      ++data;
    }
  }
}

static inline void
gst_bit_writer_align_bytes_unchecked (GstBitWriter * bitwriter,
    guint8 trailing_bit)
{
  guint32 bit_offset, bit_left;
  guint8 value = 0;

  bit_offset = (bitwriter->bit_size & 0x07);
  if (!bit_offset)
    return;

  bit_left = 8 - bit_offset;
  if (trailing_bit)
    value = _gst_bit_writer_bit_filling_mask[bit_left];
  gst_bit_writer_put_bits_uint8_unchecked (bitwriter, value, bit_left);
}

#define __GST_BIT_WRITER_WRITE_BITS_INLINE(bits) \
static inline gboolean \
_gst_bit_writer_put_bits_uint##bits##_inline( \
    GstBitWriter *bitwriter, \
    guint##bits value, \
    guint nbits \
) \
{ \
    g_return_val_if_fail(bitwriter != NULL, FALSE); \
    g_return_val_if_fail(nbits != 0, FALSE); \
    g_return_val_if_fail(nbits <= bits, FALSE); \
    \
    if (!_gst_bit_writer_check_remaining(bitwriter, nbits)) \
        return FALSE; \
    gst_bit_writer_put_bits_uint##bits##_unchecked(bitwriter, value, nbits); \
    return TRUE; \
}

__GST_BIT_WRITER_WRITE_BITS_INLINE (8)
__GST_BIT_WRITER_WRITE_BITS_INLINE (16)
__GST_BIT_WRITER_WRITE_BITS_INLINE (32)
__GST_BIT_WRITER_WRITE_BITS_INLINE (64)
#undef __GST_BIT_WRITER_WRITE_BITS_INLINE

static inline guint
_gst_bit_writer_get_size_inline (const GstBitWriter * bitwriter)
{
  g_return_val_if_fail (bitwriter != NULL, 0);

  return gst_bit_writer_get_size_unchecked (bitwriter);
}

static inline guint8 *
_gst_bit_writer_get_data_inline (const GstBitWriter * bitwriter)
{
  g_return_val_if_fail (bitwriter != NULL, NULL);

  return gst_bit_writer_get_data_unchecked (bitwriter);
}

static inline gboolean
_gst_bit_writer_set_pos_inline (GstBitWriter * bitwriter, guint pos)
{
  g_return_val_if_fail (bitwriter != NULL, FALSE);
  g_return_val_if_fail (pos <= bitwriter->bit_capacity, FALSE);

  return gst_bit_writer_set_pos_unchecked (bitwriter, pos);
}

static inline guint
_gst_bit_writer_get_remaining_inline (const GstBitWriter * bitwriter)
{
  g_return_val_if_fail (bitwriter != NULL, 0);
  g_return_val_if_fail (bitwriter->bit_size < bitwriter->bit_capacity, 0);

  return gst_bit_writer_get_remaining_unchecked (bitwriter);
}

static inline gboolean
_gst_bit_writer_put_bytes_inline (GstBitWriter * bitwriter,
    const guint8 * data, guint nbytes)
{
  g_return_val_if_fail (bitwriter != NULL, FALSE);
  g_return_val_if_fail (data != NULL, FALSE);
  g_return_val_if_fail (nbytes, FALSE);

  if (!_gst_bit_writer_check_remaining (bitwriter, nbytes * 8))
    return FALSE;

  gst_bit_writer_put_bytes_unchecked (bitwriter, data, nbytes);
  return TRUE;
}

static inline gboolean
_gst_bit_writer_align_bytes_inline (GstBitWriter * bitwriter,
    guint8 trailing_bit)
{
  g_return_val_if_fail (bitwriter != NULL, FALSE);
  g_return_val_if_fail ((trailing_bit == 0 || trailing_bit == 1), FALSE);
  g_return_val_if_fail (((bitwriter->bit_size + 7) & (~7)) <=
      bitwriter->bit_capacity, FALSE);

  gst_bit_writer_align_bytes_unchecked (bitwriter, trailing_bit);
  return TRUE;
}

#ifndef GST_BIT_WRITER_DISABLE_INLINES
#define gst_bit_writer_get_size(bitwriter) \
    _gst_bit_writer_get_size_inline(bitwriter)
#define gst_bit_writer_get_data(bitwriter) \
    _gst_bit_writer_get_data_inline(bitwriter)
#define gst_bit_writer_set_pos(bitwriter, pos) \
    G_LIKELY (_gst_bit_writer_set_pos_inline (bitwriter, pos))
#define gst_bit_writer_get_remaining(bitwriter) \
    _gst_bit_writer_get_remaining_inline(bitwriter)

#define gst_bit_writer_put_bits_uint8(bitwriter, value, nbits) \
    G_LIKELY (_gst_bit_writer_put_bits_uint8_inline (bitwriter, value, nbits))
#define gst_bit_writer_put_bits_uint16(bitwriter, value, nbits) \
    G_LIKELY (_gst_bit_writer_put_bits_uint16_inline (bitwriter, value, nbits))
#define gst_bit_writer_put_bits_uint32(bitwriter, value, nbits) \
    G_LIKELY (_gst_bit_writer_put_bits_uint32_inline (bitwriter, value, nbits))
#define gst_bit_writer_put_bits_uint64(bitwriter, value, nbits) \
    G_LIKELY (_gst_bit_writer_put_bits_uint64_inline (bitwriter, value, nbits))

#define gst_bit_writer_put_bytes(bitwriter, data, nbytes) \
    G_LIKELY (_gst_bit_writer_put_bytes_inline (bitwriter, data, nbytes))

#define gst_bit_writer_align_bytes(bitwriter, trailing_bit) \
    G_LIKELY (_gst_bit_writer_align_bytes_inline(bitwriter, trailing_bit))
#endif

G_END_DECLS

#endif /* GST_BIT_WRITER_H */
