/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *                    2002 Thomas Vander Stichele <thomas@apestaart.org>
 *
 * gstutils.h: Header for various utility functions
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


#ifndef __GST_UTILS_H__
#define __GST_UTILS_H__

#include <glib.h>
#include <gst/gstconfig.h>
#include <gst/gstbin.h>
#include <gst/gstparse.h>

G_BEGIN_DECLS

GST_API
void            gst_util_set_value_from_string  (GValue *value, const gchar *value_str);

GST_API
void            gst_util_set_object_arg         (GObject *object, const gchar *name, const gchar *value);

GST_API
gboolean        gst_util_set_object_array       (GObject * object, const gchar * name,
                                                 const GValueArray * array);
GST_API
gboolean        gst_util_get_object_array       (GObject * object, const gchar * name,
                                                 GValueArray ** array);
GST_API
void            gst_util_dump_mem               (const guchar *mem, guint size);

GST_API
void            gst_util_dump_buffer            (GstBuffer * buf);

GST_API
guint64         gst_util_gdouble_to_guint64     (gdouble value)  G_GNUC_CONST;

GST_API
gdouble         gst_util_guint64_to_gdouble     (guint64 value)  G_GNUC_CONST;

/**
 * gst_guint64_to_gdouble:
 * @value: the #guint64 value to convert
 *
 * Convert @value to a gdouble.
 *
 * Returns: @value converted to a #gdouble.
 */

/**
 * gst_gdouble_to_guint64:
 * @value: the #gdouble value to convert
 *
 * Convert @value to a guint64.
 *
 * Returns: @value converted to a #guint64.
 */
#ifdef WIN32
#define         gst_gdouble_to_guint64(value)   gst_util_gdouble_to_guint64(value)
#define         gst_guint64_to_gdouble(value)   gst_util_guint64_to_gdouble(value)
#else
#define         gst_gdouble_to_guint64(value)   ((guint64) (value))
#define         gst_guint64_to_gdouble(value)   ((gdouble) (value))
#endif

GST_API
guint64         gst_util_uint64_scale           (guint64 val, guint64 num, guint64 denom);

GST_API
guint64         gst_util_uint64_scale_round     (guint64 val, guint64 num, guint64 denom);

GST_API
guint64         gst_util_uint64_scale_ceil      (guint64 val, guint64 num, guint64 denom);

GST_API
guint64         gst_util_uint64_scale_int       (guint64 val, gint num, gint denom);

GST_API
guint64         gst_util_uint64_scale_int_round (guint64 val, gint num, gint denom);

GST_API
guint64         gst_util_uint64_scale_int_ceil  (guint64 val, gint num, gint denom);

/**
 * GST_SEQNUM_INVALID:
 *
 * A value which is guaranteed to never be returned by
 * gst_util_seqnum_next().
 *
 * Can be used as a default value in variables used to store seqnum.
 *
 * Since: 1.14
 */
#define GST_SEQNUM_INVALID (0)

GST_API
guint32         gst_util_seqnum_next            (void);

GST_API
gint32          gst_util_seqnum_compare         (guint32 s1, guint32 s2);

/**
 * GST_GROUP_ID_INVALID:
 *
 * A value which is guaranteed to never be returned by
 * gst_util_group_id_next().
 *
 * Can be used as a default value in variables used to store group_id.
 *
 * Since: 1.14
 */
#define GST_GROUP_ID_INVALID (0)

GST_API
guint           gst_util_group_id_next          (void);

/**
 * GST_CALL_PARENT:
 * @parent_class_cast: the name of the class cast macro for the parent type
 * @name: name of the function to call
 * @args: arguments enclosed in '( )'
 *
 * Just call the parent handler.  This assumes that there is a variable
 * named parent_class that points to the (duh!) parent class.  Note that
 * this macro is not to be used with things that return something, use
 * the _WITH_DEFAULT version for that
 */
#define GST_CALL_PARENT(parent_class_cast, name, args)                  \
        ((parent_class_cast(parent_class)->name != NULL) ?              \
         parent_class_cast(parent_class)->name args : (void) 0)

/**
 * GST_CALL_PARENT_WITH_DEFAULT:
 * @parent_class_cast: the name of the class cast macro for the parent type
 * @name: name of the function to call
 * @args: arguments enclosed in '( )'
 * @def_return: default result
 *
 * Same as GST_CALL_PARENT(), but in case there is no implementation, it
 * evaluates to @def_return.
 */
#define GST_CALL_PARENT_WITH_DEFAULT(parent_class_cast, name, args, def_return)\
        ((parent_class_cast(parent_class)->name != NULL) ?              \
         parent_class_cast(parent_class)->name args : def_return)

/* Define PUT and GET functions for unaligned memory */
#define _GST_GET(__data, __idx, __size, __shift) \
    (((guint##__size) (((const guint8 *) (__data))[__idx])) << (__shift))

#define _GST_PUT(__data, __idx, __size, __shift, __num) \
    (((guint8 *) (__data))[__idx] = (((guint##__size) (__num)) >> (__shift)) & 0xff)

#ifndef __GTK_DOC_IGNORE__
#if GST_HAVE_UNALIGNED_ACCESS
static inline guint16 __gst_fast_read16(const guint8 *v) {
  return *(const guint16*)(const void*)(v);
}
static inline guint32 __gst_fast_read32(const guint8 *v) {
  return *(const guint32*)(const void*)(v);
}
static inline guint64 __gst_fast_read64(const guint8 *v) {
  return *(const guint64*)(const void*)(v);
}
static inline guint16 __gst_fast_read_swap16(const guint8 *v) {
  return GUINT16_SWAP_LE_BE(*(const guint16*)(const void*)(v));
}
static inline guint32 __gst_fast_read_swap32(const guint8 *v) {
  return GUINT32_SWAP_LE_BE(*(const guint32*)(const void*)(v));
}
static inline guint64 __gst_fast_read_swap64(const guint8 *v) {
  return GUINT64_SWAP_LE_BE(*(const guint64*)(const void*)(v));
}
# define _GST_FAST_READ(s, d) __gst_fast_read##s((const guint8 *)(d))
# define _GST_FAST_READ_SWAP(s, d) __gst_fast_read_swap##s((const guint8 *)(d))

static inline void __gst_fast_write16 (guint8 *p, guint16 v) {
  *(guint16*)(void*)(p) = v;
}
static inline void __gst_fast_write32 (guint8 *p, guint32 v) {
  *(guint32*)(void*)(p) = v;
}
static inline void __gst_fast_write64 (guint8 *p, guint64 v) {
  *(guint64*)(void*)(p) = v;
}
static inline void __gst_fast_write_swap16 (guint8 *p, guint16 v) {
  *(guint16*)(void*)(p) = GUINT16_SWAP_LE_BE (v);
}
static inline void __gst_fast_write_swap32 (guint8 *p, guint32 v) {
  *(guint32*)(void*)(p) = GUINT32_SWAP_LE_BE (v);
}
static inline void __gst_fast_write_swap64 (guint8 *p, guint64 v) {
  *(guint64*)(void*)(p) = GUINT64_SWAP_LE_BE (v);
}
# define _GST_FAST_WRITE(s, d, v) __gst_fast_write##s((guint8 *)(d), (v))
# define _GST_FAST_WRITE_SWAP(s, d, v) __gst_fast_write_swap##s((guint8 *)(d), (v))
#endif
#endif


/**
 * GST_READ_UINT64_BE:
 * @data: memory location
 *
 * Read a 64 bit unsigned integer value in big endian format from the memory buffer.
 */

/**
 * GST_READ_UINT64_LE:
 * @data: memory location
 *
 * Read a 64 bit unsigned integer value in little endian format from the memory buffer.
 */
#if GST_HAVE_UNALIGNED_ACCESS
# if (G_BYTE_ORDER == G_BIG_ENDIAN)
#  define GST_READ_UINT64_BE(data)      _GST_FAST_READ (64, data)
#  define GST_READ_UINT64_LE(data)      _GST_FAST_READ_SWAP (64, data)
# else
#  define GST_READ_UINT64_BE(data)      _GST_FAST_READ_SWAP (64, data)
#  define GST_READ_UINT64_LE(data)      _GST_FAST_READ (64, data)
# endif
#else
#define _GST_READ_UINT64_BE(data)	(_GST_GET (data, 0, 64, 56) | \
					 _GST_GET (data, 1, 64, 48) | \
					 _GST_GET (data, 2, 64, 40) | \
					 _GST_GET (data, 3, 64, 32) | \
					 _GST_GET (data, 4, 64, 24) | \
					 _GST_GET (data, 5, 64, 16) | \
					 _GST_GET (data, 6, 64,  8) | \
					 _GST_GET (data, 7, 64,  0))

#define _GST_READ_UINT64_LE(data)	(_GST_GET (data, 7, 64, 56) | \
					 _GST_GET (data, 6, 64, 48) | \
					 _GST_GET (data, 5, 64, 40) | \
					 _GST_GET (data, 4, 64, 32) | \
					 _GST_GET (data, 3, 64, 24) | \
					 _GST_GET (data, 2, 64, 16) | \
					 _GST_GET (data, 1, 64,  8) | \
					 _GST_GET (data, 0, 64,  0))

#define GST_READ_UINT64_BE(data) __gst_slow_read64_be((const guint8 *)(data))
static inline guint64 __gst_slow_read64_be (const guint8 * data) {
  return _GST_READ_UINT64_BE (data);
}
#define GST_READ_UINT64_LE(data) __gst_slow_read64_le((const guint8 *)(data))
static inline guint64 __gst_slow_read64_le (const guint8 * data) {
  return _GST_READ_UINT64_LE (data);
}
#endif

/**
 * GST_READ_UINT32_BE:
 * @data: memory location
 *
 * Read a 32 bit unsigned integer value in big endian format from the memory buffer.
 */

/**
 * GST_READ_UINT32_LE:
 * @data: memory location
 *
 * Read a 32 bit unsigned integer value in little endian format from the memory buffer.
 */
#if GST_HAVE_UNALIGNED_ACCESS
# if (G_BYTE_ORDER == G_BIG_ENDIAN)
#  define GST_READ_UINT32_BE(data)      _GST_FAST_READ (32, data)
#  define GST_READ_UINT32_LE(data)      _GST_FAST_READ_SWAP (32, data)
# else
#  define GST_READ_UINT32_BE(data)      _GST_FAST_READ_SWAP (32, data)
#  define GST_READ_UINT32_LE(data)      _GST_FAST_READ (32, data)
# endif
#else
#define _GST_READ_UINT32_BE(data)	(_GST_GET (data, 0, 32, 24) | \
					 _GST_GET (data, 1, 32, 16) | \
					 _GST_GET (data, 2, 32,  8) | \
					 _GST_GET (data, 3, 32,  0))

#define _GST_READ_UINT32_LE(data)	(_GST_GET (data, 3, 32, 24) | \
					 _GST_GET (data, 2, 32, 16) | \
					 _GST_GET (data, 1, 32,  8) | \
					 _GST_GET (data, 0, 32,  0))

#define GST_READ_UINT32_BE(data) __gst_slow_read32_be((const guint8 *)(data))
static inline guint32 __gst_slow_read32_be (const guint8 * data) {
  return _GST_READ_UINT32_BE (data);
}
#define GST_READ_UINT32_LE(data) __gst_slow_read32_le((const guint8 *)(data))
static inline guint32 __gst_slow_read32_le (const guint8 * data) {
  return _GST_READ_UINT32_LE (data);
}
#endif

/**
 * GST_READ_UINT24_BE:
 * @data: memory location
 *
 * Read a 24 bit unsigned integer value in big endian format from the memory buffer.
 */
#define _GST_READ_UINT24_BE(data)       (_GST_GET (data, 0, 32, 16) | \
                                         _GST_GET (data, 1, 32,  8) | \
                                         _GST_GET (data, 2, 32,  0))

#define GST_READ_UINT24_BE(data) __gst_slow_read24_be((const guint8 *)(data))
static inline guint32 __gst_slow_read24_be (const guint8 * data) {
  return _GST_READ_UINT24_BE (data);
}

/**
 * GST_READ_UINT24_LE:
 * @data: memory location
 *
 * Read a 24 bit unsigned integer value in little endian format from the memory buffer.
 */
#define _GST_READ_UINT24_LE(data)       (_GST_GET (data, 2, 32, 16) | \
                                         _GST_GET (data, 1, 32,  8) | \
                                         _GST_GET (data, 0, 32,  0))

#define GST_READ_UINT24_LE(data) __gst_slow_read24_le((const guint8 *)(data))
static inline guint32 __gst_slow_read24_le (const guint8 * data) {
  return _GST_READ_UINT24_LE (data);
}

/**
 * GST_READ_UINT16_BE:
 * @data: memory location
 *
 * Read a 16 bit unsigned integer value in big endian format from the memory buffer.
 */
/**
 * GST_READ_UINT16_LE:
 * @data: memory location
 *
 * Read a 16 bit unsigned integer value in little endian format from the memory buffer.
 */
#if GST_HAVE_UNALIGNED_ACCESS
# if (G_BYTE_ORDER == G_BIG_ENDIAN)
#  define GST_READ_UINT16_BE(data)      _GST_FAST_READ (16, data)
#  define GST_READ_UINT16_LE(data)      _GST_FAST_READ_SWAP (16, data)
# else
#  define GST_READ_UINT16_BE(data)      _GST_FAST_READ_SWAP (16, data)
#  define GST_READ_UINT16_LE(data)      _GST_FAST_READ (16, data)
# endif
#else
#define _GST_READ_UINT16_BE(data)	(_GST_GET (data, 0, 16,  8) | \
					 _GST_GET (data, 1, 16,  0))

#define _GST_READ_UINT16_LE(data)	(_GST_GET (data, 1, 16,  8) | \
					 _GST_GET (data, 0, 16,  0))

#define GST_READ_UINT16_BE(data) __gst_slow_read16_be((const guint8 *)(data))
static inline guint16 __gst_slow_read16_be (const guint8 * data) {
  return _GST_READ_UINT16_BE (data);
}
#define GST_READ_UINT16_LE(data) __gst_slow_read16_le((const guint8 *)(data))
static inline guint16 __gst_slow_read16_le (const guint8 * data) {
  return _GST_READ_UINT16_LE (data);
}
#endif

/**
 * GST_READ_UINT8:
 * @data: memory location
 *
 * Read an 8 bit unsigned integer value from the memory buffer.
 */
#define GST_READ_UINT8(data)            (_GST_GET (data, 0,  8,  0))

/**
 * GST_WRITE_UINT64_BE:
 * @data: memory location
 * @val: value to store
 *
 * Store a 64 bit unsigned integer value in big endian format into the memory buffer.
 */
/**
 * GST_WRITE_UINT64_LE:
 * @data: memory location
 * @val: value to store
 *
 * Store a 64 bit unsigned integer value in little endian format into the memory buffer.
 */
#if GST_HAVE_UNALIGNED_ACCESS
# if (G_BYTE_ORDER == G_BIG_ENDIAN)
#  define GST_WRITE_UINT64_BE(data,val) _GST_FAST_WRITE(64,data,val)
#  define GST_WRITE_UINT64_LE(data,val) _GST_FAST_WRITE_SWAP(64,data,val)
# else
#  define GST_WRITE_UINT64_BE(data,val) _GST_FAST_WRITE_SWAP(64,data,val)
#  define GST_WRITE_UINT64_LE(data,val) _GST_FAST_WRITE(64,data,val)
# endif
#else
#define GST_WRITE_UINT64_BE(data,val)   do { \
                                          gpointer __put_data = data; \
                                          guint64 __put_val = val; \
                                          _GST_PUT (__put_data, 0, 64, 56, __put_val); \
                                          _GST_PUT (__put_data, 1, 64, 48, __put_val); \
                                          _GST_PUT (__put_data, 2, 64, 40, __put_val); \
                                          _GST_PUT (__put_data, 3, 64, 32, __put_val); \
                                          _GST_PUT (__put_data, 4, 64, 24, __put_val); \
                                          _GST_PUT (__put_data, 5, 64, 16, __put_val); \
                                          _GST_PUT (__put_data, 6, 64,  8, __put_val); \
                                          _GST_PUT (__put_data, 7, 64,  0, __put_val); \
                                        } while (0)

#define GST_WRITE_UINT64_LE(data,val)   do { \
                                          gpointer __put_data = data; \
                                          guint64 __put_val = val; \
                                          _GST_PUT (__put_data, 0, 64,  0, __put_val); \
                                          _GST_PUT (__put_data, 1, 64,  8, __put_val); \
                                          _GST_PUT (__put_data, 2, 64, 16, __put_val); \
                                          _GST_PUT (__put_data, 3, 64, 24, __put_val); \
                                          _GST_PUT (__put_data, 4, 64, 32, __put_val); \
                                          _GST_PUT (__put_data, 5, 64, 40, __put_val); \
                                          _GST_PUT (__put_data, 6, 64, 48, __put_val); \
                                          _GST_PUT (__put_data, 7, 64, 56, __put_val); \
                                        } while (0)
#endif /* !GST_HAVE_UNALIGNED_ACCESS */

/**
 * GST_WRITE_UINT32_BE:
 * @data: memory location
 * @val: value to store
 *
 * Store a 32 bit unsigned integer value in big endian format into the memory buffer.
 */
/**
 * GST_WRITE_UINT32_LE:
 * @data: memory location
 * @val: value to store
 *
 * Store a 32 bit unsigned integer value in little endian format into the memory buffer.
 */
#if GST_HAVE_UNALIGNED_ACCESS
# if (G_BYTE_ORDER == G_BIG_ENDIAN)
#  define GST_WRITE_UINT32_BE(data,val) _GST_FAST_WRITE(32,data,val)
#  define GST_WRITE_UINT32_LE(data,val) _GST_FAST_WRITE_SWAP(32,data,val)
# else
#  define GST_WRITE_UINT32_BE(data,val) _GST_FAST_WRITE_SWAP(32,data,val)
#  define GST_WRITE_UINT32_LE(data,val) _GST_FAST_WRITE(32,data,val)
# endif
#else
#define GST_WRITE_UINT32_BE(data,val)   do { \
                                          gpointer __put_data = data; \
                                          guint32 __put_val = val; \
                                          _GST_PUT (__put_data, 0, 32, 24, __put_val); \
                                          _GST_PUT (__put_data, 1, 32, 16, __put_val); \
                                          _GST_PUT (__put_data, 2, 32,  8, __put_val); \
                                          _GST_PUT (__put_data, 3, 32,  0, __put_val); \
                                        } while (0)

#define GST_WRITE_UINT32_LE(data,val)   do { \
                                          gpointer __put_data = data; \
                                          guint32 __put_val = val; \
                                          _GST_PUT (__put_data, 0, 32,  0, __put_val); \
                                          _GST_PUT (__put_data, 1, 32,  8, __put_val); \
                                          _GST_PUT (__put_data, 2, 32, 16, __put_val); \
                                          _GST_PUT (__put_data, 3, 32, 24, __put_val); \
                                        } while (0)
#endif /* !GST_HAVE_UNALIGNED_ACCESS */

/**
 * GST_WRITE_UINT24_BE:
 * @data: memory location
 * @num: value to store
 *
 * Store a 24 bit unsigned integer value in big endian format into the memory buffer.
 */
#define GST_WRITE_UINT24_BE(data, num)  do { \
                                          gpointer __put_data = data; \
                                          guint32 __put_val = num; \
                                          _GST_PUT (__put_data, 0, 32,  16, __put_val); \
                                          _GST_PUT (__put_data, 1, 32,  8, __put_val); \
                                          _GST_PUT (__put_data, 2, 32,  0, __put_val); \
                                        } while (0)

/**
 * GST_WRITE_UINT24_LE:
 * @data: memory location
 * @num: value to store
 *
 * Store a 24 bit unsigned integer value in little endian format into the memory buffer.
 */
#define GST_WRITE_UINT24_LE(data, num)  do { \
                                          gpointer __put_data = data; \
                                          guint32 __put_val = num; \
                                          _GST_PUT (__put_data, 0, 32,  0, __put_val); \
                                          _GST_PUT (__put_data, 1, 32,  8, __put_val); \
                                          _GST_PUT (__put_data, 2, 32,  16, __put_val); \
                                        } while (0)

/**
 * GST_WRITE_UINT16_BE:
 * @data: memory location
 * @val: value to store
 *
 * Store a 16 bit unsigned integer value in big endian format into the memory buffer.
 */
/**
 * GST_WRITE_UINT16_LE:
 * @data: memory location
 * @val: value to store
 *
 * Store a 16 bit unsigned integer value in little endian format into the memory buffer.
 */
#if GST_HAVE_UNALIGNED_ACCESS
# if (G_BYTE_ORDER == G_BIG_ENDIAN)
#  define GST_WRITE_UINT16_BE(data,val) _GST_FAST_WRITE(16,data,val)
#  define GST_WRITE_UINT16_LE(data,val) _GST_FAST_WRITE_SWAP(16,data,val)
# else
#  define GST_WRITE_UINT16_BE(data,val) _GST_FAST_WRITE_SWAP(16,data,val)
#  define GST_WRITE_UINT16_LE(data,val) _GST_FAST_WRITE(16,data,val)
# endif
#else
#define GST_WRITE_UINT16_BE(data,val)   do { \
                                          gpointer __put_data = data; \
                                          guint16 __put_val = val; \
                                          _GST_PUT (__put_data, 0, 16,  8, __put_val); \
                                          _GST_PUT (__put_data, 1, 16,  0, __put_val); \
                                        } while (0)

#define GST_WRITE_UINT16_LE(data,val)   do { \
                                          gpointer __put_data = data; \
                                          guint16 __put_val = val; \
                                          _GST_PUT (__put_data, 0, 16,  0, __put_val); \
                                          _GST_PUT (__put_data, 1, 16,  8, __put_val); \
                                        } while (0)
#endif /* !GST_HAVE_UNALIGNED_ACCESS */

/**
 * GST_WRITE_UINT8:
 * @data: memory location
 * @num: value to store
 *
 * Store an 8 bit unsigned integer value into the memory buffer.
 */
#define GST_WRITE_UINT8(data, num)      do { \
                                          _GST_PUT (data, 0,  8,  0, num); \
                                        } while (0)

/* Float endianness conversion macros */
#ifndef __GI_SCANNER__

/* FIXME: Remove this once we depend on a GLib version with this */
#ifndef GFLOAT_FROM_LE
/**
 * GFLOAT_SWAP_LE_BE: (skip)
 * @in: input value
 *
 * Swap byte order of a 32-bit floating point value (float).
 *
 * Returns: @in byte-swapped.
 */
static inline gfloat
GFLOAT_SWAP_LE_BE(gfloat in)
{
  union
  {
    guint32 i;
    gfloat f;
  } u;

  u.f = in;
  u.i = GUINT32_SWAP_LE_BE (u.i);
  return u.f;
}

/**
 * GDOUBLE_SWAP_LE_BE: (skip)
 * @in: input value
 *
 * Swap byte order of a 64-bit floating point value (double).
 *
 * Returns: @in byte-swapped.
 */
static inline gdouble
GDOUBLE_SWAP_LE_BE(gdouble in)
{
  union
  {
    guint64 i;
    gdouble d;
  } u;

  u.d = in;
  u.i = GUINT64_SWAP_LE_BE (u.i);
  return u.d;
}

/**
 * GDOUBLE_TO_LE: (skip)
 * @val: value
 *
 * Convert 64-bit floating point value (double) from native byte order into
 * little endian byte order.
 */
/**
 * GDOUBLE_TO_BE: (skip)
 * @val: value
 *
 * Convert 64-bit floating point value (double) from native byte order into
 * big endian byte order.
 */
/**
 * GDOUBLE_FROM_LE: (skip)
 * @val: value
 *
 * Convert 64-bit floating point value (double) from little endian byte order
 * into native byte order.
 */
/**
 * GDOUBLE_FROM_BE: (skip)
 * @val: value
 *
 * Convert 64-bit floating point value (double) from big endian byte order
 * into native byte order.
 */

/**
 * GFLOAT_TO_LE: (skip)
 * @val: value
 *
 * Convert 32-bit floating point value (float) from native byte order into
 * little endian byte order.
 */
/**
 * GFLOAT_TO_BE: (skip)
 * @val: value
 *
 * Convert 32-bit floating point value (float) from native byte order into
 * big endian byte order.
 */
/**
 * GFLOAT_FROM_LE: (skip)
 * @val: value
 *
 * Convert 32-bit floating point value (float) from little endian byte order
 * into native byte order.
 */
/**
 * GFLOAT_FROM_BE: (skip)
 * @val: value
 *
 * Convert 32-bit floating point value (float) from big endian byte order
 * into native byte order.
 */

#if G_BYTE_ORDER == G_LITTLE_ENDIAN
#define GFLOAT_TO_LE(val)    ((gfloat) (val))
#define GFLOAT_TO_BE(val)    (GFLOAT_SWAP_LE_BE (val))
#define GDOUBLE_TO_LE(val)   ((gdouble) (val))
#define GDOUBLE_TO_BE(val)   (GDOUBLE_SWAP_LE_BE (val))

#elif G_BYTE_ORDER == G_BIG_ENDIAN
#define GFLOAT_TO_LE(val)    (GFLOAT_SWAP_LE_BE (val))
#define GFLOAT_TO_BE(val)    ((gfloat) (val))
#define GDOUBLE_TO_LE(val)   (GDOUBLE_SWAP_LE_BE (val))
#define GDOUBLE_TO_BE(val)   ((gdouble) (val))

#else /* !G_LITTLE_ENDIAN && !G_BIG_ENDIAN */
#error unknown ENDIAN type
#endif /* !G_LITTLE_ENDIAN && !G_BIG_ENDIAN */

#define GFLOAT_FROM_LE(val)  (GFLOAT_TO_LE (val))
#define GFLOAT_FROM_BE(val)  (GFLOAT_TO_BE (val))
#define GDOUBLE_FROM_LE(val) (GDOUBLE_TO_LE (val))
#define GDOUBLE_FROM_BE(val) (GDOUBLE_TO_BE (val))

#endif /* !defined(GFLOAT_FROM_LE) */

#endif /* !__GI_SCANNER__ */

/**
 * GST_READ_FLOAT_LE:
 * @data: memory location
 *
 * Read a 32 bit float value in little endian format from the memory buffer.
 *
 * Returns: The floating point value read from @data
 */
static inline gfloat
GST_READ_FLOAT_LE(const guint8 *data)
{
  union
  {
    guint32 i;
    gfloat f;
  } u;

  u.i = GST_READ_UINT32_LE (data);
  return u.f;
}

/**
 * GST_READ_FLOAT_BE:
 * @data: memory location
 *
 * Read a 32 bit float value in big endian format from the memory buffer.
 *
 * Returns: The floating point value read from @data
 */
static inline gfloat
GST_READ_FLOAT_BE(const guint8 *data)
{
  union
  {
    guint32 i;
    gfloat f;
  } u;

  u.i = GST_READ_UINT32_BE (data);
  return u.f;
}

/**
 * GST_READ_DOUBLE_LE:
 * @data: memory location
 *
 * Read a 64 bit double value in little endian format from the memory buffer.
 *
 * Returns: The double-precision floating point value read from @data
 */
static inline gdouble
GST_READ_DOUBLE_LE(const guint8 *data)
{
  union
  {
    guint64 i;
    gdouble d;
  } u;

  u.i = GST_READ_UINT64_LE (data);
  return u.d;
}

/**
 * GST_READ_DOUBLE_BE:
 * @data: memory location
 *
 * Read a 64 bit double value in big endian format from the memory buffer.
 *
 * Returns: The double-precision floating point value read from @data
 */
static inline gdouble
GST_READ_DOUBLE_BE(const guint8 *data)
{
  union
  {
    guint64 i;
    gdouble d;
  } u;

  u.i = GST_READ_UINT64_BE (data);
  return u.d;
}

/**
 * GST_WRITE_FLOAT_LE:
 * @data: memory location
 * @num: value to store
 *
 * Store a 32 bit float value in little endian format into the memory buffer.
 */
static inline void
GST_WRITE_FLOAT_LE(guint8 *data, gfloat num)
{
  union
  {
    guint32 i;
    gfloat f;
  } u;

  u.f = num;
  GST_WRITE_UINT32_LE (data, u.i);
}

/**
 * GST_WRITE_FLOAT_BE:
 * @data: memory location
 * @num: value to store
 *
 * Store a 32 bit float value in big endian format into the memory buffer.
 */
static inline void
GST_WRITE_FLOAT_BE(guint8 *data, gfloat num)
{
  union
  {
    guint32 i;
    gfloat f;
  } u;

  u.f = num;
  GST_WRITE_UINT32_BE (data, u.i);
}

/**
 * GST_WRITE_DOUBLE_LE:
 * @data: memory location
 * @num: value to store
 *
 * Store a 64 bit double value in little endian format into the memory buffer.
 */
static inline void
GST_WRITE_DOUBLE_LE(guint8 *data, gdouble num)
{
  union
  {
    guint64 i;
    gdouble d;
  } u;

  u.d = num;
  GST_WRITE_UINT64_LE (data, u.i);
}

/**
 * GST_WRITE_DOUBLE_BE:
 * @data: memory location
 * @num: value to store
 *
 * Store a 64 bit double value in big endian format into the memory buffer.
 */
static inline void
GST_WRITE_DOUBLE_BE(guint8 *data, gdouble num)
{
  union
  {
    guint64 i;
    gdouble d;
  } u;

  u.d = num;
  GST_WRITE_UINT64_BE (data, u.i);
}

/* Miscellaneous utility macros */

/**
 * GST_ROUND_UP_2:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 2.
 */
#define GST_ROUND_UP_2(num)  (((num)+1)&~1)
/**
 * GST_ROUND_UP_4:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 4.
 */
#define GST_ROUND_UP_4(num)  (((num)+3)&~3)
/**
 * GST_ROUND_UP_8:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 8.
 */
#define GST_ROUND_UP_8(num)  (((num)+7)&~7)
/**
 * GST_ROUND_UP_16:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 16.
 */
#define GST_ROUND_UP_16(num) (((num)+15)&~15)
/**
 * GST_ROUND_UP_32:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 32.
 */
#define GST_ROUND_UP_32(num) (((num)+31)&~31)
/**
 * GST_ROUND_UP_64:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 64.
 */
#define GST_ROUND_UP_64(num) (((num)+63)&~63)
/**
 * GST_ROUND_UP_128:
 * @num: integer value to round up
 *
 * Rounds an integer value up to the next multiple of 128.
 * Since: 1.4
 */
#define GST_ROUND_UP_128(num) (((num)+127)&~127)
/**
 * GST_ROUND_UP_N:
 * @num: integrer value to round up
 * @align: a power of two to round up to
 *
 * Rounds an integer value up to the next multiple of @align. @align MUST be a
 * power of two.
 */
#define GST_ROUND_UP_N(num,align) ((((num) + ((align) - 1)) & ~((align) - 1)))


/**
 * GST_ROUND_DOWN_2:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 2.
 */
#define GST_ROUND_DOWN_2(num)  ((num)&(~1))
/**
 * GST_ROUND_DOWN_4:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 4.
 */
#define GST_ROUND_DOWN_4(num)  ((num)&(~3))
/**
 * GST_ROUND_DOWN_8:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 8.
 */
#define GST_ROUND_DOWN_8(num)  ((num)&(~7))
/**
 * GST_ROUND_DOWN_16:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 16.
 */
#define GST_ROUND_DOWN_16(num) ((num)&(~15))
/**
 * GST_ROUND_DOWN_32:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 32.
 */
#define GST_ROUND_DOWN_32(num) ((num)&(~31))
/**
 * GST_ROUND_DOWN_64:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 64.
 */
#define GST_ROUND_DOWN_64(num) ((num)&(~63))
/**
 * GST_ROUND_DOWN_128:
 * @num: integer value to round down
 *
 * Rounds an integer value down to the next multiple of 128.
 * Since: 1.4
 */
#define GST_ROUND_DOWN_128(num) ((num)&(~127))
/**
 * GST_ROUND_DOWN_N:
 * @num: integrer value to round down
 * @align: a power of two to round down to
 *
 * Rounds an integer value down to the next multiple of @align. @align MUST be a
 * power of two.
 */
#define GST_ROUND_DOWN_N(num,align) (((num) & ~((align) - 1)))


GST_API
void                    gst_object_default_error        (GstObject    * source,
                                                         const GError * error,
                                                         const gchar  * debug);

/* element functions */

GST_API
void                    gst_element_create_all_pads     (GstElement *element);

GST_API
GstPad*                 gst_element_get_compatible_pad  (GstElement *element, GstPad *pad,
                                                         GstCaps *caps);
GST_API
GstPadTemplate*         gst_element_get_compatible_pad_template (GstElement *element, GstPadTemplate *compattempl);

GST_API
const gchar*            gst_element_state_get_name      (GstState state);

GST_API
const gchar *           gst_element_state_change_return_get_name (GstStateChangeReturn state_ret);

GST_API
const gchar *           gst_state_change_get_name       (GstStateChange transition);

GST_API
gboolean                gst_element_link                (GstElement *src, GstElement *dest);

GST_API
gboolean                gst_element_link_many           (GstElement *element_1,
                                                         GstElement *element_2, ...) G_GNUC_NULL_TERMINATED;
GST_API
gboolean                gst_element_link_filtered       (GstElement * src,
                                                         GstElement * dest,
                                                         GstCaps *filter);
GST_API
void                    gst_element_unlink              (GstElement *src, GstElement *dest);

GST_API
void                    gst_element_unlink_many         (GstElement *element_1,
                                                         GstElement *element_2, ...) G_GNUC_NULL_TERMINATED;
GST_API
gboolean                gst_element_link_pads           (GstElement *src, const gchar *srcpadname,
                                                         GstElement *dest, const gchar *destpadname);
GST_API
gboolean                gst_element_link_pads_full      (GstElement *src, const gchar *srcpadname,
                                                         GstElement *dest, const gchar *destpadname,
                                                         GstPadLinkCheck flags);
GST_API
void                    gst_element_unlink_pads         (GstElement *src, const gchar *srcpadname,
                                                         GstElement *dest, const gchar *destpadname);
GST_API
gboolean                gst_element_link_pads_filtered  (GstElement * src, const gchar * srcpadname,
                                                         GstElement * dest, const gchar * destpadname,
                                                         GstCaps *filter);
GST_API
gboolean                gst_element_seek_simple         (GstElement   *element,
                                                         GstFormat     format,
                                                         GstSeekFlags  seek_flags,
                                                         gint64        seek_pos);

/* util elementfactory functions */

GST_API
gboolean gst_element_factory_can_sink_all_caps (GstElementFactory *factory, const GstCaps *caps);

GST_API
gboolean gst_element_factory_can_src_all_caps  (GstElementFactory *factory, const GstCaps *caps);

GST_API
gboolean gst_element_factory_can_sink_any_caps (GstElementFactory *factory, const GstCaps *caps);

GST_API
gboolean gst_element_factory_can_src_any_caps  (GstElementFactory *factory, const GstCaps *caps);

/* util query functions */

GST_API
gboolean                gst_element_query_position      (GstElement *element, GstFormat format, gint64 *cur);

GST_API
gboolean                gst_element_query_duration      (GstElement *element, GstFormat format, gint64 *duration);

GST_API
gboolean                gst_element_query_convert       (GstElement *element, GstFormat src_format, gint64 src_val,
                                                         GstFormat dest_format, gint64 *dest_val);

/* pad functions */

GST_API
void                    gst_pad_use_fixed_caps          (GstPad *pad);

GST_API
GstElement*             gst_pad_get_parent_element      (GstPad *pad);

/* util query functions */

GST_API
gboolean                gst_pad_proxy_query_accept_caps (GstPad *pad, GstQuery *query);

GST_API
gboolean                gst_pad_proxy_query_caps        (GstPad *pad, GstQuery *query);

GST_API
gboolean                gst_pad_query_position          (GstPad *pad, GstFormat format, gint64 *cur);

GST_API
gboolean                gst_pad_query_duration          (GstPad *pad, GstFormat format, gint64 *duration);

GST_API
gboolean                gst_pad_query_convert           (GstPad *pad, GstFormat src_format, gint64 src_val,
                                                         GstFormat dest_format, gint64 *dest_val);
GST_API
GstCaps *               gst_pad_query_caps              (GstPad *pad, GstCaps *filter);

GST_API
gboolean                gst_pad_query_accept_caps       (GstPad *pad, GstCaps *caps);

GST_API
gboolean                gst_pad_link_maybe_ghosting      (GstPad            *src,
                                                          GstPad            *sink);
GST_API
gboolean                gst_pad_link_maybe_ghosting_full (GstPad            *src,
                                                          GstPad            *sink,
                                                          GstPadLinkCheck   flags);
GST_API
gboolean                gst_pad_peer_query_position     (GstPad *pad, GstFormat format, gint64 *cur);

GST_API
gboolean                gst_pad_peer_query_duration     (GstPad *pad, GstFormat format, gint64 *duration);

GST_API
gboolean                gst_pad_peer_query_convert      (GstPad *pad, GstFormat src_format, gint64 src_val,
                                                         GstFormat dest_format, gint64 *dest_val);
GST_API
GstCaps *               gst_pad_peer_query_caps         (GstPad * pad, GstCaps *filter);

GST_API
gboolean                gst_pad_peer_query_accept_caps  (GstPad * pad, GstCaps *caps);

GST_API
gchar *                 gst_pad_create_stream_id               (GstPad * pad, GstElement * parent, const gchar *stream_id) G_GNUC_MALLOC;

GST_API
gchar *                 gst_pad_create_stream_id_printf        (GstPad * pad, GstElement * parent, const gchar *stream_id, ...) G_GNUC_PRINTF (3, 4) G_GNUC_MALLOC;

GST_API
gchar *                 gst_pad_create_stream_id_printf_valist (GstPad * pad, GstElement * parent, const gchar *stream_id, va_list var_args) G_GNUC_PRINTF (3, 0) G_GNUC_MALLOC;

GST_API
gchar *                 gst_pad_get_stream_id           (GstPad * pad);

GST_API
GstStream *             gst_pad_get_stream              (GstPad * pad);

/* bin functions */

GST_API
void                    gst_bin_add_many                (GstBin *bin, GstElement *element_1, ...) G_GNUC_NULL_TERMINATED;

GST_API
void                    gst_bin_remove_many             (GstBin *bin, GstElement *element_1, ...) G_GNUC_NULL_TERMINATED;

GST_API
GstPad *                gst_bin_find_unlinked_pad       (GstBin *bin, GstPadDirection direction);

GST_API
gboolean                gst_bin_sync_children_states    (GstBin *bin);

/* parse utility functions */

GST_API
GstElement *            gst_parse_bin_from_description      (const gchar     * bin_description,
                                                             gboolean          ghost_unlinked_pads,
                                                             GError         ** err);
GST_API
GstElement *            gst_parse_bin_from_description_full (const gchar     * bin_description,
                                                             gboolean          ghost_unlinked_pads,
                                                             GstParseContext * context,
                                                             GstParseFlags     flags,
                                                             GError         ** err);
GST_API
GstClockTime            gst_util_get_timestamp          (void);

/**
 * GstSearchMode:
 * @GST_SEARCH_MODE_EXACT : Only search for exact matches.
 * @GST_SEARCH_MODE_BEFORE: Search for an exact match or the element just before.
 * @GST_SEARCH_MODE_AFTER : Search for an exact match or the element just after.
 *
 * The different search modes.
 */
typedef enum {
  GST_SEARCH_MODE_EXACT = 0,
  GST_SEARCH_MODE_BEFORE,
  GST_SEARCH_MODE_AFTER
} GstSearchMode;

/**
 * GstPluginAPIFlags:
 * @GST_PLUGIN_API_FLAG_IGNORE_ENUM_MEMBERS: Ignore enum members when generating
 *   the plugins cache. This is useful if the members of the enum are generated
 *   dynamically, in order not to expose incorrect documentation to the end user.
 *
 * Since: 1.18
 */
typedef enum {
  GST_PLUGIN_API_FLAG_IGNORE_ENUM_MEMBERS = (1 << 0),
} GstPluginAPIFlags;

GST_API
gpointer      gst_util_array_binary_search      (gpointer array, guint num_elements,
                                                 gsize element_size, GCompareDataFunc search_func,
                                                 GstSearchMode mode, gconstpointer search_data,
                                                 gpointer user_data);

/* fraction operations */

GST_API
gint          gst_util_greatest_common_divisor  (gint a, gint b);

GST_API
gint64        gst_util_greatest_common_divisor_int64 (gint64 a, gint64 b);

GST_API
void          gst_util_fraction_to_double       (gint src_n, gint src_d, gdouble *dest);

GST_API
void          gst_util_double_to_fraction       (gdouble src, gint *dest_n, gint *dest_d);

GST_API
gboolean      gst_util_fraction_multiply        (gint a_n, gint a_d, gint b_n, gint b_d,
                                                 gint *res_n, gint *res_d);
GST_API
gboolean      gst_util_fraction_add             (gint a_n, gint a_d, gint b_n, gint b_d,
                                                 gint *res_n, gint *res_d);
GST_API
gint          gst_util_fraction_compare         (gint a_n, gint a_d, gint b_n, gint b_d);

GST_API
gboolean      gst_calculate_linear_regression   (const GstClockTime * xy,
                                                 GstClockTime * temp, guint n,
                                                 GstClockTime * m_num, GstClockTime * m_denom,
                                                 GstClockTime * b, GstClockTime * xbase,
                                                 gdouble * r_squared);

GST_API
void          gst_type_mark_as_plugin_api       (GType type, GstPluginAPIFlags flags);

GST_API
gboolean      gst_type_is_plugin_api            (GType type, GstPluginAPIFlags *flags);

G_END_DECLS

#endif /* __GST_UTILS_H__ */
