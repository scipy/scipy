/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *
 * gstbuffer.h: Header for GstBuffer object
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


#ifndef __GST_BUFFER_H__
#define __GST_BUFFER_H__

#include <gst/gstminiobject.h>
#include <gst/gstclock.h>
#include <gst/gstallocator.h>
#include <gst/gstcaps.h>

G_BEGIN_DECLS

GST_API GType _gst_buffer_type;

typedef struct _GstBuffer GstBuffer;
typedef struct _GstBufferPool GstBufferPool;

#include <gst/gstmeta.h>

#define GST_TYPE_BUFFER                         (_gst_buffer_type)
#define GST_IS_BUFFER(obj)                      (GST_IS_MINI_OBJECT_TYPE(obj, GST_TYPE_BUFFER))
#define GST_BUFFER_CAST(obj)                    ((GstBuffer *)(obj))
#define GST_BUFFER(obj)                         (GST_BUFFER_CAST(obj))

/**
 * GST_BUFFER_FLAGS:
 * @buf: a #GstBuffer.
 *
 * Returns a flags word containing #GstBufferFlags flags set on this buffer.
 */
#define GST_BUFFER_FLAGS(buf)                   GST_MINI_OBJECT_FLAGS(buf)
/**
 * GST_BUFFER_FLAG_IS_SET:
 * @buf: a #GstBuffer.
 * @flag: the #GstBufferFlags flag to check.
 *
 * Gives the status of a specific flag on a buffer.
 */
#define GST_BUFFER_FLAG_IS_SET(buf,flag)        GST_MINI_OBJECT_FLAG_IS_SET (buf, flag)
/**
 * GST_BUFFER_FLAG_SET:
 * @buf: a #GstBuffer.
 * @flag: the #GstBufferFlags flag to set.
 *
 * Sets a buffer flag on a buffer.
 */
#define GST_BUFFER_FLAG_SET(buf,flag)           GST_MINI_OBJECT_FLAG_SET (buf, flag)
/**
 * GST_BUFFER_FLAG_UNSET:
 * @buf: a #GstBuffer.
 * @flag: the #GstBufferFlags flag to clear.
 *
 * Clears a buffer flag.
 */
#define GST_BUFFER_FLAG_UNSET(buf,flag)         GST_MINI_OBJECT_FLAG_UNSET (buf, flag)


/**
 * GST_BUFFER_PTS:
 * @buf: a #GstBuffer.:
 *
 * Gets the presentation timestamp (pts) in nanoseconds (as a #GstClockTime)
 * of the data in the buffer. This is the timestamp when the media should be
 * presented to the user.
 *
 * Value will be %GST_CLOCK_TIME_NONE if the pts is unknown.
 */
#define GST_BUFFER_PTS(buf)                     (GST_BUFFER_CAST(buf)->pts)
/**
 * GST_BUFFER_DTS:
 * @buf: a #GstBuffer.
 *
 * Gets the decoding timestamp (dts) in nanoseconds (as a #GstClockTime)
 * of the data in the buffer. This is the timestamp when the media should be
 * decoded or processed otherwise.
 *
 * Value will be %GST_CLOCK_TIME_NONE if the dts is unknown.
 */
#define GST_BUFFER_DTS(buf)                     (GST_BUFFER_CAST(buf)->dts)
/**
 * GST_BUFFER_DTS_OR_PTS:
 * @buf: a #GstBuffer.
 *
 * Returns the buffer decoding timestamp (dts) if valid, else the buffer
 * presentation time (pts)
 *
 * Since: 1.8
 */
#define GST_BUFFER_DTS_OR_PTS(buf)              (GST_BUFFER_DTS_IS_VALID(buf) ? GST_BUFFER_DTS(buf) : GST_BUFFER_PTS (buf))
/**
 * GST_BUFFER_DURATION:
 * @buf: a #GstBuffer.
 *
 * Gets the duration in nanoseconds (as a #GstClockTime) of the data in the buffer.
 *
 * Value will be %GST_CLOCK_TIME_NONE if the duration is unknown.
 */
#define GST_BUFFER_DURATION(buf)                (GST_BUFFER_CAST(buf)->duration)
/**
 * GST_BUFFER_OFFSET:
 * @buf: a #GstBuffer.
 *
 * Gets the offset in the source file of the beginning of this buffer.
 */
#define GST_BUFFER_OFFSET(buf)                  (GST_BUFFER_CAST(buf)->offset)
/**
 * GST_BUFFER_OFFSET_END:
 * @buf: a #GstBuffer.
 *
 * Gets the offset in the source file of the end of this buffer.
 */
#define GST_BUFFER_OFFSET_END(buf)              (GST_BUFFER_CAST(buf)->offset_end)

/**
 * GST_BUFFER_OFFSET_NONE:
 *
 * Constant for no-offset return results.
 */
#define GST_BUFFER_OFFSET_NONE  ((guint64)-1)

/**
 * GST_BUFFER_DURATION_IS_VALID:
 * @buffer: a #GstBuffer
 *
 * Tests if the duration is known.
 */
#define GST_BUFFER_DURATION_IS_VALID(buffer)    (GST_CLOCK_TIME_IS_VALID (GST_BUFFER_DURATION (buffer)))
/**
 * GST_BUFFER_PTS_IS_VALID:
 * @buffer: a #GstBuffer
 *
 * Tests if the pts is known.
 */
#define GST_BUFFER_PTS_IS_VALID(buffer)   (GST_CLOCK_TIME_IS_VALID (GST_BUFFER_PTS (buffer)))
/**
 * GST_BUFFER_DTS_IS_VALID:
 * @buffer: a #GstBuffer
 *
 * Tests if the dts is known.
 */
#define GST_BUFFER_DTS_IS_VALID(buffer)   (GST_CLOCK_TIME_IS_VALID (GST_BUFFER_DTS (buffer)))
/**
 * GST_BUFFER_OFFSET_IS_VALID:
 * @buffer: a #GstBuffer
 *
 * Tests if the start offset is known.
 */
#define GST_BUFFER_OFFSET_IS_VALID(buffer)      (GST_BUFFER_OFFSET (buffer) != GST_BUFFER_OFFSET_NONE)
/**
 * GST_BUFFER_OFFSET_END_IS_VALID:
 * @buffer: a #GstBuffer
 *
 * Tests if the end offset is known.
 */
#define GST_BUFFER_OFFSET_END_IS_VALID(buffer)  (GST_BUFFER_OFFSET_END (buffer) != GST_BUFFER_OFFSET_NONE)
/**
 * GST_BUFFER_IS_DISCONT:
 * @buffer: a #GstBuffer
 *
 * Tests if the buffer marks a discontinuity in the stream.
 */
#define GST_BUFFER_IS_DISCONT(buffer)   (GST_BUFFER_FLAG_IS_SET (buffer, GST_BUFFER_FLAG_DISCONT))

/**
 * GstBufferFlags:
 * @GST_BUFFER_FLAG_LIVE:          the buffer is live data and should be discarded in
 *                                 the PAUSED state.
 * @GST_BUFFER_FLAG_DECODE_ONLY:   the buffer contains data that should be dropped
 *                                 because it will be clipped against the segment
 *                                 boundaries or because it does not contain data
 *                                 that should be shown to the user.
 * @GST_BUFFER_FLAG_DISCONT:       the buffer marks a data discontinuity in the stream.
 *                                 This typically occurs after a seek or a dropped buffer
 *                                 from a live or network source.
 * @GST_BUFFER_FLAG_RESYNC:        the buffer timestamps might have a discontinuity
 *                                 and this buffer is a good point to resynchronize.
 * @GST_BUFFER_FLAG_CORRUPTED:     the buffer data is corrupted.
 * @GST_BUFFER_FLAG_MARKER:        the buffer contains a media specific marker. for
 *                                 video this is the end of a frame boundary, for audio
 *                                 this is the start of a talkspurt. for RTP
 *                                 packets this matches the marker flag in the
 *                                 RTP packet header.
 * @GST_BUFFER_FLAG_HEADER:        the buffer contains header information that is
 *                                 needed to decode the following data.
 * @GST_BUFFER_FLAG_GAP:           the buffer has been created to fill a gap in the
 *                                 stream and contains media neutral data (elements can
 *                                 switch to optimized code path that ignores the buffer
 *                                 content).
 * @GST_BUFFER_FLAG_DROPPABLE:     the buffer can be dropped without breaking the
 *                                 stream, for example to reduce bandwidth.
 * @GST_BUFFER_FLAG_DELTA_UNIT:    this unit cannot be decoded independently.
 * @GST_BUFFER_FLAG_TAG_MEMORY:    this flag is set when memory of the buffer
 *                                 is added/removed
 * @GST_BUFFER_FLAG_LAST:          additional media specific flags can be added starting from
 *                                 this flag.
 *
 * A set of buffer flags used to describe properties of a #GstBuffer.
 */
typedef enum {
  GST_BUFFER_FLAG_LIVE          = (GST_MINI_OBJECT_FLAG_LAST << 0),
  GST_BUFFER_FLAG_DECODE_ONLY   = (GST_MINI_OBJECT_FLAG_LAST << 1),
  GST_BUFFER_FLAG_DISCONT       = (GST_MINI_OBJECT_FLAG_LAST << 2),
  GST_BUFFER_FLAG_RESYNC        = (GST_MINI_OBJECT_FLAG_LAST << 3),
  GST_BUFFER_FLAG_CORRUPTED     = (GST_MINI_OBJECT_FLAG_LAST << 4),
  GST_BUFFER_FLAG_MARKER        = (GST_MINI_OBJECT_FLAG_LAST << 5),
  GST_BUFFER_FLAG_HEADER        = (GST_MINI_OBJECT_FLAG_LAST << 6),
  GST_BUFFER_FLAG_GAP           = (GST_MINI_OBJECT_FLAG_LAST << 7),
  GST_BUFFER_FLAG_DROPPABLE     = (GST_MINI_OBJECT_FLAG_LAST << 8),
  GST_BUFFER_FLAG_DELTA_UNIT    = (GST_MINI_OBJECT_FLAG_LAST << 9),
  GST_BUFFER_FLAG_TAG_MEMORY    = (GST_MINI_OBJECT_FLAG_LAST << 10),

  /**
   * GST_BUFFER_FLAG_SYNC_AFTER:
   *
   * Elements which write to disk or permanent storage should ensure the data
   * is synced after writing the contents of this buffer.
   *
   * Since: 1.6
   */
  GST_BUFFER_FLAG_SYNC_AFTER    = (GST_MINI_OBJECT_FLAG_LAST << 11),

  /**
   * GST_BUFFER_FLAG_NON_DROPPABLE:
   *
   * This buffer is important and should not be dropped.
   *
   * This can be used to mark important buffers, e.g. to flag RTP packets
   * carrying keyframes or codec setup data for RTP Forward Error Correction
   * purposes, or to prevent still video frames from being dropped by elements
   * due to QoS.
   *
   * Since: 1.14
   */
  GST_BUFFER_FLAG_NON_DROPPABLE = (GST_MINI_OBJECT_FLAG_LAST << 12),

  GST_BUFFER_FLAG_LAST          = (GST_MINI_OBJECT_FLAG_LAST << 16)
} GstBufferFlags;

/**
 * GstBuffer:
 * @mini_object: the parent structure
 * @pool: pointer to the pool owner of the buffer
 * @pts: presentation timestamp of the buffer, can be #GST_CLOCK_TIME_NONE when the
 *     pts is not known or relevant. The pts contains the timestamp when the
 *     media should be presented to the user.
 * @dts: decoding timestamp of the buffer, can be #GST_CLOCK_TIME_NONE when the
 *     dts is not known or relevant. The dts contains the timestamp when the
 *     media should be processed.
 * @duration: duration in time of the buffer data, can be #GST_CLOCK_TIME_NONE
 *     when the duration is not known or relevant.
 * @offset: a media specific offset for the buffer data.
 *     For video frames, this is the frame number of this buffer.
 *     For audio samples, this is the offset of the first sample in this buffer.
 *     For file data or compressed data this is the byte offset of the first
 *       byte in this buffer.
 * @offset_end: the last offset contained in this buffer. It has the same
 *     format as @offset.
 *
 * The structure of a #GstBuffer. Use the associated macros to access the public
 * variables.
 */
struct _GstBuffer {
  GstMiniObject          mini_object;

  /*< public >*/ /* with COW */
  GstBufferPool         *pool;

  /* timestamp */
  GstClockTime           pts;
  GstClockTime           dts;
  GstClockTime           duration;

  /* media specific offset */
  guint64                offset;
  guint64                offset_end;
};

GST_API
GType       gst_buffer_get_type            (void);

GST_API
guint       gst_buffer_get_max_memory      (void);

/* allocation */

GST_API
GstBuffer * gst_buffer_new                 (void);

GST_API
GstBuffer * gst_buffer_new_allocate        (GstAllocator * allocator, gsize size,
                                            GstAllocationParams * params);
GST_API
GstBuffer * gst_buffer_new_wrapped_full    (GstMemoryFlags flags, gpointer data, gsize maxsize,
                                            gsize offset, gsize size, gpointer user_data,
                                            GDestroyNotify notify);
GST_API
GstBuffer * gst_buffer_new_wrapped         (gpointer data, gsize size);

GST_API
GstBuffer * gst_buffer_new_wrapped_bytes   (GBytes * bytes);

GST_API
GstBuffer * gst_buffer_new_memdup           (gconstpointer data, gsize size);

/* memory blocks */

GST_API
guint       gst_buffer_n_memory             (GstBuffer *buffer);

GST_API
void        gst_buffer_insert_memory        (GstBuffer *buffer, gint idx, GstMemory *mem);

GST_API
void        gst_buffer_replace_memory_range (GstBuffer *buffer, guint idx, gint length, GstMemory *mem);

GST_API
GstMemory * gst_buffer_peek_memory          (GstBuffer *buffer, guint idx);

GST_API
GstMemory * gst_buffer_get_memory_range     (GstBuffer *buffer, guint idx, gint length);

GST_API
void        gst_buffer_remove_memory_range  (GstBuffer *buffer, guint idx, gint length);

GST_API
void        gst_buffer_prepend_memory       (GstBuffer *buffer, GstMemory *mem);

GST_API
void        gst_buffer_append_memory        (GstBuffer *buffer, GstMemory *mem);

GST_API
void        gst_buffer_replace_memory       (GstBuffer *buffer, guint idx, GstMemory *mem);

GST_API
void        gst_buffer_replace_all_memory   (GstBuffer *buffer, GstMemory *mem);

GST_API
GstMemory * gst_buffer_get_memory           (GstBuffer *buffer, guint idx);

GST_API
GstMemory * gst_buffer_get_all_memory       (GstBuffer *buffer);

GST_API
void        gst_buffer_remove_memory        (GstBuffer *buffer, guint idx);

GST_API
void        gst_buffer_remove_all_memory    (GstBuffer *buffer);

GST_API
gboolean    gst_buffer_find_memory         (GstBuffer *buffer, gsize offset, gsize size,
                                            guint *idx, guint *length, gsize *skip);
GST_API
gboolean    gst_buffer_is_memory_range_writable  (GstBuffer *buffer, guint idx, gint length);

GST_API
gboolean    gst_buffer_is_all_memory_writable    (GstBuffer *buffer);

GST_API
gsize       gst_buffer_fill                (GstBuffer *buffer, gsize offset,
                                            gconstpointer src, gsize size);
GST_API
gsize       gst_buffer_extract             (GstBuffer *buffer, gsize offset,
                                            gpointer dest, gsize size);
GST_API
gint        gst_buffer_memcmp              (GstBuffer *buffer, gsize offset,
                                            gconstpointer mem, gsize size);
GST_API
gsize       gst_buffer_memset              (GstBuffer *buffer, gsize offset,
                                            guint8 val, gsize size);
GST_API
gsize       gst_buffer_get_sizes_range     (GstBuffer *buffer, guint idx, gint length,
                                            gsize *offset, gsize *maxsize);
GST_API
gboolean    gst_buffer_resize_range        (GstBuffer *buffer, guint idx, gint length,
                                            gssize offset, gssize size);
GST_API
gsize       gst_buffer_get_sizes           (GstBuffer *buffer, gsize *offset, gsize *maxsize);

GST_API
gsize       gst_buffer_get_size            (GstBuffer *buffer);

GST_API
void        gst_buffer_resize              (GstBuffer *buffer, gssize offset, gssize size);

GST_API
void        gst_buffer_set_size            (GstBuffer *buffer, gssize size);

GST_API
gboolean    gst_buffer_map_range           (GstBuffer *buffer, guint idx, gint length,
                                            GstMapInfo *info, GstMapFlags flags);
GST_API
gboolean    gst_buffer_map                 (GstBuffer *buffer, GstMapInfo *info, GstMapFlags flags);

GST_API
void        gst_buffer_unmap               (GstBuffer *buffer, GstMapInfo *info);

GST_API
void        gst_buffer_extract_dup         (GstBuffer *buffer, gsize offset,
                                            gsize size, gpointer *dest,
                                            gsize *dest_size);
GST_API
GstBufferFlags gst_buffer_get_flags        (GstBuffer * buffer);

GST_API
gboolean       gst_buffer_has_flags        (GstBuffer * buffer, GstBufferFlags flags);

GST_API
gboolean       gst_buffer_set_flags        (GstBuffer * buffer, GstBufferFlags flags);

GST_API
gboolean       gst_buffer_unset_flags      (GstBuffer * buffer, GstBufferFlags flags);


#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* refcounting */
static inline GstBuffer *
gst_buffer_ref (GstBuffer * buf)
{
  return (GstBuffer *) gst_mini_object_ref (GST_MINI_OBJECT_CAST (buf));
}

static inline void
gst_buffer_unref (GstBuffer * buf)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (buf));
}

static inline void
gst_clear_buffer (GstBuffer ** buf_ptr)
{
  gst_clear_mini_object ((GstMiniObject **) buf_ptr);
}

/* copy buffer */
static inline GstBuffer *
gst_buffer_copy (const GstBuffer * buf)
{
  return GST_BUFFER (gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (buf)));
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstBuffer * gst_buffer_ref       (GstBuffer * buf);

GST_API
void        gst_buffer_unref     (GstBuffer * buf);

GST_API
void        gst_clear_buffer     (GstBuffer ** buf_ptr);

GST_API
GstBuffer * gst_buffer_copy      (const GstBuffer * buf);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

GST_API
GstBuffer * gst_buffer_copy_deep (const GstBuffer * buf);

/**
 * GstBufferCopyFlags:
 * @GST_BUFFER_COPY_NONE: copy nothing
 * @GST_BUFFER_COPY_FLAGS: flag indicating that buffer flags should be copied
 * @GST_BUFFER_COPY_TIMESTAMPS: flag indicating that buffer pts, dts,
 *   duration, offset and offset_end should be copied
 * @GST_BUFFER_COPY_MEMORY: flag indicating that buffer memory should be reffed
 *   and appended to already existing memory. Unless the memory is marked as
 *   NO_SHARE, no actual copy of the memory is made but it is simply reffed.
 *   Add @GST_BUFFER_COPY_DEEP to force a real copy.
 * @GST_BUFFER_COPY_MERGE: flag indicating that buffer memory should be
 *   merged
 * @GST_BUFFER_COPY_META: flag indicating that buffer meta should be
 *   copied
 *
 * A set of flags that can be provided to the gst_buffer_copy_into()
 * function to specify which items should be copied.
 */
typedef enum {
  GST_BUFFER_COPY_NONE           = 0,
  GST_BUFFER_COPY_FLAGS          = (1 << 0),
  GST_BUFFER_COPY_TIMESTAMPS     = (1 << 1),
  GST_BUFFER_COPY_META           = (1 << 2),
  GST_BUFFER_COPY_MEMORY         = (1 << 3),
  GST_BUFFER_COPY_MERGE          = (1 << 4),

  /**
   * GST_BUFFER_COPY_DEEP:
   *
   * flag indicating that memory should always be copied instead of reffed
   *
   * Since: 1.2
   */
  GST_BUFFER_COPY_DEEP           = (1 << 5)
} GstBufferCopyFlags;

/**
 * GST_BUFFER_COPY_METADATA: (value 7) (type GstBufferCopyFlags)
 *
 * Combination of all possible metadata fields that can be copied with
 * gst_buffer_copy_into().
 */
#define GST_BUFFER_COPY_METADATA       ((GstBufferCopyFlags) (GST_BUFFER_COPY_FLAGS |\
                                          GST_BUFFER_COPY_TIMESTAMPS | GST_BUFFER_COPY_META))

/**
 * GST_BUFFER_COPY_ALL: (value 15) (type GstBufferCopyFlags)
 *
 * Combination of all possible fields that can be copied with
 * gst_buffer_copy_into().
 */
#define GST_BUFFER_COPY_ALL  ((GstBufferCopyFlags)(GST_BUFFER_COPY_METADATA | GST_BUFFER_COPY_MEMORY))

/* copies memory or metadata into newly allocated buffer */

GST_API
gboolean        gst_buffer_copy_into            (GstBuffer *dest, GstBuffer *src,
                                                 GstBufferCopyFlags flags,
                                                 gsize offset, gsize size);

/**
 * gst_buffer_is_writable:
 * @buf: a #GstBuffer
 *
 * Tests if you can safely write to a buffer's metadata or its memory array.
 * It is only safe to change buffer metadata when the current reference is
 * writable, i.e. nobody can see the modifications you will make.
 */
#define         gst_buffer_is_writable(buf)     gst_mini_object_is_writable (GST_MINI_OBJECT_CAST (buf))
/**
 * gst_buffer_make_writable:
 * @buf: (transfer full): a #GstBuffer
 *
 * Returns a writable copy of @buf. If the source buffer is
 * already writable, this will simply return the same buffer.
 *
 * Use this function to ensure that a buffer can be safely modified before
 * making changes to it, including changing the metadata such as PTS/DTS.
 *
 * If the reference count of the source buffer @buf is exactly one, the caller
 * is the sole owner and this function will return the buffer object unchanged.
 *
 * If there is more than one reference on the object, a copy will be made using
 * gst_buffer_copy(). The passed-in @buf will be unreffed in that case, and the
 * caller will now own a reference to the new returned buffer object. Note
 * that this just copies the buffer structure itself, the underlying memory is
 * not copied if it can be shared amongst multiple buffers.
 *
 * In short, this function unrefs the buf in the argument and refs the buffer
 * that it returns. Don't access the argument after calling this function unless
 * you have an additional reference to it.
 *
 * Returns: (transfer full) (nullable): a writable buffer (which may or may not be the
 *     same as @buf) or %NULL if copying is required but not possible.
 */
#define         gst_buffer_make_writable(buf)   GST_BUFFER_CAST (gst_mini_object_make_writable (GST_MINI_OBJECT_CAST (buf)))

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
static inline gboolean
gst_buffer_replace (GstBuffer **obuf, GstBuffer *nbuf)
{
  return gst_mini_object_replace ((GstMiniObject **) obuf, (GstMiniObject *) nbuf);
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
gboolean        gst_buffer_replace              (GstBuffer ** obuf,
                                                 GstBuffer * nbuf);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/* creating a region */

GST_API
GstBuffer*      gst_buffer_copy_region          (GstBuffer *parent, GstBufferCopyFlags flags,
                                                 gsize offset, gsize size);

/* append two buffers */

GST_API
GstBuffer*      gst_buffer_append_region        (GstBuffer *buf1, GstBuffer *buf2,
                                                 gssize offset, gssize size);
GST_API
GstBuffer*      gst_buffer_append               (GstBuffer *buf1, GstBuffer *buf2);

/* metadata */
#include <gst/gstmeta.h>

/**
 * GstBufferForeachMetaFunc:
 * @buffer: a #GstBuffer
 * @meta: (out) (nullable): a pointer to a #GstMeta
 * @user_data: user data passed to gst_buffer_foreach_meta()
 *
 * A function that will be called from gst_buffer_foreach_meta(). The @meta
 * field will point to a the reference of the meta.
 *
 * @buffer should not be modified from this callback.
 *
 * When this function returns %TRUE, the next meta will be
 * returned. When %FALSE is returned, gst_buffer_foreach_meta() will return.
 *
 * When @meta is set to %NULL, the item will be removed from the buffer.
 *
 * Returns: %FALSE when gst_buffer_foreach_meta() should stop
 */
typedef gboolean (*GstBufferForeachMetaFunc)    (GstBuffer *buffer, GstMeta **meta,
                                                 gpointer user_data);

GST_API
GstMeta *       gst_buffer_get_meta             (GstBuffer *buffer, GType api);

GST_API
guint           gst_buffer_get_n_meta           (GstBuffer *buffer, GType api_type);

GST_API
GstMeta *       gst_buffer_add_meta             (GstBuffer *buffer, const GstMetaInfo *info,
                                                 gpointer params);
GST_API
gboolean        gst_buffer_remove_meta          (GstBuffer *buffer, GstMeta *meta);

GST_API
GstMeta *       gst_buffer_iterate_meta         (GstBuffer *buffer, gpointer *state);

GST_API
GstMeta *       gst_buffer_iterate_meta_filtered (GstBuffer * buffer,
                                                  gpointer  * state,
                                                  GType       meta_api_type);
GST_API
gboolean        gst_buffer_foreach_meta         (GstBuffer *buffer,
                                                 GstBufferForeachMetaFunc func,
                                                 gpointer user_data);

GST_API
GstCustomMeta * gst_buffer_add_custom_meta      (GstBuffer *buffer,
                                                 const gchar *name);

GST_API
GstCustomMeta * gst_buffer_get_custom_meta      (GstBuffer *buffer,
                                                 const gchar *name);

/**
 * gst_value_set_buffer:
 * @v: a #GValue to receive the data
 * @b: (transfer none): a #GstBuffer to assign to the GstValue
 *
 * Sets @b as the value of @v.  Caller retains reference to buffer.
 */
#define         gst_value_set_buffer(v,b)       g_value_set_boxed((v),(b))
/**
 * gst_value_take_buffer:
 * @v: a #GValue to receive the data
 * @b: (transfer full): a #GstBuffer to assign to the GstValue
 *
 * Sets @b as the value of @v.  Caller gives away reference to buffer.
 */
#define         gst_value_take_buffer(v,b)      g_value_take_boxed(v,(b))
/**
 * gst_value_get_buffer:
 * @v: a #GValue to query
 *
 * Receives a #GstBuffer as the value of @v. Does not return a reference to
 * the buffer, so the pointer is only valid for as long as the caller owns
 * a reference to @v.
 *
 * Returns: (transfer none): buffer
 */
#define         gst_value_get_buffer(v)         GST_BUFFER_CAST (g_value_get_boxed(v))

typedef struct _GstParentBufferMeta GstParentBufferMeta;

/**
 * GstParentBufferMeta:
 * @parent: the parent #GstMeta structure
 * @buffer: the #GstBuffer on which a reference is being held.
 *
 * The #GstParentBufferMeta is a #GstMeta which can be attached to a #GstBuffer
 * to hold a reference to another buffer that is only released when the child
 * #GstBuffer is released.
 *
 * Typically, #GstParentBufferMeta is used when the child buffer is directly
 * using the #GstMemory of the parent buffer, and wants to prevent the parent
 * buffer from being returned to a buffer pool until the #GstMemory is available
 * for re-use.
 *
 * Since: 1.6
 */
struct _GstParentBufferMeta
{
  GstMeta parent;

  /*< public >*/
  GstBuffer *buffer;
};

GST_API
GType gst_parent_buffer_meta_api_get_type (void);
#ifndef GST_DISABLE_DEPRECATED
#define GST_TYPE_PARENT_BUFFER_META_API_TYPE GST_PARENT_BUFFER_META_API_TYPE
#endif
#define GST_PARENT_BUFFER_META_API_TYPE (gst_parent_buffer_meta_api_get_type())

/**
 * gst_buffer_get_parent_buffer_meta:
 * @b: a #GstBuffer
 *
 * Finds and returns a #GstParentBufferMeta if one exists on the
 * buffer
 */
#define gst_buffer_get_parent_buffer_meta(b) \
  ((GstParentBufferMeta*)gst_buffer_get_meta((b),GST_PARENT_BUFFER_META_API_TYPE))

GST_API
const GstMetaInfo *gst_parent_buffer_meta_get_info (void);
#define GST_PARENT_BUFFER_META_INFO (gst_parent_buffer_meta_get_info())

/* implementation */

GST_API
GstParentBufferMeta *gst_buffer_add_parent_buffer_meta (GstBuffer *buffer,
    GstBuffer *ref);

typedef struct _GstReferenceTimestampMeta GstReferenceTimestampMeta;

/**
 * GstReferenceTimestampMeta:
 * @parent: the parent #GstMeta structure
 * @reference: identifier for the timestamp reference.
 * @timestamp: timestamp
 * @duration: duration, or %GST_CLOCK_TIME_NONE
 *
 * #GstReferenceTimestampMeta can be used to attach alternative timestamps and
 * possibly durations to a #GstBuffer. These are generally not according to
 * the pipeline clock and could be e.g. the NTP timestamp when the media was
 * captured.
 *
 * The reference is stored as a #GstCaps in @reference. Examples of valid
 * references would be
 *
 *  * `timestamp/x-drivername-stream`: for timestamps that are locally
 *    generated by some driver named `drivername` when generating the stream,
 *    e.g. based on a frame counter
 *  * `timestamp/x-ntp, host=pool.ntp.org, port=123`: for timestamps based on a
 *    specific NTP server. Note that the host/port parameters might not always
 *    be given.
 *  * `timestamp/x-ptp, version=IEEE1588-2008, domain=1`: for timestamps based
 *    on a given PTP clock.
 *  * `timestamp/x-unix`: for timestamps based on the UNIX epoch according to
 *    the local clock.
 *
 * Since: 1.14
 */
struct _GstReferenceTimestampMeta
{
  GstMeta parent;

  /*< public >*/
  GstCaps *reference;
  GstClockTime timestamp, duration;
};

GST_API
GType gst_reference_timestamp_meta_api_get_type (void);
#define GST_REFERENCE_TIMESTAMP_META_API_TYPE (gst_reference_timestamp_meta_api_get_type())

GST_API
const GstMetaInfo *gst_reference_timestamp_meta_get_info (void);
#define GST_REFERENCE_TIMESTAMP_META_INFO (gst_reference_timestamp_meta_get_info())

/* implementation */

GST_API
GstReferenceTimestampMeta * gst_buffer_add_reference_timestamp_meta (GstBuffer  * buffer,
                                                                     GstCaps    * reference,
                                                                     GstClockTime timestamp,
                                                                     GstClockTime duration);

GST_API
GstReferenceTimestampMeta * gst_buffer_get_reference_timestamp_meta (GstBuffer * buffer,
                                                                     GstCaps   * reference);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBuffer, gst_buffer_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBufferPool, gst_object_unref)

/**
 * GstBufferMapInfo: (skip):
 *
 * Alias for #GstMapInfo to be used with g_auto():
 * ```c
 * void my_func(GstBuffer *buf)
 * {
 *   g_auto(GstBufferMapInfo) map = GST_MAP_INFO_INIT;
 *   if (!gst_buffer_map(buf, &map, GST_MAP_READWRITE))
 *     return;
 *   ...
 *   // No need to call gst_buffer_unmap()
 * }
 * ```
 *
 * #GstMapInfo cannot be used with g_auto() because it is ambiguous whether it
 * needs to be unmapped using gst_buffer_unmap() or gst_memory_unmap().
 *
 * See also #GstMemoryMapInfo.
 *
 * Since: 1.22
 */
typedef GstMapInfo GstBufferMapInfo;

static inline void _gst_buffer_map_info_clear(GstBufferMapInfo *info)
{
  /* we need to check for NULL, it is possible that we tried to map a buffer
   * without memory and we should be able to unmap that fine */
  if (G_LIKELY (info->memory)) {
    gst_memory_unmap (info->memory, info);
    gst_memory_unref (info->memory);
  }
}

G_DEFINE_AUTO_CLEANUP_CLEAR_FUNC(GstBufferMapInfo, _gst_buffer_map_info_clear)

G_END_DECLS

#endif /* __GST_BUFFER_H__ */
