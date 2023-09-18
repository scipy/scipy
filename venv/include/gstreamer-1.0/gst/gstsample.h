/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *
 * gstsample.h: Header for GstSample object
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


#ifndef __GST_SAMPLE_H__
#define __GST_SAMPLE_H__

#include <gst/gstbuffer.h>
#include <gst/gstbufferlist.h>
#include <gst/gstcaps.h>
#include <gst/gstsegment.h>

G_BEGIN_DECLS

GST_API GType _gst_sample_type;

#define GST_TYPE_SAMPLE      (_gst_sample_type)
#define GST_IS_SAMPLE(obj)   (GST_IS_MINI_OBJECT_TYPE(obj, GST_TYPE_SAMPLE))
#define GST_SAMPLE_CAST(obj) ((GstSample *)obj)
#define GST_SAMPLE(obj)      (GST_SAMPLE_CAST(obj))

/**
 * GstSample:
 *
 * The opaque structure of a #GstSample. A sample contains a typed memory
 * block and the associated timing information. It is mainly used to
 * exchange buffers with an application.
 */
typedef struct _GstSample GstSample;

GST_API
GType                gst_sample_get_type      (void);

/* allocation */

GST_API
GstSample *          gst_sample_new           (GstBuffer          *buffer,
                                               GstCaps            *caps,
                                               const GstSegment   *segment,
                                               GstStructure       *info);
GST_API
GstBuffer *          gst_sample_get_buffer    (GstSample *sample);

GST_API
GstCaps *            gst_sample_get_caps      (GstSample *sample);

GST_API
GstSegment *         gst_sample_get_segment   (GstSample *sample);

GST_API
const GstStructure * gst_sample_get_info      (GstSample *sample);

GST_API
GstBufferList *      gst_sample_get_buffer_list (GstSample *sample);

GST_API
void                 gst_sample_set_buffer_list (GstSample *sample, GstBufferList *buffer_list);

GST_API
void                 gst_sample_set_buffer    (GstSample *sample, GstBuffer *buffer);

GST_API
void                 gst_sample_set_caps      (GstSample *sample, GstCaps *caps);

GST_API
void                 gst_sample_set_segment   (GstSample * sample, const GstSegment *segment);

GST_API
gboolean             gst_sample_set_info      (GstSample *sample, GstStructure *info);

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* refcounting */
static inline GstSample *
gst_sample_ref (GstSample * sample)
{
  return GST_SAMPLE_CAST (gst_mini_object_ref (GST_MINI_OBJECT_CAST (
      sample)));
}

static inline void
gst_sample_unref (GstSample * sample)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (sample));
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstSample * gst_sample_ref    (GstSample * sample);

GST_API
void        gst_sample_unref  (GstSample * sample);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/**
 * gst_sample_is_writable:
 * @sample: A #GstSample
 *
 * Tests if you can safely set the buffer and / or buffer list of @sample.
 *
 * Since: 1.16
 */
#define         gst_sample_is_writable(sample)     gst_mini_object_is_writable (GST_MINI_OBJECT_CAST (sample))

/**
 * gst_sample_make_writable:
 * @sample: (transfer full): A #GstSample
 *
 * Returns a writable copy of @sample. If the source sample is
 * already writable, this will simply return the same sample.
 *
 * Use this function to ensure that a sample can be safely modified before
 * making changes to it, for example before calling gst_sample_set_buffer()
 *
 * If the reference count of the source sample @sample is exactly one, the caller
 * is the sole owner and this function will return the sample object unchanged.
 *
 * If there is more than one reference on the object, a copy will be made using
 * gst_sample_copy(). The passed-in @sample will be unreffed in that case, and the
 * caller will now own a reference to the new returned sample object.
 *
 * In short, this function unrefs the sample in the argument and refs the sample
 * that it returns. Don't access the argument after calling this function unless
 * you have an additional reference to it.
 *
 * Returns: (transfer full): a writable sample which may or may not be the
 *     same as @sample
 *
 * Since: 1.16
 */
#define         gst_sample_make_writable(sample)   GST_SAMPLE_CAST (gst_mini_object_make_writable (GST_MINI_OBJECT_CAST (sample)))

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* copy sample */
static inline GstSample *
gst_sample_copy (const GstSample * buf)
{
  return GST_SAMPLE_CAST (gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (buf)));
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstSample *   gst_sample_copy(const GstSample * buf);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/**
 * gst_value_set_sample:
 * @v: a #GValue to receive the data
 * @b: (transfer none): a #GstSample to assign to the GstValue
 *
 * Sets @b as the value of @v.  Caller retains reference to sample.
 */
#define         gst_value_set_sample(v,b)       g_value_set_boxed((v),(b))
/**
 * gst_value_take_sample:
 * @v: a #GValue to receive the data
 * @b: (transfer full): a #GstSample to assign to the GstValue
 *
 * Sets @b as the value of @v.  Caller gives away reference to sample.
 */
#define         gst_value_take_sample(v,b)      g_value_take_boxed(v,(b))
/**
 * gst_value_get_sample:
 * @v: a #GValue to query
 *
 * Receives a #GstSample as the value of @v. Does not return a reference to
 * the sample, so the pointer is only valid for as long as the caller owns
 * a reference to @v.
 *
 * Returns: (transfer none): sample
 */
#define         gst_value_get_sample(v)         GST_SAMPLE_CAST (g_value_get_boxed(v))

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstSample, gst_sample_unref)

G_END_DECLS

#endif /* __GST_SAMPLE_H__ */
