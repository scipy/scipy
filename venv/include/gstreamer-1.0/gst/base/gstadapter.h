/* GStreamer
 * Copyright (C) 2004 Benjamin Otte <otte@gnome.org>
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

#include <gst/gst.h>

#ifndef __GST_ADAPTER_H__
#define __GST_ADAPTER_H__

#include <gst/base/base-prelude.h>

G_BEGIN_DECLS


#define GST_TYPE_ADAPTER \
  (gst_adapter_get_type())
#define GST_ADAPTER(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj), GST_TYPE_ADAPTER, GstAdapter))
#define GST_ADAPTER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass), GST_TYPE_ADAPTER, GstAdapterClass))
#define GST_ADAPTER_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_ADAPTER, GstAdapterClass))
#define GST_IS_ADAPTER(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj), GST_TYPE_ADAPTER))
#define GST_IS_ADAPTER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass), GST_TYPE_ADAPTER))

/**
 * GstAdapter:
 *
 * The opaque #GstAdapter data structure.
 */
typedef struct _GstAdapter GstAdapter;
typedef struct _GstAdapterClass GstAdapterClass;

GST_BASE_API
GType                   gst_adapter_get_type            (void);

GST_BASE_API
GstAdapter *            gst_adapter_new                 (void) G_GNUC_MALLOC;

GST_BASE_API
void                    gst_adapter_clear               (GstAdapter *adapter);

GST_BASE_API
void                    gst_adapter_push                (GstAdapter *adapter, GstBuffer* buf);

GST_BASE_API
gconstpointer           gst_adapter_map                 (GstAdapter *adapter, gsize size);

GST_BASE_API
void                    gst_adapter_unmap               (GstAdapter *adapter);

GST_BASE_API
void                    gst_adapter_copy                (GstAdapter *adapter, gpointer dest,
                                                         gsize offset, gsize size);
GST_BASE_API
GBytes *                gst_adapter_copy_bytes          (GstAdapter *adapter,
                                                         gsize offset, gsize size);
GST_BASE_API
void                    gst_adapter_flush               (GstAdapter *adapter, gsize flush);

GST_BASE_API
gpointer                gst_adapter_take                (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GstBuffer*              gst_adapter_take_buffer         (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GList*                  gst_adapter_take_list           (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GstBuffer *             gst_adapter_take_buffer_fast    (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GstBufferList *         gst_adapter_take_buffer_list    (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GstBuffer*              gst_adapter_get_buffer          (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GList*                  gst_adapter_get_list            (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GstBuffer *             gst_adapter_get_buffer_fast     (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
GstBufferList *         gst_adapter_get_buffer_list     (GstAdapter *adapter, gsize nbytes);

GST_BASE_API
gsize                   gst_adapter_available           (GstAdapter *adapter);

GST_BASE_API
gsize                   gst_adapter_available_fast      (GstAdapter *adapter);

GST_BASE_API
GstClockTime            gst_adapter_prev_pts            (GstAdapter *adapter, guint64 *distance);

GST_BASE_API
GstClockTime            gst_adapter_prev_dts            (GstAdapter *adapter, guint64 *distance);

GST_BASE_API
GstClockTime            gst_adapter_prev_pts_at_offset  (GstAdapter * adapter, gsize offset, guint64 * distance);

GST_BASE_API
GstClockTime            gst_adapter_prev_dts_at_offset  (GstAdapter * adapter, gsize offset, guint64 * distance);

GST_BASE_API
guint64                 gst_adapter_prev_offset         (GstAdapter *adapter, guint64 *distance);

GST_BASE_API
GstClockTime            gst_adapter_pts_at_discont      (GstAdapter *adapter);

GST_BASE_API
GstClockTime            gst_adapter_dts_at_discont      (GstAdapter *adapter);

GST_BASE_API
guint64                 gst_adapter_offset_at_discont   (GstAdapter *adapter);

GST_BASE_API
guint64                 gst_adapter_distance_from_discont (GstAdapter *adapter);

GST_BASE_API
gssize                  gst_adapter_masked_scan_uint32  (GstAdapter * adapter, guint32 mask,
                                                         guint32 pattern, gsize offset, gsize size);
GST_BASE_API
gssize                  gst_adapter_masked_scan_uint32_peek  (GstAdapter * adapter, guint32 mask,
                                                         guint32 pattern, gsize offset, gsize size, guint32 * value);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAdapter, gst_object_unref)

G_END_DECLS

#endif /* __GST_ADAPTER_H__ */
