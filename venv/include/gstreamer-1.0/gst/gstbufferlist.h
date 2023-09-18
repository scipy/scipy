/* GStreamer
 * Copyright (C) 2009 Axis Communications <dev-gstreamer at axis dot com>
 * @author Jonas Holmberg <jonas dot holmberg at axis dot com>
 *
 * gstbufferlist.h: Header for GstBufferList object
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

#ifndef __GST_BUFFER_LIST_H__
#define __GST_BUFFER_LIST_H__

#include <gst/gstbuffer.h>

G_BEGIN_DECLS

GST_API GType _gst_buffer_list_type;

#define GST_TYPE_BUFFER_LIST      (_gst_buffer_list_type)
#define GST_IS_BUFFER_LIST(obj)   (GST_IS_MINI_OBJECT_TYPE(obj, GST_TYPE_BUFFER_LIST))
#define GST_BUFFER_LIST_CAST(obj) ((GstBufferList *)obj)
#define GST_BUFFER_LIST(obj)      (GST_BUFFER_LIST_CAST(obj))

typedef struct _GstBufferList GstBufferList;

/**
 * GstBufferListFunc:
 * @buffer: (out) (nullable): pointer to the buffer
 * @idx: the index of @buffer
 * @user_data: user data passed to gst_buffer_list_foreach()
 *
 * A function that will be called from gst_buffer_list_foreach(). The @buffer
 * field will point to a the reference of the buffer at @idx.
 *
 * When this function returns %TRUE, the next buffer will be
 * returned. When %FALSE is returned, gst_buffer_list_foreach() will return.
 *
 * When @buffer is set to %NULL, the item will be removed from the bufferlist.
 * When @buffer has been made writable, the new buffer reference can be assigned
 * to @buffer. This function is responsible for unreffing the old buffer when
 * removing or modifying.
 *
 * Returns: %FALSE when gst_buffer_list_foreach() should stop
 */
typedef gboolean   (*GstBufferListFunc)   (GstBuffer **buffer, guint idx,
                                           gpointer user_data);

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* refcounting */
static inline GstBufferList *
gst_buffer_list_ref (GstBufferList * list)
{
  return GST_BUFFER_LIST_CAST (gst_mini_object_ref (GST_MINI_OBJECT_CAST (
      list)));
}

static inline void
gst_buffer_list_unref(GstBufferList* list)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (list));
}

static inline void
gst_clear_buffer_list (GstBufferList ** list_ptr)
{
  gst_clear_mini_object ((GstMiniObject **) list_ptr);
}

/* copy */
static inline GstBufferList *
gst_buffer_list_copy (const GstBufferList * list)
{
  return GST_BUFFER_LIST_CAST (gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (list)));
}

static inline gboolean
gst_buffer_list_replace (GstBufferList **old_list, GstBufferList *new_list)
{
  return gst_mini_object_replace ((GstMiniObject **) old_list,
      (GstMiniObject *) new_list);
}

static inline gboolean
gst_buffer_list_take (GstBufferList **old_list, GstBufferList *new_list)
{
  return gst_mini_object_take ((GstMiniObject **) old_list,
      (GstMiniObject *) new_list);
}
#else  /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstBufferList * gst_buffer_list_ref     (GstBufferList * list);

GST_API
void            gst_buffer_list_unref   (GstBufferList * list);

GST_API
void            gst_clear_buffer_list   (GstBufferList ** list_ptr);

GST_API
GstBufferList * gst_buffer_list_copy    (const GstBufferList * list);

GST_API
gboolean        gst_buffer_list_replace (GstBufferList ** old_list,
                                         GstBufferList * new_list);

GST_API
gboolean        gst_buffer_list_take    (GstBufferList ** old_list,
                                         GstBufferList * new_list);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/**
 * gst_buffer_list_is_writable:
 * @list: a #GstBufferList
 *
 * Tests if you can safely add buffers and groups into a buffer list.
 */
#define gst_buffer_list_is_writable(list) gst_mini_object_is_writable (GST_MINI_OBJECT_CAST (list))

/**
 * gst_buffer_list_make_writable:
 * @list: (transfer full): a #GstBufferList
 *
 * Makes a writable buffer list from the given buffer list. If the source buffer
 * list is already writable, this will simply return the same buffer list. A
 * copy will otherwise be made using gst_buffer_list_copy().
 *
 * Returns: (transfer full): a writable list, which may or may not be the
 *     same as @list
 */
#define gst_buffer_list_make_writable(list) GST_BUFFER_LIST_CAST (gst_mini_object_make_writable (GST_MINI_OBJECT_CAST (list)))

GST_API
GType                    gst_buffer_list_get_type              (void);

/* allocation */

GST_API
GstBufferList *          gst_buffer_list_new                   (void) G_GNUC_MALLOC;

GST_API
GstBufferList *          gst_buffer_list_new_sized             (guint size) G_GNUC_MALLOC;

GST_API
guint                    gst_buffer_list_length                (GstBufferList *list);

GST_API
GstBuffer *              gst_buffer_list_get                   (GstBufferList *list, guint idx);

GST_API
GstBuffer *              gst_buffer_list_get_writable          (GstBufferList *list, guint idx);

GST_API
void                     gst_buffer_list_insert                (GstBufferList *list, gint idx, GstBuffer *buffer);

GST_API
void                     gst_buffer_list_remove                (GstBufferList *list, guint idx, guint length);

GST_API
gboolean                 gst_buffer_list_foreach               (GstBufferList *list,
                                                                GstBufferListFunc func,
								gpointer user_data);
GST_API
GstBufferList *          gst_buffer_list_copy_deep             (const GstBufferList * list);

GST_API
gsize                    gst_buffer_list_calculate_size        (GstBufferList * list);

#define gst_buffer_list_add(l,b) gst_buffer_list_insert((l),-1,(b));

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBufferList, gst_buffer_list_unref)

G_END_DECLS

#endif /* __GST_BUFFER_LIST_H__ */
