/* GStreamer
 * Copyright (C) 2004 Wim Taymans <wim@fluendo.com>
 * Copyright (C) 2011 Sebastian Dröge <sebastian.droege@collabora.co.uk>
 *
 * gstiterator.h: Header for GstIterator
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

#ifndef __GST_ITERATOR_H__
#define __GST_ITERATOR_H__

#include <glib-object.h> /* for GValue in the fold */
#include <gst/gstconfig.h>

G_BEGIN_DECLS

#define GST_TYPE_ITERATOR (gst_iterator_get_type ())

/**
 * GstIteratorResult:
 * @GST_ITERATOR_DONE:   No more items in the iterator
 * @GST_ITERATOR_OK:     An item was retrieved
 * @GST_ITERATOR_RESYNC: Datastructure changed while iterating
 * @GST_ITERATOR_ERROR:  An error happened
 *
 * The result of gst_iterator_next().
 */
typedef enum {
  GST_ITERATOR_DONE     = 0,
  GST_ITERATOR_OK       = 1,
  GST_ITERATOR_RESYNC   = 2,
  GST_ITERATOR_ERROR    = 3
} GstIteratorResult;

typedef struct _GstIterator GstIterator;

/**
 * GstIteratorItem:
 * @GST_ITERATOR_ITEM_SKIP:  Skip this item
 * @GST_ITERATOR_ITEM_PASS:  Return item
 * @GST_ITERATOR_ITEM_END:   Stop after this item.
 *
 * The result of a #GstIteratorItemFunction.
 */
typedef enum {
  GST_ITERATOR_ITEM_SKIP        = 0,
  GST_ITERATOR_ITEM_PASS        = 1,
  GST_ITERATOR_ITEM_END         = 2
} GstIteratorItem;

/**
 * GstIteratorCopyFunction:
 * @it: The original iterator
 * @copy: The copied iterator
 *
 * This function will be called when creating a copy of @it and should
 * create a copy of all custom iterator fields or increase their
 * reference counts.
 */
typedef void              (*GstIteratorCopyFunction) (const GstIterator *it, GstIterator *copy);

/**
 * GstIteratorItemFunction:
 * @it: the iterator
 * @item: the item being retrieved.
 *
 * The function that will be called after the next item of the iterator
 * has been retrieved. This function can be used to skip items or stop
 * the iterator.
 *
 * The function will be called with the iterator lock held.
 *
 * Returns: the result of the operation.
 */
typedef GstIteratorItem   (*GstIteratorItemFunction)    (GstIterator *it, const GValue * item);

/**
 * GstIteratorNextFunction:
 * @it: the iterator
 * @result: a pointer to hold the next item
 *
 * The function that will be called when the next element of the iterator
 * should be retrieved.
 *
 * Implementors of a #GstIterator should implement this
 * function and pass it to the constructor of the custom iterator.
 * The function will be called with the iterator lock held.
 *
 * Returns: the result of the operation.
 */
typedef GstIteratorResult (*GstIteratorNextFunction)    (GstIterator *it, GValue *result);
/**
 * GstIteratorResyncFunction:
 * @it: the iterator
 *
 * This function will be called whenever a concurrent update happened
 * to the iterated datastructure. The implementor of the iterator should
 * restart the iterator from the beginning and clean up any state it might
 * have.
 *
 * Implementors of a #GstIterator should implement this
 * function and pass it to the constructor of the custom iterator.
 * The function will be called with the iterator lock held.
 */
typedef void              (*GstIteratorResyncFunction)  (GstIterator *it);
/**
 * GstIteratorFreeFunction:
 * @it: the iterator
 *
 * This function will be called when the iterator is freed.
 *
 * Implementors of a #GstIterator should implement this
 * function and pass it to the constructor of the custom iterator.
 * The function will be called with the iterator lock held.
 */
typedef void              (*GstIteratorFreeFunction)    (GstIterator *it);

/**
 * GstIteratorForeachFunction:
 * @item: The item
 * @user_data: User data
 *
 * A function that is called by gst_iterator_foreach() for every element.
 */
typedef void         (*GstIteratorForeachFunction)     (const GValue * item, gpointer user_data);

/**
 * GstIteratorFoldFunction:
 * @item: the item to fold
 * @ret: a #GValue collecting the result
 * @user_data: data passed to gst_iterator_fold()
 *
 * A function to be passed to gst_iterator_fold().
 *
 * Returns: %TRUE if the fold should continue, %FALSE if it should stop.
 */
typedef gboolean          (*GstIteratorFoldFunction)    (const GValue * item, GValue * ret, gpointer user_data);

/**
 * GST_ITERATOR:
 * @it: the #GstIterator value
 *
 * Macro to cast to a #GstIterator
 */
#define GST_ITERATOR(it)                ((GstIterator*)(it))
/**
 * GST_ITERATOR_LOCK:
 * @it: the #GstIterator to get the lock of
 *
 * Macro to get the lock protecting the datastructure being iterated.
 */
#define GST_ITERATOR_LOCK(it)           (GST_ITERATOR(it)->lock)
/**
 * GST_ITERATOR_COOKIE:
 * @it: the #GstIterator to get the cookie of
 *
 * Macro to get the cookie of a #GstIterator. The cookie of the
 * iterator is the value of the master cookie when the iterator
 * was created.
 * Whenever the iterator is iterated, the value is compared to the
 * value of the master cookie. If they are different, a concurrent
 * modification happened to the iterator and a resync is needed.
 */
#define GST_ITERATOR_COOKIE(it)         (GST_ITERATOR(it)->cookie)
/**
 * GST_ITERATOR_ORIG_COOKIE:
 * @it: the #GstIterator to get the master cookie of
 *
 * Macro to get a pointer to where the master cookie is stored. The
 * master cookie protects the structure being iterated and gets updated
 * whenever the datastructure changes.
 */
#define GST_ITERATOR_ORIG_COOKIE(it)    (GST_ITERATOR(it)->master_cookie)

/**
 * GstIterator:
 * @copy: The function to copy the iterator
 * @next: The function to get the next item in the iterator
 * @item: The function to be called for each item retrieved
 * @resync: The function to call when a resync is needed.
 * @free: The function to call when the iterator is freed
 * @pushed: The iterator that is currently pushed with gst_iterator_push()
 * @type: The type of the object that this iterator will return
 * @lock: The lock protecting the data structure and the cookie.
 * @cookie: The cookie; the value of the master_cookie when this iterator was
 *          created.
 * @master_cookie: A pointer to the master cookie.
 * @size: the size of the iterator
 *
 * #GstIterator base structure. The values of this structure are
 * protected for subclasses, use the methods to use the #GstIterator.
 */
struct _GstIterator {
  /*< protected >*/
  GstIteratorCopyFunction copy;
  GstIteratorNextFunction next;
  GstIteratorItemFunction item;
  GstIteratorResyncFunction resync;
  GstIteratorFreeFunction free;

  GstIterator *pushed;          /* pushed iterator */

  GType     type;
  GMutex   *lock;
  guint32   cookie;             /* cookie of the iterator */
  guint32  *master_cookie;      /* pointer to guint32 holding the cookie when this
                                   iterator was created */
  guint     size;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType                   gst_iterator_get_type           (void);

/* creating iterators */

GST_API
GstIterator*            gst_iterator_new                (guint size,
                                                         GType type,
                                                         GMutex *lock,
                                                         guint32 *master_cookie,
                                                         GstIteratorCopyFunction copy,
                                                         GstIteratorNextFunction next,
                                                         GstIteratorItemFunction item,
                                                         GstIteratorResyncFunction resync,
                                                         GstIteratorFreeFunction free) G_GNUC_MALLOC;
GST_API
GstIterator*            gst_iterator_new_list           (GType type,
                                                         GMutex *lock,
                                                         guint32 *master_cookie,
                                                         GList **list,
                                                         GObject * owner,
                                                         GstIteratorItemFunction item) G_GNUC_MALLOC;
GST_API
GstIterator*            gst_iterator_new_single         (GType type,
                                                         const GValue * object) G_GNUC_MALLOC;
GST_API
GstIterator*            gst_iterator_copy               (const GstIterator *it) G_GNUC_MALLOC;

/* using iterators */

GST_API
GstIteratorResult       gst_iterator_next               (GstIterator *it, GValue * elem);

GST_API
void                    gst_iterator_resync             (GstIterator *it);

GST_API
void                    gst_iterator_free               (GstIterator *it);

GST_API
void                    gst_iterator_push               (GstIterator *it, GstIterator *other);

/* higher-order functions that operate on iterators */

GST_API
GstIterator*            gst_iterator_filter             (GstIterator *it, GCompareFunc func,
                                                         const GValue * user_data) G_GNUC_MALLOC;
GST_API
GstIteratorResult       gst_iterator_fold               (GstIterator *it,
                                                         GstIteratorFoldFunction func,
                                                         GValue *ret, gpointer user_data);
GST_API
GstIteratorResult       gst_iterator_foreach            (GstIterator *it,
                                                         GstIteratorForeachFunction func, gpointer user_data);
GST_API
gboolean                gst_iterator_find_custom        (GstIterator *it, GCompareFunc func,
                                                         GValue *elem, gpointer user_data);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstIterator, gst_iterator_free)

G_END_DECLS

#endif /* __GST_ITERATOR_H__ */
