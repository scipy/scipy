/* GStreamer
 * Copyright (C) 2005 David Schleef <ds@schleef.org>
 *
 * gstminiobject.h: Header for GstMiniObject
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


#ifndef __GST_MINI_OBJECT_H__
#define __GST_MINI_OBJECT_H__

#include <gst/gstconfig.h>

#include <glib-object.h>

G_BEGIN_DECLS

GST_API GType _gst_mini_object_type;

#define GST_TYPE_MINI_OBJECT               (_gst_mini_object_type)

#define GST_IS_MINI_OBJECT_TYPE(obj,type)  ((obj) && GST_MINI_OBJECT_TYPE(obj) == (type))
#define GST_MINI_OBJECT_CAST(obj)          ((GstMiniObject*)(obj))
#define GST_MINI_OBJECT_CONST_CAST(obj)    ((const GstMiniObject*)(obj))
#define GST_MINI_OBJECT(obj)               (GST_MINI_OBJECT_CAST(obj))

typedef struct _GstMiniObject GstMiniObject;

GST_API
GType           gst_mini_object_get_type   (void);

/**
 * GstMiniObjectCopyFunction:
 * @obj: MiniObject to copy
 *
 * Function prototype for methods to create copies of instances.
 *
 * Returns: reference to cloned instance.
 */
typedef GstMiniObject * (*GstMiniObjectCopyFunction) (const GstMiniObject *obj);
/**
 * GstMiniObjectDisposeFunction:
 * @obj: MiniObject to dispose
 *
 * Function prototype for when a miniobject has lost its last refcount.
 * Implementation of the mini object are allowed to revive the
 * passed object by doing a gst_mini_object_ref(). If the object is not
 * revived after the dispose function, the function should return %TRUE
 * and the memory associated with the object is freed.
 *
 * Returns: %TRUE if the object should be cleaned up.
 */
typedef gboolean (*GstMiniObjectDisposeFunction) (GstMiniObject *obj);
/**
 * GstMiniObjectFreeFunction:
 * @obj: MiniObject to free
 *
 * Virtual function prototype for methods to free resources used by
 * mini-objects.
 */
typedef void (*GstMiniObjectFreeFunction) (GstMiniObject *obj);

 /**
 * GstMiniObjectNotify:
 * @user_data: data that was provided when the notify was added
 * @obj: the mini object
 *
 * A #GstMiniObjectNotify function can be added to a mini object as a
 * callback that gets triggered when gst_mini_object_unref() drops the
 * last ref and @obj is about to be freed.
 */
typedef void (*GstMiniObjectNotify) (gpointer user_data, GstMiniObject * obj);

/**
 * GST_MINI_OBJECT_TYPE:
 * @obj: MiniObject to return type for.
 *
 * This macro returns the type of the mini-object.
 */
#define GST_MINI_OBJECT_TYPE(obj)  (GST_MINI_OBJECT_CAST(obj)->type)

/**
 * GST_MINI_OBJECT_FLAGS:
 * @obj: MiniObject to return flags for.
 *
 * This macro returns the entire set of flags for the mini-object.
 */
#define GST_MINI_OBJECT_FLAGS(obj)  (GST_MINI_OBJECT_CAST(obj)->flags)
/**
 * GST_MINI_OBJECT_FLAG_IS_SET:
 * @obj: MiniObject to check for flags.
 * @flag: Flag to check for
 *
 * This macro checks to see if the given flag is set.
 */
#define GST_MINI_OBJECT_FLAG_IS_SET(obj,flag)        !!(GST_MINI_OBJECT_FLAGS (obj) & (flag))
/**
 * GST_MINI_OBJECT_FLAG_SET:
 * @obj: MiniObject to set flag in.
 * @flag: Flag to set, can by any number of bits in guint32.
 *
 * This macro sets the given bits.
 */
#define GST_MINI_OBJECT_FLAG_SET(obj,flag)           (GST_MINI_OBJECT_FLAGS (obj) |= (flag))
/**
 * GST_MINI_OBJECT_FLAG_UNSET:
 * @obj: MiniObject to unset flag in.
 * @flag: Flag to set, must be a single bit in guint32.
 *
 * This macro unsets the given bits.
 */
#define GST_MINI_OBJECT_FLAG_UNSET(obj,flag)         (GST_MINI_OBJECT_FLAGS (obj) &= ~(flag))

/**
 * GstMiniObjectFlags:
 * @GST_MINI_OBJECT_FLAG_LOCKABLE: the object can be locked and unlocked with
 * gst_mini_object_lock() and gst_mini_object_unlock().
 * @GST_MINI_OBJECT_FLAG_LOCK_READONLY: the object is permanently locked in
 * READONLY mode. Only read locks can be performed on the object.
 * @GST_MINI_OBJECT_FLAG_MAY_BE_LEAKED: the object is expected to stay alive
 * even after gst_deinit() has been called and so should be ignored by leak
 * detection tools. (Since: 1.10)
 * @GST_MINI_OBJECT_FLAG_LAST: first flag that can be used by subclasses.
 *
 * Flags for the mini object
 */
typedef enum
{
  GST_MINI_OBJECT_FLAG_LOCKABLE      = (1 << 0),
  GST_MINI_OBJECT_FLAG_LOCK_READONLY = (1 << 1),
  GST_MINI_OBJECT_FLAG_MAY_BE_LEAKED = (1 << 2),
  /* padding */
  GST_MINI_OBJECT_FLAG_LAST          = (1 << 4)
} GstMiniObjectFlags;

/**
 * GST_MINI_OBJECT_IS_LOCKABLE:
 * @obj: a #GstMiniObject
 *
 * Check if @obj is lockable. A lockable object can be locked and unlocked with
 * gst_mini_object_lock() and gst_mini_object_unlock().
 */
#define GST_MINI_OBJECT_IS_LOCKABLE(obj)  GST_MINI_OBJECT_FLAG_IS_SET(obj, GST_MINI_OBJECT_FLAG_LOCKABLE)

/**
 * GstLockFlags:
 * @GST_LOCK_FLAG_READ: lock for read access
 * @GST_LOCK_FLAG_WRITE: lock for write access
 * @GST_LOCK_FLAG_EXCLUSIVE: lock for exclusive access
 * @GST_LOCK_FLAG_LAST: first flag that can be used for custom purposes
 *
 * Flags used when locking miniobjects
 */
typedef enum {
  GST_LOCK_FLAG_READ      = (1 << 0),
  GST_LOCK_FLAG_WRITE     = (1 << 1),
  GST_LOCK_FLAG_EXCLUSIVE = (1 << 2),

  GST_LOCK_FLAG_LAST      = (1 << 8)
} GstLockFlags;

/**
 * GST_LOCK_FLAG_READWRITE: (value 3) (type GstLockFlags)
 *
 * GstLockFlags value alias for GST_LOCK_FLAG_READ | GST_LOCK_FLAG_WRITE
 */
#define GST_LOCK_FLAG_READWRITE  ((GstLockFlags) (GST_LOCK_FLAG_READ | GST_LOCK_FLAG_WRITE))

/**
 * GST_MINI_OBJECT_REFCOUNT:
 * @obj: a #GstMiniObject
 *
 * Get access to the reference count field of the mini-object.
 */
#define GST_MINI_OBJECT_REFCOUNT(obj)           ((GST_MINI_OBJECT_CAST(obj))->refcount)
/**
 * GST_MINI_OBJECT_REFCOUNT_VALUE:
 * @obj: a #GstMiniObject
 *
 * Get the reference count value of the mini-object.
 */
#define GST_MINI_OBJECT_REFCOUNT_VALUE(obj)     (g_atomic_int_get (&(GST_MINI_OBJECT_CAST(obj))->refcount))

/**
 * GstMiniObject: (ref-func gst_mini_object_ref) (unref-func gst_mini_object_unref) (set-value-func g_value_set_boxed) (get-value-func g_value_get_boxed)
 * @type: the GType of the object
 * @refcount: atomic refcount
 * @lockstate: atomic state of the locks
 * @flags: extra flags.
 * @copy: a copy function
 * @dispose: a dispose function
 * @free: the free function
 *
 * Base class for refcounted lightweight objects.
 */
struct _GstMiniObject {
  GType   type;

  /*< public >*/ /* with COW */
  gint    refcount;
  gint    lockstate;
  guint   flags;

  GstMiniObjectCopyFunction copy;
  GstMiniObjectDisposeFunction dispose;
  GstMiniObjectFreeFunction free;

  /* < private > */
  /* Used to keep track of parents, weak ref notifies and qdata */
  guint priv_uint;
  gpointer priv_pointer;
};

GST_API
void            gst_mini_object_init (GstMiniObject *mini_object,
                                      guint flags, GType type,
                                      GstMiniObjectCopyFunction copy_func,
                                      GstMiniObjectDisposeFunction dispose_func,
                                      GstMiniObjectFreeFunction free_func);


/* refcounting */

GST_API
GstMiniObject * gst_mini_object_ref		(GstMiniObject *mini_object);

GST_API
void            gst_mini_object_unref		(GstMiniObject *mini_object);

GST_API
void        gst_clear_mini_object (GstMiniObject **object_ptr);
#define     gst_clear_mini_object(object_ptr) g_clear_pointer ((object_ptr), gst_mini_object_unref)

GST_API
void            gst_mini_object_weak_ref        (GstMiniObject *object,
					         GstMiniObjectNotify notify,
					         gpointer data);
GST_API
void            gst_mini_object_weak_unref	(GstMiniObject *object,
					         GstMiniObjectNotify notify,
					         gpointer data);

/* locking */

GST_API
gboolean        gst_mini_object_lock            (GstMiniObject *object, GstLockFlags flags);

GST_API
void            gst_mini_object_unlock          (GstMiniObject *object, GstLockFlags flags);

GST_API
gboolean        gst_mini_object_is_writable     (const GstMiniObject *mini_object);

GST_API
GstMiniObject * gst_mini_object_make_writable	(GstMiniObject *mini_object) G_GNUC_WARN_UNUSED_RESULT;

/* copy */

GST_API
GstMiniObject * gst_mini_object_copy		(const GstMiniObject *mini_object) G_GNUC_MALLOC G_GNUC_WARN_UNUSED_RESULT;


GST_API
void            gst_mini_object_set_qdata       (GstMiniObject *object, GQuark quark,
                                                 gpointer data, GDestroyNotify destroy);
GST_API
gpointer        gst_mini_object_get_qdata       (GstMiniObject *object, GQuark quark);

GST_API
gpointer        gst_mini_object_steal_qdata     (GstMiniObject *object, GQuark quark);

GST_API
void            gst_mini_object_add_parent      (GstMiniObject *object, GstMiniObject *parent);
GST_API
void            gst_mini_object_remove_parent   (GstMiniObject *object, GstMiniObject *parent);

GST_API
gboolean        gst_mini_object_replace         (GstMiniObject **olddata, GstMiniObject *newdata);

GST_API
gboolean        gst_mini_object_take            (GstMiniObject **olddata, GstMiniObject *newdata);

GST_API
GstMiniObject * gst_mini_object_steal           (GstMiniObject **olddata) G_GNUC_WARN_UNUSED_RESULT;

/**
 * GST_DEFINE_MINI_OBJECT_TYPE:
 * @TypeName: name of the new type in CamelCase
 * @type_name: name of the new type
 *
 * Define a new mini-object type with the given name
 */
#define GST_DEFINE_MINI_OBJECT_TYPE(TypeName,type_name) \
   G_DEFINE_BOXED_TYPE(TypeName,type_name,              \
       (GBoxedCopyFunc) gst_mini_object_ref,            \
       (GBoxedFreeFunc) gst_mini_object_unref)

G_END_DECLS

#endif

