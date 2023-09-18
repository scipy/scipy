/* GStreamer
 * Copyright (C) 2013 Collabora Ltd.
 *   Author: Sebastian Dröge <sebastian.droege@collabora.co.uk>
 * Copyright (C) 2013 Sebastian Dröge <slomo@circular-chaos.org>
 *
 * gstcontext.h: Header for GstContext subsystem
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

#ifndef __GST_CONTEXT_H__
#define __GST_CONTEXT_H__

#include <glib.h>

G_BEGIN_DECLS

typedef struct _GstContext GstContext;

#include <gst/gstminiobject.h>
#include <gst/gststructure.h>

GST_API GType _gst_context_type;

#define GST_TYPE_CONTEXT                         (_gst_context_type)
#define GST_IS_CONTEXT(obj)                      (GST_IS_MINI_OBJECT_TYPE (obj, GST_TYPE_CONTEXT))
#define GST_CONTEXT_CAST(obj)                    ((GstContext*)(obj))
#define GST_CONTEXT(obj)                         (GST_CONTEXT_CAST(obj))



GST_API
GType           gst_context_get_type            (void);

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* refcounting */
static inline GstContext *
gst_context_ref (GstContext * context)
{
  return (GstContext *) gst_mini_object_ref (GST_MINI_OBJECT_CAST (context));
}

static inline void
gst_context_unref (GstContext * context)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (context));
}

/* copy context */
static inline GstContext *
gst_context_copy (const GstContext * context)
{
  return GST_CONTEXT_CAST (gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (context)));
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstContext * gst_context_ref    (GstContext * context);

GST_API
void         gst_context_unref  (GstContext * context);

GST_API
GstContext * gst_context_copy   (const GstContext * context);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

/**
 * gst_context_is_writable:
 * @context: a #GstContext
 *
 * Tests if you can safely write into a context's structure or validly
 * modify the seqnum and timestamp fields.
 */
#define         gst_context_is_writable(context)     gst_mini_object_is_writable (GST_MINI_OBJECT_CAST (context))
/**
 * gst_context_make_writable:
 * @context: (transfer full): the context to make writable
 *
 * Checks if a context is writable. If not, a writable copy is made and
 * returned.
 *
 * Returns: (transfer full): a context (possibly a duplicate) that is writable.
 *
 * MT safe
 */
#define         gst_context_make_writable(context)  GST_CONTEXT_CAST (gst_mini_object_make_writable (GST_MINI_OBJECT_CAST (context)))

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
static inline gboolean
gst_context_replace (GstContext **old_context, GstContext *new_context)
{
  return gst_mini_object_replace ((GstMiniObject **) old_context, (GstMiniObject *) new_context);
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
gboolean              gst_context_replace                  (GstContext ** old_context,
                                                            GstContext * new_context);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

GST_API
GstContext *          gst_context_new                      (const gchar * context_type,
                                                            gboolean persistent) G_GNUC_MALLOC;
GST_API
const gchar *         gst_context_get_context_type         (const GstContext * context);

GST_API
gboolean              gst_context_has_context_type         (const GstContext * context, const gchar * context_type);

GST_API
const GstStructure *  gst_context_get_structure            (const GstContext * context);

GST_API
GstStructure *        gst_context_writable_structure       (GstContext * context);

GST_API
gboolean              gst_context_is_persistent            (const GstContext * context);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstContext, gst_context_unref)

G_END_DECLS

#endif /* __GST_CONTEXT_H__ */
