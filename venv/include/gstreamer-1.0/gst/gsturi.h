/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *                    2014 David Waring, British Broadcasting Corporation
 *                             <david.waring@rd.bbc.co.uk>
 *
 * gsturi.h: Header for uri to element mappings and URI manipulation.
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


#ifndef __GST_URI_H__
#define __GST_URI_H__

#include <glib.h>
#include <glib-object.h>

G_BEGIN_DECLS

#include <gst/gstconfig.h>

/**
 * gst_uri_error_quark: (attributes doc.skip=true)
 */
GST_API
GQuark gst_uri_error_quark (void);

/**
 * GST_URI_ERROR:
 *
 * Get access to the error quark of the uri subsystem.
 */
#define GST_URI_ERROR gst_uri_error_quark ()

/**
 * GstURIError:
 * @GST_URI_ERROR_UNSUPPORTED_PROTOCOL: The protocol is not supported
 * @GST_URI_ERROR_BAD_URI: There was a problem with the URI
 * @GST_URI_ERROR_BAD_STATE: Could not set or change the URI because the
 *     URI handler was in a state where that is not possible or not permitted
 * @GST_URI_ERROR_BAD_REFERENCE: There was a problem with the entity that
 *     the URI references
 *
 * Different URI-related errors that can occur.
 */
typedef enum
{
  GST_URI_ERROR_UNSUPPORTED_PROTOCOL,
  GST_URI_ERROR_BAD_URI,
  GST_URI_ERROR_BAD_STATE,
  GST_URI_ERROR_BAD_REFERENCE
} GstURIError;

/**
 * GstURIType:
 * @GST_URI_UNKNOWN: The URI direction is unknown
 * @GST_URI_SINK: The URI is a consumer.
 * @GST_URI_SRC: The URI is a producer.
 *
 * The different types of URI direction.
 */

typedef enum {
  GST_URI_UNKNOWN,
  GST_URI_SINK,
  GST_URI_SRC
} GstURIType;

/**
 * GST_URI_TYPE_IS_VALID:
 * @type: A #GstURIType
 *
 * Tests if the type direction is valid.
 */
#define GST_URI_TYPE_IS_VALID(type) ((type) == GST_URI_SRC || (type) == GST_URI_SINK)

/* uri handler functions */
#define GST_TYPE_URI_HANDLER               (gst_uri_handler_get_type ())
#define GST_URI_HANDLER(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_URI_HANDLER, GstURIHandler))
#define GST_IS_URI_HANDLER(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_URI_HANDLER))
#define GST_URI_HANDLER_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), GST_TYPE_URI_HANDLER, GstURIHandlerInterface))

/**
 * GstURIHandler:
 *
 * Opaque #GstURIHandler structure.
 */
typedef struct _GstURIHandler GstURIHandler;
typedef struct _GstURIHandlerInterface GstURIHandlerInterface;

#include <gst/gstelement.h>
#include "gstminiobject.h"

/**
 * GstURIHandlerInterface:
 * @parent: The parent interface type
 * @get_type: Method to tell whether the element handles source or sink URI.
 * @get_protocols: Method to return the list of protocols handled by the element.
 * @get_uri: Method to return the URI currently handled by the element.
 * @set_uri: Method to set a new URI.
 *
 * Any #GstElement using this interface should implement these methods.
 */
struct _GstURIHandlerInterface {
  GTypeInterface parent;

  /* vtable */
  /*< public >*/
  /* querying capabilities */
  GstURIType             (* get_type)           (GType type);
  const gchar * const *  (* get_protocols)      (GType type);

  /* using the interface */
  gchar *                (* get_uri)            (GstURIHandler * handler);
  gboolean               (* set_uri)            (GstURIHandler * handler,
                                                 const gchar   * uri,
                                                 GError       ** error);
};

/* general URI functions */

GST_API
gboolean        gst_uri_protocol_is_valid       (const gchar * protocol);

GST_API
gboolean        gst_uri_protocol_is_supported   (const GstURIType type,
                                                 const gchar *protocol);
GST_API
gboolean        gst_uri_is_valid                (const gchar * uri);

GST_API
gchar *         gst_uri_get_protocol            (const gchar * uri) G_GNUC_MALLOC;

GST_API
gboolean        gst_uri_has_protocol            (const gchar * uri,
                                                 const gchar * protocol);
GST_API
gchar *         gst_uri_get_location            (const gchar * uri) G_GNUC_MALLOC;

GST_DEPRECATED_FOR(gst_uri_new)
gchar *         gst_uri_construct               (const gchar * protocol,
                                                 const gchar * location) G_GNUC_MALLOC;
GST_API
gchar *         gst_filename_to_uri             (const gchar * filename,
                                                 GError     ** error) G_GNUC_MALLOC;
GST_API
GstElement *    gst_element_make_from_uri       (const GstURIType type,
                                                 const gchar *    uri,
                                                 const gchar *    elementname,
                                                 GError      **   error) G_GNUC_MALLOC;

/* accessing the interface */

GST_API
GType                 gst_uri_handler_get_type (void);

GST_API
GstURIType            gst_uri_handler_get_uri_type  (GstURIHandler * handler);

GST_API
const gchar * const * gst_uri_handler_get_protocols (GstURIHandler * handler);

GST_API
gchar *               gst_uri_handler_get_uri       (GstURIHandler * handler) G_GNUC_MALLOC;

GST_API
gboolean              gst_uri_handler_set_uri       (GstURIHandler * handler,
                                                     const gchar   * uri,
                                                     GError       ** error);

/*
 * GstUri Type macros.
 */
#define GST_TYPE_URI        (gst_uri_get_type ())
#define GST_IS_URI(obj)	    (GST_IS_MINI_OBJECT_TYPE (obj, GST_TYPE_URI))
#define GST_URI_CAST(obj)   ((GstUri *)(obj))
#define GST_URI_CONST_CAST(obj) ((const GstUri *)(obj))
#define GST_URI(obj)        (GST_URI_CAST(obj))

/**
 * GstUri:
 *
 * This is a private structure that holds the various parts of a parsed URI.
 */
struct _GstUri;
typedef struct _GstUri GstUri;

/**
 * GST_URI_NO_PORT:
 *
 * Value for #GstUri<!-- -->.port to indicate no port number.
 */
#define GST_URI_NO_PORT 0

/* used by GST_TYPE_URI */

GST_API
GType gst_uri_get_type (void);

/*
 * Method definitions.
 */

GST_API
GstUri * gst_uri_new                   (const gchar * scheme,
                                        const gchar * userinfo,
                                        const gchar * host,
                                        guint port,
                                        const gchar * path,
                                        const gchar * query,
                                        const gchar * fragment) G_GNUC_MALLOC;
GST_API
GstUri * gst_uri_new_with_base         (GstUri * base,
                                        const gchar * scheme,
                                        const gchar * userinfo,
                                        const gchar * host,
                                        guint port,
                                        const gchar * path,
                                        const gchar * query,
                                        const gchar * fragment) G_GNUC_MALLOC;
GST_API
GstUri * gst_uri_from_string           (const gchar * uri) G_GNUC_MALLOC;

GST_API
GstUri * gst_uri_from_string_escaped   (const gchar * uri) G_GNUC_MALLOC;

GST_API
GstUri * gst_uri_from_string_with_base (GstUri * base,
                                        const gchar * uri) G_GNUC_MALLOC;
GST_API
gboolean gst_uri_equal                 (const GstUri * first,
                                        const GstUri * second);
GST_API
GstUri * gst_uri_join                  (GstUri * base_uri,
                                        GstUri * ref_uri) G_GNUC_WARN_UNUSED_RESULT;
GST_API
gchar *  gst_uri_join_strings          (const gchar * base_uri,
                                        const gchar * ref_uri) G_GNUC_MALLOC;
GST_API
gboolean gst_uri_is_writable           (const GstUri * uri);

GST_API
GstUri * gst_uri_make_writable         (GstUri * uri) G_GNUC_WARN_UNUSED_RESULT;

GST_API
gchar *  gst_uri_to_string             (const GstUri * uri) G_GNUC_MALLOC;

GST_API
gboolean gst_uri_is_normalized         (const GstUri * uri);

GST_API
gboolean gst_uri_normalize             (GstUri * uri);

GST_API
const gchar * gst_uri_get_scheme       (const GstUri * uri);

GST_API
gboolean gst_uri_set_scheme            (GstUri * uri, const gchar * scheme);

GST_API
const gchar * gst_uri_get_userinfo     (const GstUri * uri);

GST_API
gboolean gst_uri_set_userinfo          (GstUri * uri, const gchar * userinfo);

GST_API
const gchar * gst_uri_get_host         (const GstUri * uri);

GST_API
gboolean gst_uri_set_host              (GstUri * uri, const gchar * host);

GST_API
guint gst_uri_get_port                 (const GstUri * uri);

GST_API
gboolean gst_uri_set_port              (GstUri * uri, guint port);

GST_API
gchar * gst_uri_get_path               (const GstUri * uri);

GST_API
gboolean gst_uri_set_path              (GstUri * uri, const gchar * path);

GST_API
gchar * gst_uri_get_path_string        (const GstUri * uri);

GST_API
gboolean gst_uri_set_path_string       (GstUri * uri, const gchar * path);

GST_API
GList * gst_uri_get_path_segments      (const GstUri * uri);

GST_API
gboolean gst_uri_set_path_segments     (GstUri * uri, GList * path_segments);

GST_API
gboolean gst_uri_append_path           (GstUri * uri,
                                        const gchar * relative_path);
GST_API
gboolean gst_uri_append_path_segment   (GstUri * uri,
                                        const gchar * path_segment);
GST_API
gchar * gst_uri_get_query_string       (const GstUri * uri);

GST_API
gboolean gst_uri_set_query_string      (GstUri * uri, const gchar * query);

GST_API
GHashTable * gst_uri_get_query_table   (const GstUri * uri);

GST_API
gboolean gst_uri_set_query_table       (GstUri * uri,
                                        GHashTable * query_table);
GST_API
gboolean gst_uri_set_query_value       (GstUri * uri, const gchar * query_key,
                                        const gchar * query_value);
GST_API
gboolean gst_uri_remove_query_key      (GstUri * uri, const gchar * query_key);

GST_API
gboolean gst_uri_query_has_key         (const GstUri * uri,
                                        const gchar * query_key);

GST_API
const gchar * gst_uri_get_query_value  (const GstUri * uri,
                                        const gchar * query_key);

GST_API
GList * gst_uri_get_query_keys         (const GstUri * uri);

GST_API
const gchar * gst_uri_get_fragment     (const GstUri * uri);

GST_API
gboolean gst_uri_set_fragment          (GstUri * uri, const gchar * fragment);

GST_API
GHashTable * gst_uri_get_media_fragment_table  (const GstUri * uri);

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
static inline GstUri *
gst_uri_copy (const GstUri * uri)
{
  return GST_URI_CAST (gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (uri)));
}

static inline GstUri *
gst_uri_ref (GstUri * uri)
{
  return GST_URI_CAST (gst_mini_object_ref (GST_MINI_OBJECT_CAST (uri)));
}

static inline void
gst_uri_unref (GstUri * uri)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (uri));
}

static inline void
gst_clear_uri (GstUri ** uri_ptr)
{
  gst_clear_mini_object ((GstMiniObject **) uri_ptr);
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstUri * gst_uri_copy  (const GstUri * uri);

GST_API
GstUri * gst_uri_ref   (GstUri * uri);

GST_API
void     gst_uri_unref (GstUri * uri);

GST_API
void     gst_clear_uri (GstUri ** uri_ptr);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstUri, gst_uri_unref)

G_END_DECLS

#endif /* __GST_URI_H__ */
