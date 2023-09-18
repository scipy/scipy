/* GStreamer
 * Copyright (C) <2005,2009> Wim Taymans <wim.taymans@gmail.com>
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
/*
 * Unless otherwise indicated, Source Code is licensed under MIT license.
 * See further explanation attached in License Statement (distributed in the file
 * LICENSE).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __GST_RTSP_CONNECTION_H__
#define __GST_RTSP_CONNECTION_H__

#include <glib.h>

#include <gst/gstconfig.h>
#include <gst/rtsp/gstrtspdefs.h>
#include <gst/rtsp/gstrtspurl.h>
#include <gst/rtsp/gstrtspmessage.h>
#include <gio/gio.h>

G_BEGIN_DECLS

/**
 * GstRTSPConnection:
 *
 * Opaque RTSP connection object.
 */
typedef struct _GstRTSPConnection GstRTSPConnection;

/* opening/closing a connection */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_create         (const GstRTSPUrl *url, GstRTSPConnection **conn);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_create_from_socket (GSocket * socket,
                                                       const gchar * ip,
                                                       guint16 port,
                                                       const gchar * initial_buffer,
                                                       GstRTSPConnection ** conn);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_accept                 (GSocket * socket, GstRTSPConnection ** conn, GCancellable * cancellable);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_connect_usec           (GstRTSPConnection * conn, gint64 timeout);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_connect_with_response_usec (GstRTSPConnection * conn, gint64 timeout, GstRTSPMessage * response);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_close                  (GstRTSPConnection *conn);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_free                   (GstRTSPConnection *conn);

/* TLS connections */

GST_RTSP_API
GTlsConnection *     gst_rtsp_connection_get_tls                  (GstRTSPConnection * conn, GError ** error);

GST_RTSP_API
gboolean             gst_rtsp_connection_set_tls_validation_flags (GstRTSPConnection * conn, GTlsCertificateFlags flags);

GST_RTSP_API
GTlsCertificateFlags gst_rtsp_connection_get_tls_validation_flags (GstRTSPConnection * conn);

GST_RTSP_API
void                 gst_rtsp_connection_set_tls_database (GstRTSPConnection * conn, GTlsDatabase * database);

GST_RTSP_API
GTlsDatabase *       gst_rtsp_connection_get_tls_database (GstRTSPConnection * conn);

GST_RTSP_API
void                 gst_rtsp_connection_set_tls_interaction (GstRTSPConnection * conn, GTlsInteraction * interaction);

GST_RTSP_API
GTlsInteraction *    gst_rtsp_connection_get_tls_interaction (GstRTSPConnection * conn);

typedef gboolean (*GstRTSPConnectionAcceptCertificateFunc) (GTlsConnection *conn,
                                                            GTlsCertificate *peer_cert,
                                                            GTlsCertificateFlags errors,
                                                            gpointer user_data);
GST_RTSP_API
void                 gst_rtsp_connection_set_accept_certificate_func (GstRTSPConnection * conn,
                                                                      GstRTSPConnectionAcceptCertificateFunc func,
                                                                      gpointer user_data,
                                                                      GDestroyNotify destroy_notify);

/* sending/receiving raw bytes */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_read_usec      (GstRTSPConnection * conn, guint8 * data,
                                                       guint size, gint64 timeout);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_write_usec     (GstRTSPConnection * conn, const guint8 * data,
                                                       guint size, gint64 timeout);

/* sending/receiving messages */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_send_usec      (GstRTSPConnection *conn, GstRTSPMessage *message,
                                                       gint64 timeout);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_send_messages_usec (GstRTSPConnection *conn, GstRTSPMessage *messages, guint n_messages,
                                                       gint64 timeout);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_receive_usec    (GstRTSPConnection *conn, GstRTSPMessage *message,
                                                       gint64 timeout);

/* status management */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_poll_usec      (GstRTSPConnection *conn, GstRTSPEvent events,
                                                       GstRTSPEvent *revents, gint64 timeout);

/* reset the timeout */

GST_RTSP_API
gint64             gst_rtsp_connection_next_timeout_usec (GstRTSPConnection *conn);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_reset_timeout  (GstRTSPConnection *conn);

/* flushing state */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_flush          (GstRTSPConnection *conn, gboolean flush);

/* HTTP proxy support */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_set_proxy      (GstRTSPConnection *conn,
                                                       const gchar *host, guint port);

/* configure authentication data */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_set_auth       (GstRTSPConnection *conn, GstRTSPAuthMethod method,
                                                       const gchar *user, const gchar *pass);

GST_RTSP_API
void               gst_rtsp_connection_set_auth_param    (GstRTSPConnection *conn,
                                                          const gchar * param,
                                                          const gchar *value);

GST_RTSP_API
void               gst_rtsp_connection_clear_auth_params (GstRTSPConnection *conn);

/* configure DSCP */

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_set_qos_dscp   (GstRTSPConnection *conn,
                                                       guint qos_dscp);

/* Content-Length limit */
GST_RTSP_API
void               gst_rtsp_connection_set_content_length_limit (GstRTSPConnection *conn,
                                                                 guint limit);

/* accessors */

GST_RTSP_API
GstRTSPUrl *       gst_rtsp_connection_get_url        (const GstRTSPConnection *conn);

GST_RTSP_API
const gchar *      gst_rtsp_connection_get_ip         (const GstRTSPConnection *conn);

GST_RTSP_API
void               gst_rtsp_connection_set_ip         (GstRTSPConnection *conn, const gchar *ip);

GST_RTSP_API
GSocket *          gst_rtsp_connection_get_read_socket  (const GstRTSPConnection *conn);

GST_RTSP_API
GSocket *          gst_rtsp_connection_get_write_socket (const GstRTSPConnection *conn);

GST_RTSP_API
void               gst_rtsp_connection_set_http_mode  (GstRTSPConnection *conn,
                                                       gboolean enable);

/* tunneling */

GST_RTSP_API
void               gst_rtsp_connection_set_tunneled   (GstRTSPConnection *conn, gboolean tunneled);

GST_RTSP_API
gboolean           gst_rtsp_connection_is_tunneled    (const GstRTSPConnection *conn);

GST_RTSP_API
const gchar *      gst_rtsp_connection_get_tunnelid   (const GstRTSPConnection *conn);

GST_RTSP_API
GstRTSPResult      gst_rtsp_connection_do_tunnel      (GstRTSPConnection *conn, GstRTSPConnection *conn2);

GST_RTSP_API
void               gst_rtsp_connection_set_remember_session_id (GstRTSPConnection *conn, gboolean remember);

GST_RTSP_API
gboolean           gst_rtsp_connection_get_remember_session_id (GstRTSPConnection *conn);

GST_RTSP_API
void               gst_rtsp_connection_set_ignore_x_server_reply (GstRTSPConnection *conn, gboolean ignore);

GST_RTSP_API
gboolean           gst_rtsp_connection_get_ignore_x_server_reply (const GstRTSPConnection *conn);

/* async IO */

/**
 * GstRTSPWatch:
 *
 * Opaque RTSP watch object that can be used for asynchronous RTSP
 * operations.
 */
typedef struct _GstRTSPWatch GstRTSPWatch;

/**
 * GstRTSPWatchFuncs:
 * @message_received: callback when a message was received
 * @message_sent: callback when a message was sent
 * @closed: callback when the connection is closed
 * @error: callback when an error occurred
 * @tunnel_start: a client started a tunneled connection. The tunnelid of the
 *   connection must be saved.
 * @tunnel_complete: a client finished a tunneled connection. In this callback
 *   you usually pair the tunnelid of this connection with the saved one using
 *   gst_rtsp_connection_do_tunnel().
 * @error_full: callback when an error occurred with more information than
 *   the @error callback.
 * @tunnel_lost: callback when the post connection of a tunnel is closed.
 * @tunnel_http_response: callback when an HTTP response to the GET request
 *   is about to be sent for a tunneled connection. The response can be
 *   modified in the callback. Since: 1.4.
 *
 * Callback functions from a #GstRTSPWatch.
 */
typedef struct {
  GstRTSPResult     (*message_received) (GstRTSPWatch *watch, GstRTSPMessage *message,
                                         gpointer user_data);
  GstRTSPResult     (*message_sent)     (GstRTSPWatch *watch, guint id,
                                         gpointer user_data);
  GstRTSPResult     (*closed)           (GstRTSPWatch *watch, gpointer user_data);
  GstRTSPResult     (*error)            (GstRTSPWatch *watch, GstRTSPResult result,
                                         gpointer user_data);
  GstRTSPStatusCode (*tunnel_start)     (GstRTSPWatch *watch, gpointer user_data);
  GstRTSPResult     (*tunnel_complete)  (GstRTSPWatch *watch, gpointer user_data);
  GstRTSPResult     (*error_full)       (GstRTSPWatch *watch, GstRTSPResult result,
                                         GstRTSPMessage *message, guint id,
                                         gpointer user_data);
  GstRTSPResult     (*tunnel_lost)      (GstRTSPWatch *watch, gpointer user_data);
  GstRTSPResult     (*tunnel_http_response) (GstRTSPWatch *watch,
                                             GstRTSPMessage *request,
                                             GstRTSPMessage *response,
                                             gpointer user_data);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING-1];
} GstRTSPWatchFuncs;

GST_RTSP_API
GstRTSPWatch *     gst_rtsp_watch_new                (GstRTSPConnection *conn,
                                                      GstRTSPWatchFuncs *funcs,
                                                      gpointer user_data,
                                                      GDestroyNotify notify);

GST_RTSP_API
void               gst_rtsp_watch_reset              (GstRTSPWatch *watch);

GST_RTSP_API
void               gst_rtsp_watch_unref              (GstRTSPWatch *watch);

GST_RTSP_API
guint              gst_rtsp_watch_attach             (GstRTSPWatch *watch,
                                                      GMainContext *context);

GST_RTSP_API
void               gst_rtsp_watch_set_send_backlog  (GstRTSPWatch *watch,
                                                     gsize bytes, guint messages);

GST_RTSP_API
void               gst_rtsp_watch_get_send_backlog  (GstRTSPWatch *watch,
                                                     gsize *bytes, guint *messages);

GST_RTSP_API
GstRTSPResult      gst_rtsp_watch_write_data         (GstRTSPWatch *watch,
                                                      const guint8 *data,
                                                      guint size, guint *id);

GST_RTSP_API
GstRTSPResult      gst_rtsp_watch_send_message       (GstRTSPWatch *watch,
                                                      GstRTSPMessage *message,
                                                      guint *id);

GST_RTSP_API
GstRTSPResult      gst_rtsp_watch_send_messages      (GstRTSPWatch *watch,
                                                      GstRTSPMessage *messages,
                                                      guint n_messages,
                                                      guint *id);

GST_RTSP_API
GstRTSPResult      gst_rtsp_watch_wait_backlog_usec  (GstRTSPWatch * watch,
                                                      gint64 timeout);

GST_RTSP_API
void               gst_rtsp_watch_set_flushing       (GstRTSPWatch * watch,
                                                      gboolean flushing);

#ifndef GST_DISABLE_DEPRECATED

/* Deprecated */
G_GNUC_BEGIN_IGNORE_DEPRECATIONS

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_connect_usec)
GstRTSPResult      gst_rtsp_connection_connect                (GstRTSPConnection * conn, GTimeVal * timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_connect_with_response_usec)
GstRTSPResult      gst_rtsp_connection_connect_with_response  (GstRTSPConnection * conn, GTimeVal * timeout, GstRTSPMessage * response);


GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_read_usec)
GstRTSPResult      gst_rtsp_connection_read           (GstRTSPConnection * conn, guint8 * data,
                                                       guint size, GTimeVal * timeout);


GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_write_usec)
GstRTSPResult      gst_rtsp_connection_write          (GstRTSPConnection * conn, const guint8 * data,
                                                       guint size, GTimeVal * timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_send_usec)
GstRTSPResult      gst_rtsp_connection_send           (GstRTSPConnection *conn, GstRTSPMessage *message,
                                                       GTimeVal *timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_send_messages_usec)
GstRTSPResult      gst_rtsp_connection_send_messages  (GstRTSPConnection *conn, GstRTSPMessage *messages, guint n_messages,
                                                       GTimeVal *timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_receive_usec)
GstRTSPResult      gst_rtsp_connection_receive        (GstRTSPConnection *conn, GstRTSPMessage *message,
                                                       GTimeVal *timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_poll_usec)
GstRTSPResult      gst_rtsp_connection_poll           (GstRTSPConnection *conn, GstRTSPEvent events,
                                                       GstRTSPEvent *revents, GTimeVal *timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_connection_next_timeout_usec)
GstRTSPResult      gst_rtsp_connection_next_timeout   (GstRTSPConnection *conn, GTimeVal *timeout);

GST_RTSP_DEPRECATED_FOR (gst_rtsp_watch_wait_backlog_usec)
GstRTSPResult      gst_rtsp_watch_wait_backlog       (GstRTSPWatch * watch,
                                                      GTimeVal *timeout);

G_GNUC_END_IGNORE_DEPRECATIONS

#endif /* GST_DISABLE_DEPRECATED */

G_END_DECLS

#endif /* __GST_RTSP_CONNECTION_H__ */
