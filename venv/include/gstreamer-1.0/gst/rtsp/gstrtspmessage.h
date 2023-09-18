/* GStreamer
 * Copyright (C) <2005,2006> Wim Taymans <wim@fluendo.com>
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

#ifndef __GST_RTSP_MESSAGE_H__
#define __GST_RTSP_MESSAGE_H__

#include <gst/gst.h>
#include <gst/rtsp/gstrtspdefs.h>

G_BEGIN_DECLS

/**
 * GstRTSPMsgType:
 * @GST_RTSP_MESSAGE_INVALID: invalid message type
 * @GST_RTSP_MESSAGE_REQUEST: RTSP request message
 * @GST_RTSP_MESSAGE_RESPONSE: RTSP response message
 * @GST_RTSP_MESSAGE_HTTP_REQUEST: HTTP request message.
 * @GST_RTSP_MESSAGE_HTTP_RESPONSE: HTTP response message.
 * @GST_RTSP_MESSAGE_DATA: data message
 *
 * The type of a message.
 */
typedef enum
{
  GST_RTSP_MESSAGE_INVALID,
  GST_RTSP_MESSAGE_REQUEST,
  GST_RTSP_MESSAGE_RESPONSE,
  GST_RTSP_MESSAGE_HTTP_REQUEST,
  GST_RTSP_MESSAGE_HTTP_RESPONSE,
  GST_RTSP_MESSAGE_DATA
} GstRTSPMsgType;

typedef struct _GstRTSPMessage GstRTSPMessage;

/**
 * GstRTSPMessage:
 * @type: the message type
 *
 * An RTSP message containing request, response or data messages. Depending on
 * the @type, the appropriate structure may be accessed.
 */
struct _GstRTSPMessage
{
  GstRTSPMsgType    type;

  union {
    struct {
      GstRTSPMethod      method;
      gchar             *uri;
      GstRTSPVersion     version;
    } request;
    struct {
      GstRTSPStatusCode  code;
      gchar             *reason;
      GstRTSPVersion     version;
    } response;
    struct {
      guint8             channel;
    } data;
  } type_data;

  /*< private >*/
  GArray        *hdr_fields;

  guint8        *body;
  guint          body_size;

  GstBuffer     *body_buffer;
  gpointer _gst_reserved[GST_PADDING-1];
};

GST_RTSP_API
GType                   gst_rtsp_msg_get_type            (void);

#define GST_TYPE_RTSP_MESSAGE           (gst_rtsp_msg_get_type())
#define GST_RTSP_MESSAGE_CAST(object)   ((GstRTSPMessage *)(object))
#define GST_RTSP_MESSAGE(object)        (GST_RTSP_MESSAGE_CAST(object))

/* memory management */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_new             (GstRTSPMessage **msg);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_init            (GstRTSPMessage *msg);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_unset           (GstRTSPMessage *msg);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_free            (GstRTSPMessage *msg);
GST_RTSP_API
GstRTSPResult      gst_rtsp_message_copy            (const GstRTSPMessage *msg,
                                                     GstRTSPMessage **copy);

GST_RTSP_API
GstRTSPMsgType     gst_rtsp_message_get_type        (GstRTSPMessage *msg);

/* request */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_new_request     (GstRTSPMessage **msg,
                                                     GstRTSPMethod method,
                                                     const gchar *uri);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_init_request    (GstRTSPMessage *msg,
                                                     GstRTSPMethod method,
                                                     const gchar *uri);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_parse_request   (GstRTSPMessage *msg,
                                                     GstRTSPMethod *method,
                                                     const gchar **uri,
                                                     GstRTSPVersion *version);

/* response */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_new_response    (GstRTSPMessage **msg,
                                                     GstRTSPStatusCode code,
                                                     const gchar *reason,
                                                     const GstRTSPMessage *request);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_init_response   (GstRTSPMessage *msg,
                                                     GstRTSPStatusCode code,
                                                     const gchar *reason,
                                                     const GstRTSPMessage *request);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_parse_response  (GstRTSPMessage *msg,
                                                     GstRTSPStatusCode *code,
                                                     const gchar **reason,
                                                     GstRTSPVersion *version);

/* data */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_new_data        (GstRTSPMessage **msg,
                                                     guint8 channel);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_init_data       (GstRTSPMessage *msg,
                                                     guint8 channel);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_parse_data      (GstRTSPMessage *msg,
                                                     guint8 *channel);

/* headers */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_add_header      (GstRTSPMessage *msg,
                                                     GstRTSPHeaderField field,
                                                     const gchar *value);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_take_header     (GstRTSPMessage *msg,
                                                     GstRTSPHeaderField field,
                                                     gchar *value);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_remove_header   (GstRTSPMessage *msg,
                                                     GstRTSPHeaderField field,
                                                     gint indx);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_get_header      (const GstRTSPMessage *msg,
                                                     GstRTSPHeaderField field,
                                                     gchar **value,
                                                     gint indx);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_add_header_by_name    (GstRTSPMessage * msg,
                                                           const gchar    * header,
                                                           const gchar    * value);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_take_header_by_name   (GstRTSPMessage * msg,
                                                           const gchar    * header,
                                                           gchar          * value);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_remove_header_by_name (GstRTSPMessage * msg,
                                                           const gchar    * header,
                                                           gint             index);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_get_header_by_name    (GstRTSPMessage * msg,
                                                           const gchar    * header,
                                                           gchar         ** value,
                                                           gint             index);

/* header serialization */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_append_headers  (const GstRTSPMessage *msg,
                                                     GString *str);

/* handling the body */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_set_body        (GstRTSPMessage *msg,
                                                     const guint8 *data,
                                                     guint size);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_take_body       (GstRTSPMessage *msg,
                                                     guint8 *data,
                                                     guint size);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_get_body        (const GstRTSPMessage *msg,
                                                     guint8 **data,
                                                     guint *size);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_steal_body      (GstRTSPMessage *msg,
                                                     guint8 **data,
                                                     guint *size);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_set_body_buffer (GstRTSPMessage *msg,
                                                     GstBuffer * buffer);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_take_body_buffer(GstRTSPMessage *msg,
                                                     GstBuffer * buffer);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_get_body_buffer (const GstRTSPMessage *msg,
                                                     GstBuffer ** buffer);

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_steal_body_buffer(GstRTSPMessage *msg,
                                                      GstBuffer ** buffer);

GST_RTSP_API
gboolean           gst_rtsp_message_has_body_buffer(const GstRTSPMessage *msg);

typedef struct _GstRTSPAuthCredential GstRTSPAuthCredential;
typedef struct _GstRTSPAuthParam GstRTSPAuthParam;

/**
 * GstRTSPAuthCredential:
 * @scheme: a #GstRTSPAuthMethod
 * @params: A NULL-terminated array of #GstRTSPAuthParam
 * @authorization: The authorization for the basic schem
 *
 * RTSP Authentication credentials
 *
 * Since: 1.12
 */
struct _GstRTSPAuthCredential {
  GstRTSPAuthMethod scheme;

  /* For Basic/Digest WWW-Authenticate and Digest
   * Authorization */
  GstRTSPAuthParam **params; /* NULL terminated */

  /* For Basic Authorization */
  gchar *authorization;
};

/**
 * GstRTSPAuthParam:
 * @name: The name of the parameter
 * @value: The value of the parameter
 *
 * RTSP Authentication parameter
 *
 * Since: 1.12
 */
struct _GstRTSPAuthParam {
  gchar *name;
  gchar *value;
};

GST_RTSP_API
GstRTSPAuthParam *       gst_rtsp_auth_param_copy (GstRTSPAuthParam * param);
GST_RTSP_API
void                     gst_rtsp_auth_param_free (GstRTSPAuthParam * param);

GST_RTSP_API
GstRTSPAuthCredential ** gst_rtsp_message_parse_auth_credentials (GstRTSPMessage * msg, GstRTSPHeaderField field);

GST_RTSP_API
void                     gst_rtsp_auth_credentials_free (GstRTSPAuthCredential ** credentials);

#define GST_TYPE_RTSP_AUTH_CREDENTIAL gst_rtsp_auth_credential_get_type()

GST_RTSP_API
GType                    gst_rtsp_auth_credential_get_type (void);

#define GST_TYPE_RTSP_AUTH_PARAM gst_rtsp_auth_param_get_type()

GST_RTSP_API
GType                    gst_rtsp_auth_param_get_type (void);

/* debug */

GST_RTSP_API
GstRTSPResult      gst_rtsp_message_dump            (GstRTSPMessage *msg);

G_END_DECLS

#endif /* __GST_RTSP_MESSAGE_H__ */
