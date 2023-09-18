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

#ifndef __GST_RTSP_URL_H__
#define __GST_RTSP_URL_H__

#include <glib.h>
#include <glib-object.h>

#include <gst/rtsp/gstrtspdefs.h>
#include <gst/rtsp/gstrtsptransport.h>

G_BEGIN_DECLS

/**
 * GST_RTSP_DEFAULT_PORT:
 *
 * The default RTSP port to connect to.
 */
#define GST_RTSP_DEFAULT_PORT       554

#define GST_TYPE_RTSP_URL  (gst_rtsp_url_get_type())

typedef struct _GstRTSPUrl GstRTSPUrl;

/**
 * GstRTSPUrl:
 * @transports: the transports allowed
 * @family: the family
 * @user: the user
 * @passwd: the password
 * @host: the host
 * @port: the port
 * @abspath: the absolute path
 * @query: additional query parameters
 *
 * This structure contains the result of a parsed RTSP URL
 */
struct _GstRTSPUrl {
  GstRTSPLowerTrans  transports;
  GstRTSPFamily      family;
  gchar             *user;
  gchar             *passwd;
  gchar             *host;
  guint16            port;
  gchar             *abspath;
  gchar             *query;
};

GST_RTSP_API
GType gst_rtsp_url_get_type (void);

GST_RTSP_API
GstRTSPResult      gst_rtsp_url_parse           (const gchar *urlstr, GstRTSPUrl **url);

GST_RTSP_API
GstRTSPUrl*        gst_rtsp_url_copy            (const GstRTSPUrl *url);

GST_RTSP_API
void               gst_rtsp_url_free            (GstRTSPUrl *url);

GST_RTSP_API
gchar*             gst_rtsp_url_get_request_uri (const GstRTSPUrl *url);

GST_RTSP_API
gchar *            gst_rtsp_url_get_request_uri_with_control (const GstRTSPUrl * url,
                                                              const gchar * control_path);

GST_RTSP_API
gchar**            gst_rtsp_url_decode_path_components
                                                (const GstRTSPUrl *url);
GST_RTSP_API
GstRTSPResult      gst_rtsp_url_set_port        (GstRTSPUrl *url, guint16 port);

GST_RTSP_API
GstRTSPResult      gst_rtsp_url_get_port        (const GstRTSPUrl *url, guint16 *port);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstRTSPUrl, gst_rtsp_url_free)

G_END_DECLS

#endif /* __GST_RTSP_URL_H__ */
