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

#ifndef __GST_RTSP_DEFS_H__
#define __GST_RTSP_DEFS_H__

#include <gst/gst.h>
#include <gst/rtsp/rtsp-prelude.h>

G_BEGIN_DECLS

/**
 * GST_RTSP_CHECK:
 * @stmt: a statement
 * @label: a label
 *
 * Macro that checks the return value of @stmt and jumps to @label when it does
 * not equal #GST_RTSP_OK.
 */
#define GST_RTSP_CHECK(stmt, label)  \
G_STMT_START { \
  if (G_UNLIKELY ((res = (stmt)) != GST_RTSP_OK)) \
    goto label; \
} G_STMT_END

/**
 * GstRTSPResult:
 * @GST_RTSP_OK: no error
 * @GST_RTSP_ERROR: some unspecified error occurred
 * @GST_RTSP_EINVAL: invalid arguments were provided to a function
 * @GST_RTSP_EINTR: an operation was canceled
 * @GST_RTSP_ENOMEM: no memory was available for the operation
 * @GST_RTSP_ERESOLV: a host resolve error occurred
 * @GST_RTSP_ENOTIMPL: function not implemented
 * @GST_RTSP_ESYS: a system error occurred, errno contains more details
 * @GST_RTSP_EPARSE: a parsing error occurred
 * @GST_RTSP_EWSASTART: windows networking could not start
 * @GST_RTSP_EWSAVERSION: windows networking stack has wrong version
 * @GST_RTSP_EEOF: end-of-file was reached
 * @GST_RTSP_ENET: a network problem occurred, h_errno contains more details
 * @GST_RTSP_ENOTIP: the host is not an IP host
 * @GST_RTSP_ETIMEOUT: a timeout occurred
 * @GST_RTSP_ETGET: the tunnel GET request has been performed
 * @GST_RTSP_ETPOST: the tunnel POST request has been performed
 * @GST_RTSP_ELAST: last error
 *
 * Result codes from the RTSP functions.
 */
typedef enum {
  GST_RTSP_OK          =  0,
  /* errors */
  GST_RTSP_ERROR       = -1,
  GST_RTSP_EINVAL      = -2,
  GST_RTSP_EINTR       = -3,
  GST_RTSP_ENOMEM      = -4,
  GST_RTSP_ERESOLV     = -5,
  GST_RTSP_ENOTIMPL    = -6,
  GST_RTSP_ESYS        = -7,
  GST_RTSP_EPARSE      = -8,
  GST_RTSP_EWSASTART   = -9,
  GST_RTSP_EWSAVERSION = -10,
  GST_RTSP_EEOF        = -11,
  GST_RTSP_ENET        = -12,
  GST_RTSP_ENOTIP      = -13,
  GST_RTSP_ETIMEOUT    = -14,
  GST_RTSP_ETGET       = -15,
  GST_RTSP_ETPOST      = -16,

  GST_RTSP_ELAST       = -17
} GstRTSPResult;

/**
 * GstRTSPEvent:
 * @GST_RTSP_EV_READ: connection is readable
 * @GST_RTSP_EV_WRITE: connection is writable
 *
 * The possible events for the connection.
 */
typedef enum {
  GST_RTSP_EV_READ  = (1 << 0),
  GST_RTSP_EV_WRITE = (1 << 1)
} GstRTSPEvent;

/**
 * GstRTSPFamily:
 * @GST_RTSP_FAM_NONE: unknown network family
 * @GST_RTSP_FAM_INET: internet
 * @GST_RTSP_FAM_INET6: internet V6
 *
 * The possible network families.
 */
typedef enum {
  GST_RTSP_FAM_NONE,
  GST_RTSP_FAM_INET,
  GST_RTSP_FAM_INET6
} GstRTSPFamily;

/**
 * GstRTSPState:
 * @GST_RTSP_STATE_INVALID: invalid state
 * @GST_RTSP_STATE_INIT: initializing
 * @GST_RTSP_STATE_READY: ready for operation
 * @GST_RTSP_STATE_SEEKING: seeking in progress
 * @GST_RTSP_STATE_PLAYING: playing
 * @GST_RTSP_STATE_RECORDING: recording
 *
 * The different RTSP states.
 */
typedef enum {
  GST_RTSP_STATE_INVALID,
  GST_RTSP_STATE_INIT,
  GST_RTSP_STATE_READY,
  GST_RTSP_STATE_SEEKING,
  GST_RTSP_STATE_PLAYING,
  GST_RTSP_STATE_RECORDING
} GstRTSPState;

/**
 * GstRTSPVersion:
 * @GST_RTSP_VERSION_INVALID: unknown/invalid version
 * @GST_RTSP_VERSION_1_0: version 1.0
 * @GST_RTSP_VERSION_1_1: version 1.1.
 * @GST_RTSP_VERSION_2_0: version 2.0.
 *
 * The supported RTSP versions.
 */
typedef enum {
  GST_RTSP_VERSION_INVALID = 0x00,
  GST_RTSP_VERSION_1_0     = 0x10,
  GST_RTSP_VERSION_1_1     = 0x11,
  GST_RTSP_VERSION_2_0     = 0x20
} GstRTSPVersion;

/**
 * GstRTSPMethod:
 * @GST_RTSP_INVALID: invalid method
 * @GST_RTSP_DESCRIBE: the DESCRIBE method
 * @GST_RTSP_ANNOUNCE: the ANNOUNCE method
 * @GST_RTSP_GET_PARAMETER: the GET_PARAMETER method
 * @GST_RTSP_OPTIONS: the OPTIONS method
 * @GST_RTSP_PAUSE: the PAUSE method
 * @GST_RTSP_PLAY: the PLAY method
 * @GST_RTSP_RECORD: the RECORD method
 * @GST_RTSP_REDIRECT: the REDIRECT method
 * @GST_RTSP_SETUP: the SETUP method
 * @GST_RTSP_SET_PARAMETER: the SET_PARAMETER method
 * @GST_RTSP_TEARDOWN: the TEARDOWN method
 * @GST_RTSP_GET: the GET method (HTTP).
 * @GST_RTSP_POST: the POST method (HTTP).
 *
 * The different supported RTSP methods.
 */
typedef enum {
  GST_RTSP_INVALID          = 0,
  GST_RTSP_DESCRIBE         = (1 <<  0),
  GST_RTSP_ANNOUNCE         = (1 <<  1),
  GST_RTSP_GET_PARAMETER    = (1 <<  2),
  GST_RTSP_OPTIONS          = (1 <<  3),
  GST_RTSP_PAUSE            = (1 <<  4),
  GST_RTSP_PLAY             = (1 <<  5),
  GST_RTSP_RECORD           = (1 <<  6),
  GST_RTSP_REDIRECT         = (1 <<  7),
  GST_RTSP_SETUP            = (1 <<  8),
  GST_RTSP_SET_PARAMETER    = (1 <<  9),
  GST_RTSP_TEARDOWN         = (1 << 10),
  GST_RTSP_GET              = (1 << 11),
  GST_RTSP_POST             = (1 << 12)
} GstRTSPMethod;

/**
 * GstRTSPAuthMethod:
 * @GST_RTSP_AUTH_NONE: no authentication
 * @GST_RTSP_AUTH_BASIC: basic authentication
 * @GST_RTSP_AUTH_DIGEST: digest authentication
 *
 * Authentication methods, ordered by strength
 */
typedef enum {
  GST_RTSP_AUTH_NONE    = 0x00,
  GST_RTSP_AUTH_BASIC   = 0x01,
  GST_RTSP_AUTH_DIGEST  = 0x02
} GstRTSPAuthMethod;

/**
 * GST_RTSP_AUTH_MAX:
 *
 * Strongest available authentication method
 */
#define GST_RTSP_AUTH_MAX GST_RTSP_AUTH_DIGEST

/**
 * GstRTSPHeaderField:
 *
 * Enumeration of rtsp header fields
 */
typedef enum {
  GST_RTSP_HDR_INVALID,

  /*
   * R = Request
   * r = response
   * g = general
   * e = entity
   */
  GST_RTSP_HDR_ACCEPT,              /* Accept               R      opt.      entity */
  GST_RTSP_HDR_ACCEPT_ENCODING,     /* Accept-Encoding      R      opt.      entity */
  GST_RTSP_HDR_ACCEPT_LANGUAGE,     /* Accept-Language      R      opt.      all */
  GST_RTSP_HDR_ALLOW,               /* Allow                r      opt.      all */
  GST_RTSP_HDR_AUTHORIZATION,       /* Authorization        R      opt.      all */
  GST_RTSP_HDR_BANDWIDTH,           /* Bandwidth            R      opt.      all */
  GST_RTSP_HDR_BLOCKSIZE,           /* Blocksize            R      opt.      all but OPTIONS, TEARDOWN */
  GST_RTSP_HDR_CACHE_CONTROL,       /* Cache-Control        g      opt.      SETUP */
  GST_RTSP_HDR_CONFERENCE,          /* Conference           R      opt.      SETUP */
  GST_RTSP_HDR_CONNECTION,          /* Connection           g      req.      all */
  GST_RTSP_HDR_CONTENT_BASE,        /* Content-Base         e      opt.      entity */
  GST_RTSP_HDR_CONTENT_ENCODING,    /* Content-Encoding     e      req.      SET_PARAMETER, DESCRIBE, ANNOUNCE */
  GST_RTSP_HDR_CONTENT_LANGUAGE,    /* Content-Language     e      req.      DESCRIBE, ANNOUNCE */
  GST_RTSP_HDR_CONTENT_LENGTH,      /* Content-Length       e      req.      SET_PARAMETER, ANNOUNCE, entity */
  GST_RTSP_HDR_CONTENT_LOCATION,    /* Content-Location     e      opt.      entity */
  GST_RTSP_HDR_CONTENT_TYPE,        /* Content-Type         e      req.      SET_PARAMETER, ANNOUNCE, entity */
  GST_RTSP_HDR_CSEQ,                /* CSeq                 g      req.      all */
  GST_RTSP_HDR_DATE,                /* Date                 g      opt.      all */
  GST_RTSP_HDR_EXPIRES,             /* Expires              e      opt.      DESCRIBE, ANNOUNCE */
  GST_RTSP_HDR_FROM,                /* From                 R      opt.      all */
  GST_RTSP_HDR_IF_MODIFIED_SINCE,   /* If-Modified-Since    R      opt.      DESCRIBE, SETUP */
  GST_RTSP_HDR_LAST_MODIFIED,       /* Last-Modified        e      opt.      entity */
  GST_RTSP_HDR_PROXY_AUTHENTICATE,  /* Proxy-Authenticate */
  GST_RTSP_HDR_PROXY_REQUIRE,       /* Proxy-Require        R      req.      all */
  GST_RTSP_HDR_PUBLIC,              /* Public               r      opt.      all */
  GST_RTSP_HDR_RANGE,               /* Range                Rr     opt.      PLAY, PAUSE, RECORD */
  GST_RTSP_HDR_REFERER,             /* Referrer              R      opt.      all */
  GST_RTSP_HDR_REQUIRE,             /* Require              R      req.      all */
  GST_RTSP_HDR_RETRY_AFTER,         /* Retry-After          r      opt.      all */
  GST_RTSP_HDR_RTP_INFO,            /* RTP-Info             r      req.      PLAY */
  GST_RTSP_HDR_SCALE,               /* Scale                Rr     opt.      PLAY, RECORD */
  GST_RTSP_HDR_SESSION,             /* Session              Rr     req.      all but SETUP, OPTIONS */
  GST_RTSP_HDR_SERVER,              /* Server               r      opt.      all */
  GST_RTSP_HDR_SPEED,               /* Speed                Rr     opt.      PLAY */
  GST_RTSP_HDR_TRANSPORT,           /* Transport            Rr     req.      SETUP */
  GST_RTSP_HDR_UNSUPPORTED,         /* Unsupported          r      req.      all */
  GST_RTSP_HDR_USER_AGENT,          /* User-Agent           R      opt.      all */
  GST_RTSP_HDR_VIA,                 /* Via                  g      opt.      all */
  GST_RTSP_HDR_WWW_AUTHENTICATE,    /* WWW-Authenticate     r      opt.      all */

  /* Real extensions */
  GST_RTSP_HDR_CLIENT_CHALLENGE,    /* ClientChallenge */
  GST_RTSP_HDR_REAL_CHALLENGE1,     /* RealChallenge1 */
  GST_RTSP_HDR_REAL_CHALLENGE2,     /* RealChallenge2 */
  GST_RTSP_HDR_REAL_CHALLENGE3,     /* RealChallenge3 */
  GST_RTSP_HDR_SUBSCRIBE,           /* Subscribe */
  GST_RTSP_HDR_ALERT,               /* Alert */
  GST_RTSP_HDR_CLIENT_ID,           /* ClientID */
  GST_RTSP_HDR_COMPANY_ID,          /* CompanyID */
  GST_RTSP_HDR_GUID,                /* GUID */
  GST_RTSP_HDR_REGION_DATA,         /* RegionData */
  GST_RTSP_HDR_MAX_ASM_WIDTH,       /* SupportsMaximumASMBandwidth */
  GST_RTSP_HDR_LANGUAGE,            /* Language */
  GST_RTSP_HDR_PLAYER_START_TIME,   /* PlayerStarttime */

  GST_RTSP_HDR_LOCATION,            /* Location */

  GST_RTSP_HDR_ETAG,                /* ETag */
  GST_RTSP_HDR_IF_MATCH,            /* If-Match */

  /* WM extensions [MS-RTSP] */
  GST_RTSP_HDR_ACCEPT_CHARSET,      /* Accept-Charset */
  GST_RTSP_HDR_SUPPORTED,           /* Supported */
  GST_RTSP_HDR_VARY,                /* Vary */
  GST_RTSP_HDR_X_ACCELERATE_STREAMING,    /* X-Accelerate-Streaming */
  GST_RTSP_HDR_X_ACCEPT_AUTHENT,    /* X-Accept-Authentication */
  GST_RTSP_HDR_X_ACCEPT_PROXY_AUTHENT,    /* X-Accept-Proxy-Authentication */
  GST_RTSP_HDR_X_BROADCAST_ID,      /* X-Broadcast-Id */
  GST_RTSP_HDR_X_BURST_STREAMING,   /* X-Burst-Streaming */
  GST_RTSP_HDR_X_NOTICE,            /* X-Notice */
  GST_RTSP_HDR_X_PLAYER_LAG_TIME,   /* X-Player-Lag-Time */
  GST_RTSP_HDR_X_PLAYLIST,          /* X-Playlist */
  GST_RTSP_HDR_X_PLAYLIST_CHANGE_NOTICE,  /* X-Playlist-Change-Notice */
  GST_RTSP_HDR_X_PLAYLIST_GEN_ID,   /* X-Playlist-Gen-Id */
  GST_RTSP_HDR_X_PLAYLIST_SEEK_ID,  /* X-Playlist-Seek-Id */
  GST_RTSP_HDR_X_PROXY_CLIENT_AGENT,      /* X-Proxy-Client-Agent */
  GST_RTSP_HDR_X_PROXY_CLIENT_VERB, /* X-Proxy-Client-Verb */
  GST_RTSP_HDR_X_RECEDING_PLAYLISTCHANGE, /* X-Receding-PlaylistChange */
  GST_RTSP_HDR_X_RTP_INFO,          /* X-RTP-Info */
  GST_RTSP_HDR_X_STARTUPPROFILE,    /* X-StartupProfile */

  GST_RTSP_HDR_TIMESTAMP,           /* Timestamp */

  GST_RTSP_HDR_AUTHENTICATION_INFO, /* Authentication-Info */
  GST_RTSP_HDR_HOST,                /* Host */
  GST_RTSP_HDR_PRAGMA,              /* Pragma */
  GST_RTSP_HDR_X_SERVER_IP_ADDRESS, /* X-Server-IP-Address */
  GST_RTSP_HDR_X_SESSIONCOOKIE,     /* x-sessioncookie */

  GST_RTSP_HDR_RTCP_INTERVAL,       /* RTCP-Interval */

  /* Since: 1.4 */
  GST_RTSP_HDR_KEYMGMT,             /* KeyMgmt */

  /* Since: 1.14 */
  GST_RTSP_HDR_PIPELINED_REQUESTS,  /* Pipelined-Requests Rr     opt.      SETUP */
  GST_RTSP_HDR_MEDIA_PROPERTIES,    /* Media-Properties   Rr     opt.      SETUP */
  GST_RTSP_HDR_SEEK_STYLE,          /* Seek-Type          Rr     opt.      PLAY */
  GST_RTSP_HDR_ACCEPT_RANGES,       /* Accept-Ranges      Rr     opt.      SETUP, GET_PARAMETER */

  /* Onvif extensions */
  GST_RTSP_HDR_FRAMES,              /* Frames             Rr     opt.      PLAY */
  GST_RTSP_HDR_RATE_CONTROL,        /* Rate-control       Rr     opt.      PLAY */

  GST_RTSP_HDR_LAST
} GstRTSPHeaderField;

/**
 * GstRTSPStatusCode:
 *
 * Enumeration of rtsp status codes
 */
typedef enum {
  GST_RTSP_STS_INVALID                              = 0,
  GST_RTSP_STS_CONTINUE                             = 100,
  GST_RTSP_STS_OK                                   = 200,
  GST_RTSP_STS_CREATED                              = 201,
  GST_RTSP_STS_LOW_ON_STORAGE                       = 250,
  GST_RTSP_STS_MULTIPLE_CHOICES                     = 300,
  GST_RTSP_STS_MOVED_PERMANENTLY                    = 301,
  GST_RTSP_STS_MOVE_TEMPORARILY                     = 302,
  GST_RTSP_STS_SEE_OTHER                            = 303,
  GST_RTSP_STS_NOT_MODIFIED                         = 304,
  GST_RTSP_STS_USE_PROXY                            = 305,
  GST_RTSP_STS_BAD_REQUEST                          = 400,
  GST_RTSP_STS_UNAUTHORIZED                         = 401,
  GST_RTSP_STS_PAYMENT_REQUIRED                     = 402,
  GST_RTSP_STS_FORBIDDEN                            = 403,
  GST_RTSP_STS_NOT_FOUND                            = 404,
  GST_RTSP_STS_METHOD_NOT_ALLOWED                   = 405,
  GST_RTSP_STS_NOT_ACCEPTABLE                       = 406,
  GST_RTSP_STS_PROXY_AUTH_REQUIRED                  = 407,
  GST_RTSP_STS_REQUEST_TIMEOUT                      = 408,
  GST_RTSP_STS_GONE                                 = 410,
  GST_RTSP_STS_LENGTH_REQUIRED                      = 411,
  GST_RTSP_STS_PRECONDITION_FAILED                  = 412,
  GST_RTSP_STS_REQUEST_ENTITY_TOO_LARGE             = 413,
  GST_RTSP_STS_REQUEST_URI_TOO_LARGE                = 414,
  GST_RTSP_STS_UNSUPPORTED_MEDIA_TYPE               = 415,
  GST_RTSP_STS_PARAMETER_NOT_UNDERSTOOD             = 451,
  GST_RTSP_STS_CONFERENCE_NOT_FOUND                 = 452,
  GST_RTSP_STS_NOT_ENOUGH_BANDWIDTH                 = 453,
  GST_RTSP_STS_SESSION_NOT_FOUND                    = 454,
  GST_RTSP_STS_METHOD_NOT_VALID_IN_THIS_STATE       = 455,
  GST_RTSP_STS_HEADER_FIELD_NOT_VALID_FOR_RESOURCE  = 456,
  GST_RTSP_STS_INVALID_RANGE                        = 457,
  GST_RTSP_STS_PARAMETER_IS_READONLY                = 458,
  GST_RTSP_STS_AGGREGATE_OPERATION_NOT_ALLOWED      = 459,
  GST_RTSP_STS_ONLY_AGGREGATE_OPERATION_ALLOWED     = 460,
  GST_RTSP_STS_UNSUPPORTED_TRANSPORT                = 461,
  GST_RTSP_STS_DESTINATION_UNREACHABLE              = 462,
  GST_RTSP_STS_KEY_MANAGEMENT_FAILURE               = 463, /* since 1.4 */
  GST_RTSP_STS_INTERNAL_SERVER_ERROR                = 500,
  GST_RTSP_STS_NOT_IMPLEMENTED                      = 501,
  GST_RTSP_STS_BAD_GATEWAY                          = 502,
  GST_RTSP_STS_SERVICE_UNAVAILABLE                  = 503,
  GST_RTSP_STS_GATEWAY_TIMEOUT                      = 504,
  GST_RTSP_STS_RTSP_VERSION_NOT_SUPPORTED           = 505,
  GST_RTSP_STS_OPTION_NOT_SUPPORTED                 = 551
} GstRTSPStatusCode;

GST_RTSP_API
gchar*             gst_rtsp_strresult          (GstRTSPResult result);

GST_RTSP_API
const gchar*       gst_rtsp_method_as_text     (GstRTSPMethod method);

GST_RTSP_API
const gchar*       gst_rtsp_version_as_text    (GstRTSPVersion version);

GST_RTSP_API
const gchar*       gst_rtsp_header_as_text     (GstRTSPHeaderField field);

GST_RTSP_API
const gchar*       gst_rtsp_status_as_text     (GstRTSPStatusCode code);

GST_RTSP_API
gchar*             gst_rtsp_options_as_text    (GstRTSPMethod options);

GST_RTSP_API
GstRTSPMethod      gst_rtsp_options_from_text  (const gchar *options);

GST_RTSP_API
GstRTSPHeaderField gst_rtsp_find_header_field  (const gchar *header);

GST_RTSP_API
GstRTSPMethod      gst_rtsp_find_method        (const gchar *method);

GST_RTSP_API
gboolean           gst_rtsp_header_allow_multiple (GstRTSPHeaderField field);

GST_RTSP_API
gchar *            gst_rtsp_generate_digest_auth_response (const gchar *algorithm,
                                                           const gchar *method,
                                                           const gchar *realm,
                                                           const gchar *username,
                                                           const gchar *password,
                                                           const gchar *uri,
                                                           const gchar *nonce);

GST_RTSP_API
gchar *            gst_rtsp_generate_digest_auth_response_from_md5 (const gchar *algorithm,
                                                                    const gchar * method,
                                                                    const gchar * md5,
                                                                    const gchar * uri,
                                                                    const gchar * nonce);

G_END_DECLS

#endif /* __GST_RTSP_DEFS_H__ */
