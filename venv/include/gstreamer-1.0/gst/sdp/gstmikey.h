/* GStreamer
 * Copyright (C) <2014> Wim Taymans <wim.taymans@gmail.com>
 *
 * gstmikey.h: various helper functions to manipulate mikey messages
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

#ifndef __GST_MIKEY_H__
#define __GST_MIKEY_H__

#include <gst/gst.h>
#include <gst/sdp/sdp-prelude.h>

G_BEGIN_DECLS

GST_SDP_API
GType gst_mikey_message_get_type(void);
#define GST_TYPE_MIKEY_MESSAGE (gst_mikey_message_get_type())

typedef struct _GstMIKEYMessage GstMIKEYMessage;
typedef struct _GstMIKEYEncryptInfo GstMIKEYEncryptInfo;
typedef struct _GstMIKEYDecryptInfo GstMIKEYDecryptInfo;

/**
 * GST_MIKEY_VERSION:
 *
 * The supported MIKEY version 1.
 */
#define GST_MIKEY_VERSION 1

/**
 * GstMIKEYType:
 * @GST_MIKEY_TYPE_INVALID: Invalid type
 * @GST_MIKEY_TYPE_PSK_INIT: Initiator's pre-shared key message
 * @GST_MIKEY_TYPE_PSK_VERIFY: Verification message of a Pre-shared key message
 * @GST_MIKEY_TYPE_PK_INIT: Initiator's public-key transport message
 * @GST_MIKEY_TYPE_PK_VERIFY: Verification message of a public-key message
 * @GST_MIKEY_TYPE_DH_INIT: Initiator's DH exchange message
 * @GST_MIKEY_TYPE_DH_RESP: Responder's DH exchange message
 * @GST_MIKEY_TYPE_ERROR: Error message
 *
 * Different MIKEY data types.
 */
typedef enum
{
  GST_MIKEY_TYPE_INVALID    = -1,
  GST_MIKEY_TYPE_PSK_INIT   = 0,
  GST_MIKEY_TYPE_PSK_VERIFY = 1,
  GST_MIKEY_TYPE_PK_INIT    = 2,
  GST_MIKEY_TYPE_PK_VERIFY  = 3,
  GST_MIKEY_TYPE_DH_INIT    = 4,
  GST_MIKEY_TYPE_DH_RESP    = 5,
  GST_MIKEY_TYPE_ERROR      = 6
} GstMIKEYType;

/**
 * GstMIKEYPayloadType:
 * @GST_MIKEY_PT_LAST: Last payload
 * @GST_MIKEY_PT_KEMAC: Key data transport payload
 * @GST_MIKEY_PT_PKE: Envelope data payload
 * @GST_MIKEY_PT_DH: DH data payload
 * @GST_MIKEY_PT_SIGN: Signature payload
 * @GST_MIKEY_PT_T: Timestamp payload
 * @GST_MIKEY_PT_ID: ID payload
 * @GST_MIKEY_PT_CERT: Certificate Payload
 * @GST_MIKEY_PT_CHASH: Cert hash payload
 * @GST_MIKEY_PT_V: Verification message payload
 * @GST_MIKEY_PT_SP: Security Policy payload
 * @GST_MIKEY_PT_RAND: RAND payload
 * @GST_MIKEY_PT_ERR: Error payload
 * @GST_MIKEY_PT_KEY_DATA: Key data sub-payload
 * @GST_MIKEY_PT_GEN_EXT: General Extension Payload

 * Different MIKEY Payload types.
 */
typedef enum
{
  GST_MIKEY_PT_LAST      = 0,
  GST_MIKEY_PT_KEMAC     = 1,
  GST_MIKEY_PT_PKE       = 2,
  GST_MIKEY_PT_DH        = 3,
  GST_MIKEY_PT_SIGN      = 4,
  GST_MIKEY_PT_T         = 5,
  GST_MIKEY_PT_ID        = 6,
  GST_MIKEY_PT_CERT      = 7,
  GST_MIKEY_PT_CHASH     = 8,
  GST_MIKEY_PT_V         = 9,
  GST_MIKEY_PT_SP        = 10,
  GST_MIKEY_PT_RAND      = 11,
  GST_MIKEY_PT_ERR       = 12,
  GST_MIKEY_PT_KEY_DATA  = 20,
  GST_MIKEY_PT_GEN_EXT   = 21
} GstMIKEYPayloadType;

/**
 * GstMIKEYPRFFunc:
 * @GST_MIKEY_PRF_MIKEY_1: MIKEY-1 PRF function
 *
 * The PRF function that has been/will be used for key derivation
 */
typedef enum
{
  GST_MIKEY_PRF_MIKEY_1  = 0
} GstMIKEYPRFFunc;

/**
 * GstMIKEYMapType:
 * @GST_MIKEY_MAP_TYPE_SRTP: SRTP
 *
 * Specifies the method of uniquely mapping Crypto Sessions to the security
 * protocol sessions.
 */
typedef enum
{
  GST_MIKEY_MAP_TYPE_SRTP  = 0
} GstMIKEYMapType;

/**
 * GstMIKEYMapSRTP:
 * @policy: The security policy applied for the stream with @ssrc
 * @ssrc: the SSRC that must be used for the stream
 * @roc: current rollover counter
 *
 * The Security policy Map item for SRTP
 */
typedef struct {
  guint8  policy;
  guint32 ssrc;
  guint32 roc;
} GstMIKEYMapSRTP;

typedef struct _GstMIKEYPayload GstMIKEYPayload;

GST_SDP_API
GType gst_mikey_payload_get_type(void);
#define GST_TYPE_MIKEY_PAYLOAD (gst_mikey_payload_get_type())

/**
 * GstMIKEYPayload:
 * @type: the payload type
 * @len: length of the payload
 *
 * Hold the common fields for all payloads
 */
struct _GstMIKEYPayload {
  /* < private > */
  GstMiniObject mini_object;

  /* < public > */
  GstMIKEYPayloadType type;
  guint len;
};

GST_SDP_API
GstMIKEYPayload *   gst_mikey_payload_new      (GstMIKEYPayloadType type);

/**
 * gst_mikey_payload_ref:
 * @payload: The payload to refcount
 *
 * Increase the refcount of this payload.
 *
 * Returns: (transfer full): @payload (for convenience when doing assignments)
 *
 * Since: 1.4
 */
static inline GstMIKEYPayload *
gst_mikey_payload_ref (GstMIKEYPayload * payload)
{
  return (GstMIKEYPayload *) gst_mini_object_ref (GST_MINI_OBJECT_CAST (payload));
}

/**
 * gst_mikey_payload_unref:
 * @payload: (transfer full): the payload to refcount
 *
 * Decrease the refcount of an payload, freeing it if the refcount reaches 0.
 *
 * Since: 1.4
 */
static inline void
gst_mikey_payload_unref (GstMIKEYPayload * payload)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (payload));
}

/**
 * gst_mikey_payload_copy:
 * @payload: a #GstMIKEYPayload.
 *
 * Create a copy of the given payload.
 *
 * Returns: (transfer full): a new copy of @payload.
 *
 * Since: 1.4
 */
static inline GstMIKEYPayload *
gst_mikey_payload_copy (const GstMIKEYPayload * payload)
{
  return (GstMIKEYPayload *) gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (payload));
}

/**
 * GstMIKEYEncAlg:
 * @GST_MIKEY_ENC_NULL: no encryption
 * @GST_MIKEY_ENC_AES_CM_128: AES-CM using a 128-bit key
 * @GST_MIKEY_ENC_AES_KW_128: AES Key Wrap using a 128-bit key
 * @GST_MIKEY_ENC_AES_GCM_128: AES-GCM using a 128-bit key (Since: 1.16)
 *
 * The encryption algorithm used to encrypt the Encr data field
 */
typedef enum
{
  GST_MIKEY_ENC_NULL        = 0,
  GST_MIKEY_ENC_AES_CM_128  = 1,
  GST_MIKEY_ENC_AES_KW_128  = 2,
  GST_MIKEY_ENC_AES_GCM_128 = 6
} GstMIKEYEncAlg;

/**
 * GstMIKEYMacAlg:
 * @GST_MIKEY_MAC_NULL: no authentication
 * @GST_MIKEY_MAC_HMAC_SHA_1_160: HMAC-SHA-1-160
 *
 * Specifies the authentication algorithm used
 */
typedef enum
{
  GST_MIKEY_MAC_NULL            = 0,
  GST_MIKEY_MAC_HMAC_SHA_1_160  = 1
} GstMIKEYMacAlg;

/**
 * GstMIKEYPayloadKEMAC:
 * @pt: the common #GstMIKEYPayload
 * @enc_alg: the #GstMIKEYEncAlg
 * @mac_alg: the #GstMIKEYMacAlg
 * @subpayloads: the subpayloads
 *
 * A structure holding the KEMAC payload
 */
typedef struct {
  GstMIKEYPayload pt;

  GstMIKEYEncAlg  enc_alg;
  GstMIKEYMacAlg  mac_alg;
  GArray *subpayloads;
} GstMIKEYPayloadKEMAC;

GST_SDP_API
gboolean                gst_mikey_payload_kemac_set        (GstMIKEYPayload *payload,
                                                            GstMIKEYEncAlg enc_alg,
                                                            GstMIKEYMacAlg mac_alg);

GST_SDP_API
guint                   gst_mikey_payload_kemac_get_n_sub  (const GstMIKEYPayload *payload);

GST_SDP_API
const GstMIKEYPayload * gst_mikey_payload_kemac_get_sub    (const GstMIKEYPayload *payload, guint idx);

GST_SDP_API
gboolean                gst_mikey_payload_kemac_remove_sub (GstMIKEYPayload *payload, guint idx);

GST_SDP_API
gboolean                gst_mikey_payload_kemac_add_sub    (GstMIKEYPayload *payload,
                                                            GstMIKEYPayload *newpay);

/**
 * GstMIKEYCacheType:
 * @GST_MIKEY_CACHE_NONE: The envelope key MUST NOT be cached
 * @GST_MIKEY_CACHE_ALWAYS: The envelope key MUST be cached
 * @GST_MIKEY_CACHE_FOR_CSB: The envelope key MUST be cached, but only
 *                           to be used for the specific CSB.
 *
 * The different cache types
 */
typedef enum
{
  GST_MIKEY_CACHE_NONE       = 0,
  GST_MIKEY_CACHE_ALWAYS     = 1,
  GST_MIKEY_CACHE_FOR_CSB    = 2
} GstMIKEYCacheType;

/**
 * GstMIKEYPayloadPKE:
 * @pt: the common #GstMIKEYPayload
 * @C: envelope key cache indicator
 * @data_len: length of @data
 * @data: the encrypted envelope key
 *
 * The Envelope data payload contains the encrypted envelope key that is
 * used in the public-key transport to protect the data in the Key data
 * transport payload.  The encryption algorithm used is implicit from
 * the certificate/public key used.
 */
typedef struct {
  GstMIKEYPayload pt;

  GstMIKEYCacheType C;
  guint16           data_len;
  guint8           *data;
} GstMIKEYPayloadPKE;

GST_SDP_API
gboolean               gst_mikey_payload_pke_set     (GstMIKEYPayload *payload,
                                                      GstMIKEYCacheType C,
                                                      guint16 data_len, const guint8 *data);


/**
 * GstMIKEYTSType:
 * @GST_MIKEY_TS_TYPE_NTP_UTC: an NTP time in UTC timezone
 * @GST_MIKEY_TS_TYPE_NTP: an NTP time
 * @GST_MIKEY_TS_TYPE_COUNTER: a counter
 *
 * Specifies the timestamp type.
 */
typedef enum
{
  GST_MIKEY_TS_TYPE_NTP_UTC  = 0,
  GST_MIKEY_TS_TYPE_NTP      = 1,
  GST_MIKEY_TS_TYPE_COUNTER  = 2
} GstMIKEYTSType;

/**
 * GstMIKEYPayloadT:
 * @pt: the payload header
 * @type: a #GstMIKEYTSType
 * @ts_value: the timestamp value
 *
 * The timestamp payload carries the timestamp information
 */
typedef struct {
  GstMIKEYPayload pt;

  GstMIKEYTSType  type;
  guint8         *ts_value;
} GstMIKEYPayloadT;

GST_SDP_API
gboolean   gst_mikey_payload_t_set   (GstMIKEYPayload *payload,
                                      GstMIKEYTSType type, const guint8 *ts_value);

/**
 * GstMIKEYPayloadSPParam:
 * @type: specifies the type of the parameter
 * @len: specifies the length of @val
 * @val: specifies the value of the parameter
 *
 * A Type/Length/Value field for security parameters
 */
typedef struct {
  guint8  type;
  guint8  len;
  guint8 *val;
} GstMIKEYPayloadSPParam;

/**
 * GstMIKEYSecProto:
 * @GST_MIKEY_SEC_PROTO_SRTP: SRTP
 *
 * Specifies the security protocol
 */
typedef enum
{
  GST_MIKEY_SEC_PROTO_SRTP  = 0
} GstMIKEYSecProto;

/**
 * GstMIKEYSecSRTP:
 * @GST_MIKEY_SP_SRTP_ENC_ALG: Encryption algorithm
 * @GST_MIKEY_SP_SRTP_ENC_KEY_LEN: Session Encr. key length
 * @GST_MIKEY_SP_SRTP_AUTH_ALG: Authentication algorithm
 * @GST_MIKEY_SP_SRTP_AUTH_KEY_LEN: Session Auth. key length
 * @GST_MIKEY_SP_SRTP_SALT_KEY_LEN: Session Salt key length
 * @GST_MIKEY_SP_SRTP_PRF: SRTP Pseudo Random Function
 * @GST_MIKEY_SP_SRTP_KEY_DERIV_RATE: Key derivation rate
 * @GST_MIKEY_SP_SRTP_SRTP_ENC: SRTP encryption off/on, 0 if off, 1 if on
 * @GST_MIKEY_SP_SRTP_SRTCP_ENC: SRTCP encryption off/on, 0 if off, 1 if on
 * @GST_MIKEY_SP_SRTP_FEC_ORDER: sender's FEC order
 * @GST_MIKEY_SP_SRTP_SRTP_AUTH: SRTP authentication off/on, 0 if off, 1 if on
 * @GST_MIKEY_SP_SRTP_AUTH_TAG_LEN: Authentication tag length
 * @GST_MIKEY_SP_SRTP_SRTP_PREFIX_LEN: SRTP prefix length
 * @GST_MIKEY_SP_SRTP_AEAD_AUTH_TAG_LEN: AEAD authentication tag length (Since: 1.16)
 *
 * This policy specifies the parameters for SRTP and SRTCP
 */
typedef enum
{
  GST_MIKEY_SP_SRTP_ENC_ALG         =    0,
  GST_MIKEY_SP_SRTP_ENC_KEY_LEN     =    1,
  GST_MIKEY_SP_SRTP_AUTH_ALG        =    2,
  GST_MIKEY_SP_SRTP_AUTH_KEY_LEN    =    3,
  GST_MIKEY_SP_SRTP_SALT_KEY_LEN    =    4,
  GST_MIKEY_SP_SRTP_PRF             =    5,
  GST_MIKEY_SP_SRTP_KEY_DERIV_RATE  =    6,
  GST_MIKEY_SP_SRTP_SRTP_ENC        =    7,
  GST_MIKEY_SP_SRTP_SRTCP_ENC       =    8,
  GST_MIKEY_SP_SRTP_FEC_ORDER       =    9,
  GST_MIKEY_SP_SRTP_SRTP_AUTH       =   10,
  GST_MIKEY_SP_SRTP_AUTH_TAG_LEN    =   11,
  GST_MIKEY_SP_SRTP_SRTP_PREFIX_LEN =   12,
  GST_MIKEY_SP_SRTP_AEAD_AUTH_TAG_LEN = 20
} GstMIKEYSecSRTP;

/**
 * GstMIKEYPayloadSP:
 * @pt: the payload header
 * @policy: the policy number
 * @proto: the security protocol
 * @params: array of #GstMIKEYPayloadSPParam
 *
 * The Security Policy payload defines a set of policies that apply to a
 * specific security protocol
 */
typedef struct {
  GstMIKEYPayload pt;

  guint policy;
  GstMIKEYSecProto proto;
  GArray *params;
} GstMIKEYPayloadSP;

GST_SDP_API
gboolean            gst_mikey_payload_sp_set          (GstMIKEYPayload *payload,
                                                       guint policy, GstMIKEYSecProto proto);
GST_SDP_API
guint               gst_mikey_payload_sp_get_n_params (const GstMIKEYPayload *payload);

GST_SDP_API
const GstMIKEYPayloadSPParam *
                    gst_mikey_payload_sp_get_param    (const GstMIKEYPayload *payload, guint idx);

GST_SDP_API
gboolean            gst_mikey_payload_sp_remove_param (GstMIKEYPayload *payload, guint idx);

GST_SDP_API
gboolean            gst_mikey_payload_sp_add_param    (GstMIKEYPayload *payload,
                                                       guint8 type, guint8 len, const guint8 *val);

/**
 * GstMIKEYPayloadRAND:
 * @pt: the payload header
 * @len: the length of @rand
 * @rand: random values
 *
 * The RAND payload consists of a (pseudo-)random bit-string
 */
typedef struct {
  GstMIKEYPayload pt;

  guint8  len;
  guint8 *rand;
} GstMIKEYPayloadRAND;

GST_SDP_API
gboolean   gst_mikey_payload_rand_set     (GstMIKEYPayload *payload,
                                           guint8 len, const guint8 *rand);

/**
 * GstMIKEYKeyDataType:
 * @GST_MIKEY_KD_TGK: a TEK Generation Key
 * @GST_MIKEY_KD_TEK: Traffic-Encrypting Key
 *
 * The type of key.
 */
typedef enum
{
  GST_MIKEY_KD_TGK      = 0,
  GST_MIKEY_KD_TEK      = 2,
} GstMIKEYKeyDataType;

/**
 * GstMIKEYKVType:
 * @GST_MIKEY_KV_NULL: No specific usage rule
 * @GST_MIKEY_KV_SPI: The key is associated with the SPI/MKI
 * @GST_MIKEY_KV_INTERVAL: The key has a start and expiration time
 *
 * The key validity type
 */
typedef enum
{
  GST_MIKEY_KV_NULL      = 0,
  GST_MIKEY_KV_SPI       = 1,
  GST_MIKEY_KV_INTERVAL  = 2,
} GstMIKEYKVType;

/**
 * GstMIKEYPayloadKeyData:
 * @pt: the payload header
 * @key_type: the #GstMIKEYKeyDataType of @key_data
 * @key_len: length of @key_data
 * @key_data: the key data
 * @salt_len: the length of @salt_data, can be 0
 * @salt_data: salt data
 * @kv_type: the Key Validity type
 * @kv_len: length of @kv_data
 * @kv_data: key validity data
 *
 * The Key data payload contains key material. It should be added as sub
 * payload to the KEMAC.
 */
typedef struct {
  GstMIKEYPayload pt;

  GstMIKEYKeyDataType key_type;
  guint16  key_len;
  guint8  *key_data;
  guint16  salt_len;
  guint8  *salt_data;
  GstMIKEYKVType kv_type;
  guint8   kv_len[2];
  guint8  *kv_data[2];
} GstMIKEYPayloadKeyData;

GST_SDP_API
gboolean   gst_mikey_payload_key_data_set_key      (GstMIKEYPayload *payload,
                                                    GstMIKEYKeyDataType key_type,
                                                    guint16 key_len, const guint8 *key_data);

GST_SDP_API
gboolean   gst_mikey_payload_key_data_set_salt     (GstMIKEYPayload *payload,
                                                    guint16 salt_len, const guint8 *salt_data);

GST_SDP_API
gboolean   gst_mikey_payload_key_data_set_spi      (GstMIKEYPayload *payload,
                                                    guint8 spi_len, const guint8 *spi_data);

GST_SDP_API
gboolean   gst_mikey_payload_key_data_set_interval (GstMIKEYPayload *payload,
                                                    guint8 vf_len, const guint8 *vf_data,
                                                    guint8 vt_len, const guint8 *vt_data);

/**
 * GstMIKEYMessage:
 * @version: the version
 * @type: the #GstMIKEYType message type
 * @V: verify flag
 * @prf_func: a #GstMIKEYPRFFunc
 * @CSB_id: Identifies the Crypto Session Bundle
 * @map_type: a #GstMIKEYMapType
 * @map_info: map info array of type depending on @map_type
 * @payloads: the payload array of #GstMIKEYPayload
 *
 * Structure holding the information of the MIKEY message
 */
struct _GstMIKEYMessage
{
  /* < private > */
  GstMiniObject mini_object;

  /* < public > */
  guint8 version;
  GstMIKEYType type;
  gboolean V;
  GstMIKEYPRFFunc prf_func;
  guint32 CSB_id;
  GstMIKEYMapType map_type;
  GArray *map_info;
  GArray *payloads;
};


GST_SDP_API
GstMIKEYMessage *           gst_mikey_message_new               (void);

GST_SDP_API
GstMIKEYMessage *           gst_mikey_message_new_from_data     (gconstpointer data, gsize size,
                                                                 GstMIKEYDecryptInfo *info, GError **error);

GST_SDP_API
GstMIKEYMessage *           gst_mikey_message_new_from_bytes    (GBytes *bytes, GstMIKEYDecryptInfo *info,
                                                                 GError **error);

GST_SDP_API
GBytes *                    gst_mikey_message_to_bytes          (GstMIKEYMessage *msg, GstMIKEYEncryptInfo *info,
                                                                 GError **error);

GST_SDP_API
GstMIKEYMessage *           gst_mikey_message_new_from_caps     (GstCaps *caps);

GST_SDP_API
gboolean                    gst_mikey_message_to_caps           (const GstMIKEYMessage *msg, GstCaps *caps);

GST_SDP_API
gchar *                     gst_mikey_message_base64_encode     (GstMIKEYMessage* msg);

/**
 * gst_mikey_message_ref:
 * @message: The message to refcount
 *
 * Increase the refcount of this message.
 *
 * Returns: (transfer full): @message (for convenience when doing assignments)
 *
 * Since: 1.4
 */
static inline GstMIKEYMessage *
gst_mikey_message_ref (GstMIKEYMessage * message)
{
  return (GstMIKEYMessage *) gst_mini_object_ref (GST_MINI_OBJECT_CAST (message));
}

/**
 * gst_mikey_message_unref:
 * @message: (transfer full): the message to refcount
 *
 * Decrease the refcount of an message, freeing it if the refcount reaches 0.
 *
 * Since: 1.4
 */
static inline void
gst_mikey_message_unref (GstMIKEYMessage * message)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (message));
}

/**
 * gst_mikey_message_copy:
 * @message: a #GstMIKEYMessage.
 *
 * Create a copy of the given message.
 *
 * Returns: (transfer full): a new copy of @message.
 *
 * Since: 1.4
 */
static inline GstMIKEYMessage *
gst_mikey_message_copy (const GstMIKEYMessage * message)
{
  return (GstMIKEYMessage *) gst_mini_object_copy (GST_MINI_OBJECT_CONST_CAST (message));
}


GST_SDP_API
gboolean                    gst_mikey_message_set_info          (GstMIKEYMessage *msg,
                                                                 guint8 version, GstMIKEYType type, gboolean V,
                                                                 GstMIKEYPRFFunc prf_func, guint32 CSB_id,
                                                                 GstMIKEYMapType map_type);

GST_SDP_API
guint                       gst_mikey_message_get_n_cs          (const GstMIKEYMessage *msg);

/* SRTP crypto sessions */

GST_SDP_API
const GstMIKEYMapSRTP *     gst_mikey_message_get_cs_srtp       (const GstMIKEYMessage *msg, guint idx);

GST_SDP_API
gboolean                    gst_mikey_message_insert_cs_srtp    (GstMIKEYMessage *msg, gint idx,
                                                                 const GstMIKEYMapSRTP *map);

GST_SDP_API
gboolean                    gst_mikey_message_replace_cs_srtp   (GstMIKEYMessage *msg, gint idx,
                                                                 const GstMIKEYMapSRTP *map);

GST_SDP_API
gboolean                    gst_mikey_message_remove_cs_srtp    (GstMIKEYMessage *msg, gint idx);

GST_SDP_API
gboolean                    gst_mikey_message_add_cs_srtp       (GstMIKEYMessage *msg,
                                                                 guint8 policy, guint32 ssrc, guint32 roc);

/* adding/retrieving payloads */

GST_SDP_API
guint                       gst_mikey_message_get_n_payloads    (const GstMIKEYMessage *msg);

GST_SDP_API
const GstMIKEYPayload *     gst_mikey_message_get_payload       (const GstMIKEYMessage *msg, guint idx);

GST_SDP_API
const GstMIKEYPayload *     gst_mikey_message_find_payload      (const GstMIKEYMessage *msg,
                                                                 GstMIKEYPayloadType type, guint nth);

GST_SDP_API
gboolean                    gst_mikey_message_remove_payload    (GstMIKEYMessage *msg, guint idx);

GST_SDP_API
gboolean                    gst_mikey_message_insert_payload    (GstMIKEYMessage *msg, guint idx,
                                                                 GstMIKEYPayload *payload);

GST_SDP_API
gboolean                    gst_mikey_message_add_payload       (GstMIKEYMessage *msg,
                                                                 GstMIKEYPayload *payload);

GST_SDP_API
gboolean                    gst_mikey_message_replace_payload   (GstMIKEYMessage *msg, guint idx,
                                                                 GstMIKEYPayload *payload);


/* Key data transport payload (KEMAC) */
/* Envelope data payload (PKE) */

GST_SDP_API
gboolean                    gst_mikey_message_add_pke           (GstMIKEYMessage *msg,
                                                                 GstMIKEYCacheType C,
                                                                 guint16 data_len, const guint8 *data);
/* DH data payload (DH) */
/* Signature payload (SIGN) */

/* Timestamp payload (T) */

GST_SDP_API
gboolean                    gst_mikey_message_add_t             (GstMIKEYMessage *msg,
                                                                 GstMIKEYTSType type, const guint8 *ts_value);

GST_SDP_API
gboolean                    gst_mikey_message_add_t_now_ntp_utc (GstMIKEYMessage *msg);
/* ID payload (ID) */
/* Certificate Payload (CERT) */
/* Cert hash payload (CHASH)*/
/* Ver msg payload (V) */
/* Security Policy payload (SP)*/
/* RAND payload (RAND) */

GST_SDP_API
gboolean                    gst_mikey_message_add_rand          (GstMIKEYMessage *msg,
                                                                 guint8 len, const guint8 *rand);

GST_SDP_API
gboolean                    gst_mikey_message_add_rand_len      (GstMIKEYMessage *msg, guint8 len);

/* Error payload (ERR) */
/* Key data sub-payload */
/* General Extension Payload */


G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstMIKEYMessage, gst_mikey_message_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstMIKEYPayload, gst_mikey_payload_unref)

G_END_DECLS

#endif /* __GST_MIKEY_H__ */
