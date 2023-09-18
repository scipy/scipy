/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _HMACCT_H_
#define _HMACCT_H_

SEC_BEGIN_PROTOS

extern SECStatus HMAC_ConstantTime(
    unsigned char *result,
    unsigned int *resultLen,
    unsigned int maxResultLen,
    const SECHashObject *hashObj,
    const unsigned char *secret,
    unsigned int secretLen,
    const unsigned char *header,
    unsigned int headerLen,
    const unsigned char *body,
    unsigned int bodyLen,
    unsigned int bodyTotalLen);

extern SECStatus SSLv3_MAC_ConstantTime(
    unsigned char *result,
    unsigned int *resultLen,
    unsigned int maxResultLen,
    const SECHashObject *hashObj,
    const unsigned char *secret,
    unsigned int secretLen,
    const unsigned char *header,
    unsigned int headerLen,
    const unsigned char *body,
    unsigned int bodyLen,
    unsigned int bodyTotalLen);

SEC_END_PROTOS

#endif
