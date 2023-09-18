/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _SHSIGN_H_
#define _SHSIGN_H_

#define SGN_SUFFIX ".chk"
#define NSS_SIGN_CHK_MAGIC1 0xf1
#define NSS_SIGN_CHK_MAGIC2 0xc5
/* new hmac based signatures */
#define NSS_SIGN_CHK_MAJOR_VERSION 0x02
#define NSS_SIGN_CHK_MINOR_VERSION 0x01
#define NSS_SIGN_CHK_TYPE_FLAGS 0xff000000
#define NSS_SIGN_CHK_FLAG_HMAC 0x80000000

typedef struct NSSSignChkHeaderStr NSSSignChkHeader;
struct NSSSignChkHeaderStr {
    unsigned char magic1;
    unsigned char magic2;
    unsigned char majorVersion;
    unsigned char minorVersion;
    unsigned char offset[4];
    unsigned char type[4];
};
#endif /* _SHSIGN_H_ */
