/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _PKCS11NI_H_
#define _PKCS11NI_H_

/*
 * pkcs11ni.h
 *
 * This file contains softoken private exports for NSS
 */

/* softoken slot ID's */
#define SFTK_MIN_USER_SLOT_ID 4
#define SFTK_MAX_USER_SLOT_ID 100
#define SFTK_MIN_FIPS_USER_SLOT_ID 101
#define SFTK_MAX_FIPS_USER_SLOT_ID 127

#endif /* _PKCS11NI_H_ */
