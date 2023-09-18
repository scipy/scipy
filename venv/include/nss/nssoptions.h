/*
 * NSS utility functions
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 *  Include the default limits here
 */
/* SSL default limits are here so we don't have to import a private SSL header
 * file into NSS proper */

/* The minimum server key sizes accepted by the clients.
 * Not 1024 to be conservative. */
#define SSL_RSA_MIN_MODULUS_BITS 1023
/* 1023 to avoid cases where p = 2q+1 for a 512-bit q turns out to be
 * only 1023 bits and similar.  We don't have good data on whether this
 * happens because NSS used to count bit lengths incorrectly. */
#define SSL_DH_MIN_P_BITS 1023
#define SSL_DSA_MIN_P_BITS 1023
/* not really used by SSL, but define it here for consistency */
#define SSL_ECC_MIN_CURVE_BITS 255
