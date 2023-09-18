/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef __ec_h_
#define __ec_h_

#define EC_DEBUG 0

#define ANSI_X962_CURVE_OID_TOTAL_LEN 10
#define SECG_CURVE_OID_TOTAL_LEN 7
#define PKIX_NEWCURVES_OID_TOTAL_LEN 11

struct ECMethodStr {
    ECCurveName name;
    SECStatus (*mul)(SECItem *result, SECItem *scalar, SECItem *point);
    SECStatus (*validate)(const SECItem *point);
};
typedef struct ECMethodStr ECMethod;

#endif /* __ec_h_ */
