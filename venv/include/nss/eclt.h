/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/* This header holds ECC types and must not be exported publicly. */

#ifndef __eclt_h_
#define __eclt_h_

/* byte encoding of curve parameters */
struct ECCurveBytesStr {
    char *text;
    ECField field;
    size_t size;
    const PRUint8 *irr;
    const PRUint8 *curvea;
    const PRUint8 *curveb;
    const PRUint8 *genx;
    const PRUint8 *geny;
    const PRUint8 *order;
    const PRUint8 *base;
    int cofactor;
    int security;
    size_t pointSize;
    size_t scalarSize;
    unsigned int usage;
};
typedef struct ECCurveBytesStr ECCurveBytes;

#endif /* __ecl_h_ */
