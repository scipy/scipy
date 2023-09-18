/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/* Although this is not an exported header file, code which uses elliptic
 * curve point operations will need to include it. */

#ifndef __ecl_h_
#define __ecl_h_

#include "blapi.h"
#include "ecl-exp.h"
#include "mpi.h"
#include "eclt.h"

struct ECGroupStr;
typedef struct ECGroupStr ECGroup;

/* Construct ECGroup from named parameters. */
ECGroup *ECGroup_fromName(const ECCurveName name);

/* Free an allocated ECGroup. */
void ECGroup_free(ECGroup *group);

/* Elliptic curve scalar-point multiplication. Computes Q(x, y) = k * P(x,
 * y).  If x, y = NULL, then P is assumed to be the generator (base point)
 * of the group of points on the elliptic curve. Input and output values
 * are assumed to be NOT field-encoded. */
mp_err ECPoint_mul(const ECGroup *group, const mp_int *k, const mp_int *px,
                   const mp_int *py, mp_int *qx, mp_int *qy);

/* Elliptic curve scalar-point multiplication. Computes Q(x, y) = k1 * G +
 * k2 * P(x, y), where G is the generator (base point) of the group of
 * points on the elliptic curve. Input and output values are assumed to
 * be NOT field-encoded. */
mp_err ECPoints_mul(const ECGroup *group, const mp_int *k1,
                    const mp_int *k2, const mp_int *px, const mp_int *py,
                    mp_int *qx, mp_int *qy);

/* Validates an EC public key as described in Section 5.2.2 of X9.62.
 * Returns MP_YES if the public key is valid, MP_NO if the public key
 * is invalid, or an error code if the validation could not be
 * performed. */
mp_err ECPoint_validate(const ECGroup *group, const mp_int *px, const mp_int *py);

SECStatus ec_Curve25519_pt_mul(SECItem *X, SECItem *k, SECItem *P);
SECStatus ec_Curve25519_pt_validate(const SECItem *px);

#endif /* __ecl_h_ */
