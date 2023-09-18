/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "mpi.h"

#define CHECK_SEC_OK(func)         \
    if (SECSuccess != (rv = func)) \
    goto cleanup

#define CHECK_MPI_OK(func)      \
    if (MP_OKAY > (err = func)) \
    goto cleanup

#define OCTETS_TO_MPINT(oc, mp, len) \
    CHECK_MPI_OK(mp_read_unsigned_octets((mp), oc, len))

#define SECITEM_TO_MPINT(it, mp) \
    CHECK_MPI_OK(mp_read_unsigned_octets((mp), (it).data, (it).len))

#define MPINT_TO_SECITEM(mp, it, arena)                         \
    do {                                                        \
        int mpintLen = mp_unsigned_octet_size(mp);              \
        if (mpintLen <= 0) {                                    \
            err = MP_RANGE;                                     \
            goto cleanup;                                       \
        }                                                       \
        SECITEM_AllocItem(arena, (it), mpintLen);               \
        if ((it)->data == NULL) {                               \
            err = MP_MEM;                                       \
            goto cleanup;                                       \
        }                                                       \
        err = mp_to_unsigned_octets(mp, (it)->data, (it)->len); \
        if (err < 0)                                            \
            goto cleanup;                                       \
        else                                                    \
            err = MP_OKAY;                                      \
    } while (0)

#define MP_TO_SEC_ERROR(err)                          \
    switch (err) {                                    \
        case MP_MEM:                                  \
            PORT_SetError(SEC_ERROR_NO_MEMORY);       \
            break;                                    \
        case MP_RANGE:                                \
            PORT_SetError(SEC_ERROR_BAD_DATA);        \
            break;                                    \
        case MP_BADARG:                               \
            PORT_SetError(SEC_ERROR_INVALID_ARGS);    \
            break;                                    \
        default:                                      \
            PORT_SetError(SEC_ERROR_LIBRARY_FAILURE); \
            break;                                    \
    }

/* Fill the `used` digits of an mp_int with random bits */
mp_err mpp_random_secure(mp_int *a);

/* Pseudo-primality testing using `mpp_random_secure` to choose Miller-Rabin base */
mp_err mpp_pprime_secure(mp_int *a, int nt);

/* Variant of `mpp_make_prime` using `mpp_random_secure` to choose Miller-Rabin base */
mp_err mpp_make_prime_secure(mp_int *start, mp_size nBits, mp_size strong);
