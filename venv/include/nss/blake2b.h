/*
 * blake2b.h - header file for blake2b hash function
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BLAKE_H
#define BLAKE_H

#include <stddef.h>
#include <stdint.h>

struct Blake2bContextStr {
    uint64_t h[8];                     /* chained state */
    uint64_t t[2];                     /* total number of bytes */
    uint64_t f;                        /* last block flag */
    uint8_t buf[BLAKE2B_BLOCK_LENGTH]; /* input buffer */
    size_t buflen;                     /* size of remaining bytes in buf */
    size_t outlen;                     /* digest size */
};

#endif /* BLAKE_H */
