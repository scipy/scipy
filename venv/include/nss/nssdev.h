/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef NSSDEV_H
#define NSSDEV_H

/*
 * nssdev.h
 *
 * High-level methods for interaction with cryptoki devices
 */

#ifndef NSSDEVT_H
#include "nssdevt.h"
#endif /* NSSDEVT_H */

PR_BEGIN_EXTERN_C

/* NSSAlgorithmAndParameters
 *
 * NSSAlgorithmAndParameters_CreateSHA1Digest
 * NSSAlgorithmAndParameters_CreateMD5Digest
 */

NSS_EXTERN NSSAlgorithmAndParameters *
NSSAlgorithmAndParameters_CreateSHA1Digest(
    NSSArena *arenaOpt);

NSS_EXTERN NSSAlgorithmAndParameters *
NSSAlgorithmAndParameters_CreateMD5Digest(
    NSSArena *arenaOpt);

PR_END_EXTERN_C

#endif /* DEV_H */
