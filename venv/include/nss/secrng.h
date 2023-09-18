/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _SECRNG_H_
#define _SECRNG_H_
/*
 * secrng.h - public data structures and prototypes for the secure random
 *	      number generator
 */

/******************************************/
/*
** Random number generation. A cryptographically strong random number
** generator.
*/

#include "blapi.h"

/* the number of bytes to read from the system random number generator */
#define SYSTEM_RNG_SEED_COUNT 1024

SEC_BEGIN_PROTOS

/*
** The following functions are provided by the security library
** but are differently implemented for the UNIX, Win, and OS/2
** versions
*/

/*
** Get the "noisiest" information available on the system.
** The amount of data returned depends on the system implementation.
** It will not exceed maxbytes, but may be (much) less.
** Returns number of noise bytes copied into buf, or zero if error.
*/
extern size_t RNG_GetNoise(void *buf, size_t maxbytes);

/*
** RNG_SystemInfoForRNG should be called before any use of SSL. It
** gathers up the system specific information to help seed the
** state of the global random number generator.
*/
extern void RNG_SystemInfoForRNG(void);

/*
** Use the contents (and stat) of a file to help seed the
** global random number generator.
*/
extern void RNG_FileForRNG(const char *filename);

/*
** Get maxbytes bytes of random data from the system random number
** generator.
** Returns the number of bytes copied into buf -- maxbytes if success
** or zero if error.
** Errors:
**   PR_NOT_IMPLEMENTED_ERROR   There is no system RNG on the platform.
**   SEC_ERROR_NEED_RANDOM      The system RNG failed.
*/
extern size_t RNG_SystemRNG(void *buf, size_t maxbytes);

SEC_END_PROTOS

#endif /* _SECRNG_H_ */
