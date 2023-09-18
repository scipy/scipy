/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_revocationchecker.h
 *
 * RevocationChecker Object Type Definition
 *
 */

#ifndef _PKIX_REVOCATIONCHECKER_H
#define _PKIX_REVOCATIONCHECKER_H

#include "pkixt.h"
#include "certt.h"

#ifdef __cplusplus
extern "C" {
#endif

/* NOTE: nbio logistic removed. Will be replaced later. */

/*
 * All Flags are prefixed by CERT_REV_M_, where _M_ indicates
 * this is a method dependent flag.
 */

/*
 * Whether or not to use a method for revocation testing.
 * If set to "do not test", then all other flags are ignored.
 */
#define PKIX_REV_M_DO_NOT_TEST_USING_THIS_METHOD     0x00L
#define PKIX_REV_M_TEST_USING_THIS_METHOD            0x01L

/*
 * Whether or not NSS is allowed to attempt to fetch fresh information
 *         from the network.
 * (Although fetching will never happen if fresh information for the
 *           method is already locally available.)
 */
#define PKIX_REV_M_ALLOW_NETWORK_FETCHING            0x00L
#define PKIX_REV_M_FORBID_NETWORK_FETCHING           0x02L

/*
 * Example for an implicit default source:
 *         The globally configured default OCSP responder.
 * IGNORE means:
 *        ignore the implicit default source, whether it's configured or not.
 * ALLOW means:
 *       if an implicit default source is configured, 
 *          then it overrides any available or missing source in the cert.
 *       if no implicit default source is configured,
 *          then we continue to use what's available (or not available) 
 *          in the certs.
 */ 
#define PKIX_REV_M_ALLOW_IMPLICIT_DEFAULT_SOURCE     0x00L
#define PKIX_REV_M_IGNORE_IMPLICIT_DEFAULT_SOURCE    0x04L /* OCSP only */

/*
 * Defines the behavior if no fresh information is available,
 *   fetching from the network is allowed, but the source of revocation
 *   information is unknown (even after considering implicit sources,
 *   if allowed by other flags).
 * SKIPT_TEST means:
 *          We ignore that no fresh information is available and 
 *          skip this test.
 * REQUIRE_INFO means:
 *          We still require that fresh information is available.
 *          Other flags define what happens on missing fresh info.
 */

#define PKIX_REV_M_SKIP_TEST_ON_MISSING_SOURCE       0x00L
#define PKIX_REV_M_REQUIRE_INFO_ON_MISSING_SOURCE    0x08L

/*
 * Defines the behavior if we are unable to obtain fresh information.
 * INGORE means:
 *      Return "cert status unknown"
 * FAIL means:
 *      Return "cert revoked".
 */

#define PKIX_REV_M_IGNORE_MISSING_FRESH_INFO         0x00L
#define PKIX_REV_M_FAIL_ON_MISSING_FRESH_INFO        0x10L

/*
 * What should happen if we were able to find fresh information using
 * this method, and the data indicated the cert is good?
 * STOP_TESTING means:
 *              Our success is sufficient, do not continue testing
 *              other methods.
 * CONTINUE_TESTING means:
 *                  We will continue and test the next allowed
 *                  specified method.
 */

#define PKIX_REV_M_STOP_TESTING_ON_FRESH_INFO        0x00L
#define PKIX_REV_M_CONTINUE_TESTING_ON_FRESH_INFO    0x20L

/*
 * All Flags are prefixed by PKIX_REV_MI_, where _MI_ indicates
 * this is a method independent flag.
 */

/*
 * This defines the order to checking.
 * EACH_METHOD_SEPARATELY means:
 *      Do all tests related to a particular allowed method
 *      (both local information and network fetching) in a single step.
 *      Only after testing for a particular method is done,
 *      then switching to the next method will happen.
 * ALL_LOCAL_INFORMATION_FIRST means:
 *      Start by testing the information for all allowed methods
 *      which are already locally available. Only after that is done
 *      consider to fetch from the network (as allowed by other flags).
 */
#define PKIX_REV_MI_TEST_EACH_METHOD_SEPARATELY       0x00L
#define PKIX_REV_MI_TEST_ALL_LOCAL_INFORMATION_FIRST  0x01L

/*
 * Use this flag to specify that it's necessary that fresh information
 * is available for at least one of the allowed methods, but it's
 * irrelevant which of the mechanisms succeeded.
 * NO_OVERALL_INFO_REQUIREMENT means:
 *     We strictly follow the requirements for each individual method.
 * REQUIRE_SOME_FRESH_INFO_AVAILABLE means:
 *     After the individual tests have been executed, we must have
 *     been able to find fresh information using at least one method.
 *     If we were unable to find fresh info, it's a failure.
 */
#define PKIX_REV_MI_NO_OVERALL_INFO_REQUIREMENT       0x00L
#define PKIX_REV_MI_REQUIRE_SOME_FRESH_INFO_AVAILABLE 0x02L

/* Defines check time for the cert, revocation methods lists and
 * flags for leaf and chain certs revocation tests. */
struct PKIX_RevocationCheckerStruct {
    PKIX_List *leafMethodList;
    PKIX_List *chainMethodList;
    PKIX_UInt32 leafMethodListFlags;
    PKIX_UInt32 chainMethodListFlags;
};

/* see source file for function documentation */

PKIX_Error *pkix_RevocationChecker_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_REVOCATIONCHECKER_H */
