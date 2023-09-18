/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_comcertselparams.h
 *
 * ComCertSelParams Object Type Definition
 *
 */

#ifndef _PKIX_COMCERTSELPARAMS_H
#define _PKIX_COMCERTSELPARAMS_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * pathToNamesConstraint is Name Constraints generated based on the
 * pathToNames. We save a cached copy to save regeneration for each
 * check. SubjAltNames also has its cache, since SubjAltNames are
 * verified by checker, its cache copy is stored in checkerstate.
 */
struct PKIX_ComCertSelParamsStruct {
        PKIX_Int32 version;
        PKIX_Int32 minPathLength;
        PKIX_Boolean matchAllSubjAltNames;
        PKIX_PL_X500Name *subject;
        PKIX_List *policies; /* List of PKIX_PL_OID */
        PKIX_PL_Cert *cert;
        PKIX_PL_CertNameConstraints *nameConstraints;
        PKIX_List *pathToNames; /* List of PKIX_PL_GeneralNames */
        PKIX_List *subjAltNames; /* List of PKIX_PL_GeneralNames */
        PKIX_List *extKeyUsage; /* List of PKIX_PL_OID */
        PKIX_UInt32 keyUsage;
        PKIX_PL_Date *date;
        PKIX_PL_Date *certValid;
        PKIX_PL_X500Name *issuer;
        PKIX_PL_BigInt *serialNumber;
        PKIX_PL_ByteArray *authKeyId;
        PKIX_PL_ByteArray *subjKeyId;
        PKIX_PL_PublicKey *subjPubKey;
        PKIX_PL_OID *subjPKAlgId;
        PKIX_Boolean leafCertFlag;
};

/* see source file for function documentation */

PKIX_Error *pkix_ComCertSelParams_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_COMCERTSELPARAMS_H */
