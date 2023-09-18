/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _KEYI_H_
#define _KEYI_H_
#include "secerr.h"

SEC_BEGIN_PROTOS
/* NSS private functions */
/* map an oid to a keytype... actually this function and it's converse
 *  are good candidates for public functions..  */
KeyType seckey_GetKeyType(SECOidTag pubKeyOid);

/* extract the 'encryption' (could be signing) and hash oids from and
 * algorithm, key and parameters (parameters is the parameters field
 * of a algorithm ID structure (SECAlgorithmID)*/
SECStatus sec_DecodeSigAlg(const SECKEYPublicKey *key, SECOidTag sigAlg,
                           const SECItem *param, SECOidTag *encalg, SECOidTag *hashalg);

/* just get the 'encryption' oid from the combined signature oid */
SECOidTag sec_GetEncAlgFromSigAlg(SECOidTag sigAlg);

/* extract the RSA-PSS hash algorithms and salt length from
 * parameters, taking into account of the default implications.
 *
 * (parameters is the parameters field of a algorithm ID structure
 * (SECAlgorithmID)*/
SECStatus sec_DecodeRSAPSSParams(PLArenaPool *arena,
                                 const SECItem *params,
                                 SECOidTag *hashAlg,
                                 SECOidTag *maskHashAlg,
                                 unsigned long *saltLength);

/* convert the encoded RSA-PSS parameters into PKCS #11 mechanism parameters */
SECStatus sec_DecodeRSAPSSParamsToMechanism(PLArenaPool *arena,
                                            const SECItem *params,
                                            CK_RSA_PKCS_PSS_PARAMS *mech);

/* make sure the key length matches the policy for keyType */
SECStatus seckey_EnforceKeySize(KeyType keyType, unsigned keyLength,
                                SECErrorCodes error);
SEC_END_PROTOS

#endif /* _KEYHI_H_ */
