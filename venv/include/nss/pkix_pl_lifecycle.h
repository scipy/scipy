/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_lifecycle.h
 *
 * Lifecycle Definitions
 *
 */

#ifndef _PKIX_PL_LIFECYCLE_H
#define _PKIX_PL_LIFECYCLE_H

#include "pkix_pl_common.h"
#include "pkix_pl_oid.h"
#include "pkix_pl_aiamgr.h"
#include "pkix_pl_bigint.h"
#include "pkix_pl_bytearray.h"
#include "pkix_pl_hashtable.h"
#include "pkix_pl_mutex.h"
#include "pkix_pl_rwlock.h"
#include "pkix_pl_monitorlock.h"
#include "pkix_pl_string.h"
#include "pkix_pl_cert.h"
#include "pkix_pl_x500name.h"
#include "pkix_pl_generalname.h"
#include "pkix_pl_publickey.h"
#include "pkix_pl_date.h"
#include "pkix_pl_basicconstraints.h"
#include "pkix_pl_certpolicyinfo.h"
#include "pkix_pl_certpolicymap.h"
#include "pkix_pl_certpolicyqualifier.h"
#include "pkix_pl_crlentry.h"
#include "pkix_pl_crl.h"
#include "pkix_pl_colcertstore.h"
#ifndef NSS_PKIX_NO_LDAP
#include "pkix_pl_ldapcertstore.h"
#include "pkix_pl_ldapdefaultclient.h"
#include "pkix_pl_ldaprequest.h"
#include "pkix_pl_ldapresponse.h"
#endif /* !NSS_PKIX_NO_LDAP */
#include "pkix_pl_socket.h"
#include "pkix_pl_infoaccess.h"
#include "pkix_store.h"
#include "pkix_error.h"
#include "pkix_logger.h"
#include "pkix_list.h"
#include "pkix_trustanchor.h"
#include "pkix_procparams.h"
#include "pkix_valparams.h"
#include "pkix_valresult.h"
#include "pkix_verifynode.h"
#include "pkix_resourcelimits.h"
#include "pkix_certchainchecker.h"
#include "pkix_revocationchecker.h"
#include "pkix_certselector.h"
#include "pkix_comcertselparams.h"
#include "pkix_crlselector.h"
#include "pkix_comcrlselparams.h"
#include "pkix_targetcertchecker.h"
#include "pkix_basicconstraintschecker.h"
#include "pkix_policynode.h"
#include "pkix_policychecker.h"
#include "pkix_crlchecker.h"
#include "pkix_signaturechecker.h"
#include "pkix_buildresult.h"
#include "pkix_build.h"
#include "pkix_pl_nameconstraints.h"
#include "pkix_nameconstraintschecker.h"
#include "pkix_ocspchecker.h"
#include "pkix_pl_ocspcertid.h"
#include "pkix_pl_ocsprequest.h"
#include "pkix_pl_ocspresponse.h"
#include "pkix_pl_httpdefaultclient.h"
#include "pkix_pl_httpcertstore.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_InitializeParamsStruct {
        PKIX_List *loggers;
        PKIX_UInt32 majorVersion;
        PKIX_UInt32 minorVersion;
};

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_LIFECYCLE_H */
