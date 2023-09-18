/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef CK_H
#define CK_H

/*
 * ck.h
 *
 * This header file consolidates all header files needed by the source
 * files implementing the NSS Cryptoki Framework.  This makes managing
 * the source files a bit easier.
 */

/* Types */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

#ifndef NSSCKT_H
#include "nssckt.h"
#endif /* NSSCKT_H */

#ifndef NSSCKFT_H
#include "nssckft.h"
#endif /* NSSCKFT_H */

#ifndef NSSCKEPV_H
#include "nssckepv.h"
#endif /* NSSCKEPV_H */

#ifndef NSSCKFWT_H
#include "nssckfwt.h"
#endif /* NSSCKFWT_H */

#ifndef NSSCKMDT_H
#include "nssckmdt.h"
#endif /* NSSCKMDT_H */

#ifndef CKT_H
#include "ckt.h"
#endif /* CKT_H */

#ifndef CKFWTM_H
#include "ckfwtm.h"
#endif /* CKFWTM_H */

/* Prototypes */

#ifndef NSSBASE_H
#include "nssbase.h"
#endif /* NSSBASE_H */

#ifndef NSSCKG_H
#include "nssckg.h"
#endif /* NSSCKG_H */

#ifndef NSSCKFW_H
#include "nssckfw.h"
#endif /* NSSCKFW_H */

#ifndef NSSCKFWC_H
#include "nssckfwc.h"
#endif /* NSSCKFWC_H */

#ifndef CKFW_H
#include "ckfw.h"
#endif /* CKFW_H */

#ifndef CKFWM_H
#include "ckfwm.h"
#endif /* CKFWM_H */

#ifndef CKMD_H
#include "ckmd.h"
#endif /* CKMD_H */

/* NSS-private */

/* nss_ZNEW and the like.  We might want to publish the memory APIs.. */

#ifndef BASE_H
#include "base.h"
#endif /* BASE_H */

#endif /* CK_H */
