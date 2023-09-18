/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _LDAP_H_
#define _LDAP_H_

#include "certt.h"
#include "pkixt.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const SEC_ASN1Template PKIX_PL_LDAPCrossCertPairTemplate[];
SEC_ASN1_CHOOSER_DECLARE(PKIX_PL_LDAPCrossCertPairTemplate)
extern const SEC_ASN1Template PKIX_PL_LDAPMessageTemplate[];
SEC_ASN1_CHOOSER_DECLARE(PKIX_PL_LDAPMessageTemplate)
extern const SEC_ASN1Template LDAPFilterTemplate[];
SEC_ASN1_CHOOSER_DECLARE(LDAPFilterTemplate)

/* ********************************************************************** */

#define SEC_ASN1_LDAP_STRING SEC_ASN1_OCTET_STRING

#define LDAPATTR_CACERT         (1<<0)
#define LDAPATTR_USERCERT       (1<<1)
#define LDAPATTR_CROSSPAIRCERT  (1<<2)
#define LDAPATTR_CERTREVLIST    (1<<3)
#define LDAPATTR_AUTHREVLIST    (1<<4)
#define MAX_LDAPATTRS                   5
typedef PKIX_UInt32 LdapAttrMask;

typedef enum {
        SIMPLE_AUTH                     = 0,
        KRBV42LDAP_AUTH                 = 1,
        KRBV42DSA_AUTH                  = 2
} AuthType;

typedef enum {
        BASE_OBJECT                     = 0,
        SINGLE_LEVEL                    = 1,
        WHOLE_SUBTREE                   = 2
} ScopeType;

typedef enum {
        NEVER_DEREF                     = 0,
        DEREF_IN_SEARCHING              = 1,
        DEREF_FINDING_BASEOBJ           = 2,
        ALWAYS_DEREF                    = 3
} DerefType;

typedef enum {
        LDAP_INITIALSUBSTRING_TYPE      = 0,
        LDAP_ANYSUBSTRING_TYPE          = 1,
        LDAP_FINALSUBSTRING_TYPE        = 2
} LDAPSubstringFilterType;

typedef enum {
        LDAP_ANDFILTER_TYPE             = 0,
        LDAP_ORFILTER_TYPE              = 1,
        LDAP_NOTFILTER_TYPE             = 2,
        LDAP_EQUALFILTER_TYPE           = 3,
        LDAP_SUBSTRINGFILTER_TYPE       = 4,
        LDAP_GREATEROREQUALFILTER_TYPE  = 5,
        LDAP_LESSOREQUALFILTER_TYPE     = 6,
        LDAP_PRESENTFILTER_TYPE         = 7,
        LDAP_APPROXMATCHFILTER_TYPE     = 8
} LDAPFilterType;

typedef enum {
        LDAP_BIND_TYPE                  = 0,
        LDAP_BINDRESPONSE_TYPE          = 1,
        LDAP_UNBIND_TYPE                = 2,
        LDAP_SEARCH_TYPE                = 3,
        LDAP_SEARCHRESPONSEENTRY_TYPE   = 4,
        LDAP_SEARCHRESPONSERESULT_TYPE  = 5,
        LDAP_ABANDONREQUEST_TYPE        = 16
} LDAPMessageType;

typedef enum {
        SUCCESS                         = 0,
        OPERATIONSERROR                 = 1,
        PROTOCOLERROR                   = 2,
        TIMELIMITEXCEEDED               = 3,
        SIZELIMITEXCEEDED               = 4,
        COMPAREFALSE                    = 5,
        COMPARETRUE                     = 6,
        AUTHMETHODNOTSUPPORTED          = 7,
        STRONGAUTHREQUIRED              = 8,
        NOSUCHATTRIBUTE                 = 16,
        UNDEFINEDATTRIBUTETYPE          = 17,
        INAPPROPRIATEMATCHING           = 18,
        CONSTRAINTVIOLATION             = 19,
        ATTRIBUTEORVALUEEXISTS          = 20,
        INVALIDATTRIBUTESYNTAX          = 21,
        NOSUCHOBJECT                    = 32,
        ALIASPROBLEM                    = 33,
        INVALIDDNSYNTAX                 = 34,
        ISLEAF                          = 35,
        ALIASDEREFERENCINGPROBLEM       = 36,
        INAPPROPRIATEAUTHENTICATION     = 48,
        INVALIDCREDENTIALS              = 49,
        INSUFFICIENTACCESSRIGHTS        = 50,
        BUSY                            = 51,
        UNAVAILABLE                     = 52,
        UNWILLINGTOPERFORM              = 53,
        LOOPDETECT                      = 54,
        NAMINGVIOLATION                 = 64,
        OBJECTCLASSVIOLATION            = 65,
        NOTALLOWEDONNONLEAF             = 66,
        NOTALLOWEDONRDN                 = 67,
        ENTRYALREADYEXISTS              = 68,
        OBJECTCLASSMODSPROHIBITED       = 69,
        OTHER                           = 80
} LDAPResultCode;

typedef struct LDAPLocationStruct                LDAPLocation;
typedef struct LDAPCertPairStruct                LDAPCertPair;
typedef struct LDAPSimpleBindStruct              LDAPSimpleBind;
typedef struct LDAPBindAPIStruct                 LDAPBindAPI;
typedef struct LDAPBindStruct                    LDAPBind;
typedef struct LDAPResultStruct                  LDAPBindResponse;
typedef struct LDAPResultStruct                  LDAPResult;
typedef struct LDAPSearchResponseAttrStruct      LDAPSearchResponseAttr;
typedef struct LDAPSearchResponseEntryStruct     LDAPSearchResponseEntry;
typedef struct LDAPResultStruct                  LDAPSearchResponseResult;
typedef struct LDAPUnbindStruct                  LDAPUnbind;
typedef struct LDAPFilterStruct                  LDAPFilter;
typedef struct LDAPAndFilterStruct               LDAPAndFilter;
typedef struct LDAPNotFilterStruct               LDAPNotFilter;
typedef struct LDAPSubstringStruct               LDAPSubstring;
typedef struct LDAPSubstringFilterStruct         LDAPSubstringFilter;
typedef struct LDAPPresentFilterStruct           LDAPPresentFilter;
typedef struct LDAPAttributeValueAssertionStruct LDAPAttributeValueAssertion;
typedef struct LDAPNameComponentStruct           LDAPNameComponent;
typedef struct LDAPRequestParamsStruct           LDAPRequestParams;
typedef struct LDAPSearchStruct                  LDAPSearch;
typedef struct LDAPAbandonRequestStruct          LDAPAbandonRequest;
typedef struct protocolOpStruct                  LDAPProtocolOp;
typedef struct LDAPMessageStruct                 LDAPMessage;
typedef LDAPAndFilter                            LDAPOrFilter;
typedef LDAPAttributeValueAssertion              LDAPEqualFilter;
typedef LDAPAttributeValueAssertion              LDAPGreaterOrEqualFilter;
typedef LDAPAttributeValueAssertion              LDAPLessOrEqualFilter;
typedef LDAPAttributeValueAssertion              LDAPApproxMatchFilter;

struct LDAPLocationStruct {
        PLArenaPool *arena;
        void *serverSite;
        void **filterString;
        void **attrBitString;
};

struct LDAPCertPairStruct {
        SECItem forward;
        SECItem reverse;
};

struct LDAPSimpleBindStruct {
        char *bindName;
        char *authentication;
};

struct LDAPBindAPIStruct {
        AuthType selector;
        union {
                LDAPSimpleBind simple;
        } chooser;
};

struct LDAPBindStruct {
        SECItem version;
        SECItem bindName;
        SECItem authentication;
};

struct LDAPResultStruct {
        SECItem resultCode;
        SECItem matchedDN;
        SECItem errorMessage;
};

struct LDAPSearchResponseAttrStruct {
        SECItem attrType;
        SECItem **val;
};

struct LDAPSearchResponseEntryStruct {
        SECItem objectName;
        LDAPSearchResponseAttr **attributes;
};

struct LDAPUnbindStruct {
        SECItem dummy;
};

struct LDAPAndFilterStruct {
        LDAPFilter **filters;
};

struct LDAPNotFilterStruct {
        LDAPFilter *filter;
};

struct LDAPSubstringStruct {
        LDAPSubstringFilterType selector;
        SECItem item;
};

struct LDAPSubstringFilterStruct {
        SECItem attrType;
        LDAPSubstring *strings;
};

struct LDAPPresentFilterStruct {
        SECItem attrType;
};

struct LDAPAttributeValueAssertionStruct {
        SECItem attrType;
        SECItem attrValue;
};

struct LDAPFilterStruct {
        LDAPFilterType selector;
        union {
                LDAPAndFilter andFilter;
                LDAPOrFilter orFilter;
                LDAPNotFilter notFilter;
                LDAPEqualFilter equalFilter;
                LDAPSubstringFilter substringFilter;
                LDAPGreaterOrEqualFilter greaterOrEqualFilter;
                LDAPLessOrEqualFilter lessOrEqualFilter;
                LDAPPresentFilter presentFilter;
                LDAPApproxMatchFilter approxMatchFilter;
        } filter;
};

struct LDAPNameComponentStruct {
        unsigned char *attrType;
        unsigned char *attrValue;
};

struct LDAPRequestParamsStruct {
        char *baseObject;          /* e.g. "c=US" */
        ScopeType scope;
        DerefType derefAliases;
        PKIX_UInt32 sizeLimit;     /* 0 = no limit */
        PRIntervalTime timeLimit;  /* 0 = no limit */
        LDAPNameComponent **nc; /* e.g. {{"cn","xxx"},{"o","yyy"},NULL} */
        LdapAttrMask attributes;
};

struct LDAPSearchStruct {
        SECItem baseObject;
        SECItem scope;
        SECItem derefAliases;
        SECItem sizeLimit;
        SECItem timeLimit;
        SECItem attrsOnly;
        LDAPFilter filter;
        SECItem **attributes;
};

struct LDAPAbandonRequestStruct {
        SECItem messageID;
};

struct protocolOpStruct {
        LDAPMessageType selector;
        union {
                LDAPBind bindMsg;
                LDAPBindResponse bindResponseMsg;
                LDAPUnbind unbindMsg;
                LDAPSearch searchMsg;
                LDAPSearchResponseEntry searchResponseEntryMsg;
                LDAPSearchResponseResult searchResponseResultMsg;
                LDAPAbandonRequest abandonRequestMsg;
        } op;
};

struct LDAPMessageStruct {
        SECItem messageID;
        LDAPProtocolOp protocolOp;
};

typedef struct PKIX_PL_LdapClientStruct PKIX_PL_LdapClient;

typedef PKIX_Error *
(*PKIX_PL_LdapClient_InitiateFcn)(
        PKIX_PL_LdapClient *client,
        LDAPRequestParams *requestParams,
        void **pNBIO,
        PKIX_List **pResponse,
        void *plContext);

typedef PKIX_Error *
(*PKIX_PL_LdapClient_ResumeFcn)(
        PKIX_PL_LdapClient *client,
        void **pNBIO,
        PKIX_List **pResponse,
        void *plContext);

struct PKIX_PL_LdapClientStruct {
        PKIX_PL_LdapClient_InitiateFcn initiateFcn;
        PKIX_PL_LdapClient_ResumeFcn resumeFcn;
};

#ifdef __cplusplus
}
#endif

#endif
