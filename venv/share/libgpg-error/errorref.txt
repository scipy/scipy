# errorref.txt - Description of error codes
# Copyright (C) 2003-2004, 2010, 2013-2016 g10 Code GmbH
#
# This file is part of libgpg-error.
#
# libgpg-error is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
#
# libgpg-error is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, see <https://www.gnu.org/licenses/>.


GPG_ERR_UNKNOWN_PACKET          Unknown packet

    GNUPG:  - Redefined to G10ERR_UNKNOWN_PACKET in gpg.

GPG_ERR_UNKNOWN_VERSION         Unknown version in packet

    Used by GnuPG 2.1 to identify valid OpenPGP packets with an
    unknown version.

GPG_ERR_PUBKEY_ALGO             Invalid public key algorithm

    GNUPG:  - Redefined to G10ERR_PUBKEY_ALGO in gpg.
            - Public key algorithm is not allowed by OpenPGP.
    GCRYPT: - Public key algorithm is not defined or not available.
              Note that this is also the case if the algorithm has
              been disabled.
            - [version < 1.5] Checking of the RSA secret key failed
                              (consistency check).

GPG_ERR_DIGEST_ALGO             Invalid digest algorithm

    GNUPG:  - Digest algorithm is not supported.
            - Redefined to G10ERR_PUBKEY_ALGO in gpg.
            - Digest algorithm is not allowed by OpenPGP.
            - Unsupported algorithm given to "--hash=" option of
              certain Assuan server commands.
            - Signature creation or verification failed due to
              an unsupported hash algorithm.
    GCRYPT: - Digest key algorithm is not defined or not available.
              Note that this is also the case if the algorithm has
              been disabled.
            - Unsupported digest algorithm in a selftest.
            - Invalid digest algorithm used in FIPS mode.  Note that
              in enforced-FIPS mode there is no such error return.
            - Message digested or HMAC computation finished with no
              message algorithm enabled for the hash context.
            - Bad digest algorithm given to public key function.

GPG_ERR_BAD_PUBKEY              Bad public key

    GNUPG:  - Redefined to G10ERR_BAD_PUBKEY in gpg.
            - Missing public or domain parameter in an s-expression.
              If the curve name is mssing GPG_ERR_INV_CURVE may be
              used as well.

GPG_ERR_BAD_SECKEY              Bad secret key

    GNUPG:  - Invalid format of a S-expression encoded private key in
              gpg-agent.
            - Missing secret parameter in an s-expression.
            - A protected or shadowed private key was passed to the
              OpenPGP card application for storing it on the card.
            - A private key passed to the OpenPGP card application does
              not match the requirements of the card or misses required
              parameters.
            - Gpg'agents import key command is not able to convert
              the key to the internal format.
    GCRYPT: - Checking the secret key failed (consistency check).


GPG_ERR_BAD_SIGNATURE           Bad signature

    GNUPG:  - Redefined to G10ERR_BAD_SIGN in gpg.
            - The MDC check of an OpenPGP encrypted message failed.
            - A OpenPGP key signature did not verify.
            - A signature with a key flagged as "never trust" was made.
    GCRYPT: - A public key signature did not verify.

GPG_ERR_NO_PUBKEY               No public key

    GNUPG:  - Redefined to G10ERR_NO_PUBKEY in gpg.
            - A key was requested from an OpenPGP card but the key is
              not stored on the card.
            - The public key could not be retrieved from a corresponding
              certificate on a card (command READKEY in scd).
            - A requested certificate was not found or an unspecified
              error occurred while selecting a X.509 certificate in
              gpgsm.
            - The specified certificate or key was not found.  This
              does not necessary mean that the certifciate is not
              available but the specification method may not be usable
              for the given certificate.  May also happen for
              certificates somewhere in the chain while validaiting a
              certificate chain.
            - The requested encryption certificate was not found.
            - A certificate specified in a CMS message is not
              available and thus the signature could not be verified
              or details of the certificate be shown.
    GPA:    - No key was given for encryption.
            - The selected encryption protocol is not available.

GPG_ERR_CHECKSUM                Checksum error

    GNUPG:  - The checksum of an unprotected OpenPGP secret key packet
              is wrong.
    GCRYPT: - Decryption in AESWRAP mode does not match the expected IV.
    [more to come]

GPG_ERR_BAD_PASSPHRASE          Bad passphrase

    GNUPG: - The entered passphrase does not verify

GPG_ERR_CIPHER_ALGO             Invalid cipher algorithm

GPG_ERR_KEYRING_OPEN            Cannot open keyring

GPG_ERR_INV_PACKET              Invalid packet

GPG_ERR_INV_ARMOR               Invalid armor

GPG_ERR_NO_USER_ID              No user ID

GPG_ERR_NO_SECKEY               No secret key

   NTBTLS: - No private key or pre-shared key available.

GPG_ERR_WRONG_SECKEY            Wrong secret key used

GPG_ERR_BAD_KEY                         Bad session key

    GNUPG: - gpg-agent's command IMPORT_KEY or EXPORT_KEY is used
             without a prior KEYWRAP_KEY command.

    [more to come]


GPG_ERR_COMPR_ALGO              Unknown compression algorithm

GPG_ERR_NO_PRIME                Number is not prime

GPG_ERR_NO_ENCODING_METHOD      Invalid encoding method

GPG_ERR_NO_ENCRYPTION_SCHEME    Invalid encryption scheme

GPG_ERR_NO_SIGNATURE_SCHEME     Invalid signature scheme

GPG_ERR_INV_ATTR                Invalid attribute


GPG_ERR_NO_VALUE                No value

    GNUPG:  - A timestamp value is expect but there is none.
    KSBA:   - A timestamp value is expect but there is none.
            - A certificate is missing a required property.
            - A CMS object is missing a required property.
            - Converting a Distinguised Name to an RFC2253 string failed.


GPG_ERR_NOT_FOUND               Not found

    A search operation did not return a matching value.


GPG_ERR_VALUE_NOT_FOUND                 Value not found

    GNUPG:  - A keyblock or a cert object was requested but not
              found.  This might indicate an internal error here.


GPG_ERR_SYNTAX                  Syntax error

GPG_ERR_BAD_MPI                         Bad MPI value

GPG_ERR_INV_PASSPHRASE          Invalid passphrase

    GNUPG:  - Required constraints of the passphrase are not met.

GPG_ERR_SIG_CLASS               Invalid signature class

GPG_ERR_RESOURCE_LIMIT          Resources exhausted

GPG_ERR_INV_KEYRING             Invalid keyring

GPG_ERR_TRUSTDB                         Trust DB error


GPG_ERR_BAD_CERT                Bad certificate

   NTBTLS: - No subject found in the certifciate.


GPG_ERR_INV_USER_ID             Invalid user ID

    GNUPG:  - Used to indicate a bad specification of a user id.
    [more to come]


GPG_ERR_UNEXPECTED              Unexpected error

GPG_ERR_TIME_CONFLICT           Time conflict

GPG_ERR_KEYSERVER               Keyserver error


GPG_ERR_WRONG_PUBKEY_ALGO       Wrong public key algorithm

    GNUPG: - The algorithm is not expected.  For example a DSA
             algorithm is used where a non-DSA algorithm is expected
             or vice versa.  May indicate an internal error.
    NTBTLS: - Public key type mismatch.  The peer presented a
              different key type than requested.


GPG_ERR_TRIBUTE_TO_D_A          Tribute to D. A.

GPG_ERR_WEAK_KEY                Weak encryption key

GPG_ERR_INV_KEYLEN              Invalid key length

GPG_ERR_INV_ARG                         Invalid argument

    GCRYPT:  - Unsupported length of input data in encrypt or decrypt
               cipher functions.  For example not matching the block
               lengths of the algorithm.
             - Incompatible args given; e.g. two or none if exactly one
               is required.
    [more to come]


GPG_ERR_BAD_URI                         Syntax error in URI

GPG_ERR_INV_URI                         Invalid URI

GPG_ERR_NETWORK                         Network error

GPG_ERR_UNKNOWN_HOST            Unknown host

        Used instead of the non-portable EHOSTNOTFOUND which is
        returned by some systems as a mapping of h_errno's
        HOST_NOT_FOUND


GPG_ERR_SELFTEST_FAILED                 Selftest failed

GPG_ERR_NOT_ENCRYPTED           Data not encrypted

GPG_ERR_NOT_PROCESSED           Data not processed

GPG_ERR_UNUSABLE_PUBKEY                 Unusable public key

GPG_ERR_UNUSABLE_SECKEY                 Unusable secret key

GPG_ERR_INV_VALUE               Invalid value

    NTBTLS: - A DH parameter is out of range
    GnuPG:  - An Assuan server returns a status line with
              unexpected values.

GPG_ERR_BAD_CERT_CHAIN          Bad certificate chain

GPG_ERR_MISSING_CERT            Missing certificate

    NTBTLS: - The server needs to send a certifciate but none has been
              set.  See also GPG_ERR_MISSING_ISSUER_CERT and
              GPG_ERR_MISSING_CLIENT_CERT.



GPG_ERR_NO_DATA                         No data

GPG_ERR_BUG                     Bug

GPG_ERR_NOT_SUPPORTED           Not supported

        Used if a feature is currently not supported but may be
        enabled for example using a program option.  Commonly used if
        a feature has been disabled by an administrator.  See also
        GPG_ERR_NOT_ENABLED.  Sometimes also used for features which
        are not yet supported.


GPG_ERR_INV_OP                  Invalid operation code

GPG_ERR_TIMEOUT                 Timeout

    Some function or network access timed out.

GPG_ERR_INTERNAL                Internal error

GPG_ERR_EOF_GCRYPT              EOF (gcrypt)

GPG_ERR_INV_OBJ                 Invalid object

GPG_ERR_TOO_SHORT               Provided object is too short

GPG_ERR_TOO_LARGE               Provided object is too large

GPG_ERR_NO_OBJ                  Missing item in object

GPG_ERR_NOT_IMPLEMENTED                 Not implemented

    NTBTLS: - The requested feature is not implemented.

GPG_ERR_CONFLICT                Conflicting use

    NTBTLS: - Function has already been called and may not be called
              again at this protocol state.
    GNUPG:  - Returned by g13 when creating a new container on a device
              which seems to be in use.


GPG_ERR_INV_CIPHER_MODE                 Invalid cipher mode

GPG_ERR_INV_FLAG                Invalid flag

   GPGME: Used to indicate an invalid combination of flags.


GPG_ERR_INV_HANDLE              Invalid handle

GPG_ERR_TRUNCATED               Result truncated

GPG_ERR_INCOMPLETE_LINE                 Incomplete line

GPG_ERR_INV_RESPONSE            Invalid response

GPG_ERR_NO_AGENT                No agent running

GPG_ERR_AGENT                   agent error

GPG_ERR_INV_DATA                Invalid data

    GNUPG:  - Used in app-openpgp.c for a badly formatted request.
    GCRYPT: - No passphrase given for gcry_kdf_derive.
            - An opaque MPI is given to a public key function but not
              expected.

GPG_ERR_ASSUAN_SERVER_FAULT     Unspecific Assuan server fault

GPG_ERR_ASSUAN                  General Assuan error

    GNUPG: - Used by Assuan command handler if they fail to do basic
             things like an es_fdopen or es_fopencookie.


GPG_ERR_INV_SESSION_KEY                 Invalid session key

GPG_ERR_INV_SEXP                Invalid S-expression

GPG_ERR_UNSUPPORTED_ALGORITHM   Unsupported algorithm

GPG_ERR_NO_PIN_ENTRY            No pinentry

GPG_ERR_PIN_ENTRY               pinentry error

GPG_ERR_BAD_PIN                         Bad PIN

GPG_ERR_INV_NAME                Invalid name

    GNUPG:  - Formerly used in GPGSM to indicate an error in
              the specification of a user id.  Later replaced by
              GPG_ERR_INV_USER_ID.
            - In G13 to indicate a bad file name (e.g. one with
              an embedded Nul byte when given as escaped string.
            - In SCDAEMON for an unknown attribute name.

    Also used for URLs which have non-acceptable characters for the
    specific application.

    [more to come]

GPG_ERR_BAD_DATA                Bad data

GPG_ERR_INV_PARAMETER           Invalid parameter

    GNUPG:  - Returned if gpg-agent sends a new generated key with
              unknown parameter names.
            - Invalid parameter in the parameter file for key
              generation by gpgsm.

GPG_ERR_WRONG_CARD              Wrong card

GPG_ERR_NO_DIRMNGR              No dirmngr

GPG_ERR_DIRMNGR                         dirmngr error

GPG_ERR_CERT_REVOKED            Certificate revoked

GPG_ERR_NO_CRL_KNOWN            No CRL known

GPG_ERR_CRL_TOO_OLD             CRL too old

GPG_ERR_LINE_TOO_LONG           Line too long

GPG_ERR_NOT_TRUSTED             Not trusted

GPG_ERR_CANCELED                Operation cancelled

GPG_ERR_BAD_CA_CERT             Bad CA certificate

GPG_ERR_CERT_EXPIRED            Certificate expired

GPG_ERR_CERT_TOO_YOUNG          Certificate too young

GPG_ERR_UNSUPPORTED_CERT        Unsupported certificate

GPG_ERR_UNKNOWN_SEXP            Unknown S-expression

GPG_ERR_UNSUPPORTED_PROTECTION  Unsupported protection

GPG_ERR_CORRUPTED_PROTECTION    Corrupted protection

GPG_ERR_AMBIGUOUS_NAME          Ambiguous name

GPG_ERR_CARD                    Card error

GPG_ERR_CARD_RESET              Card reset required

GPG_ERR_CARD_REMOVED            Card removed

GPG_ERR_INV_CARD                Invalid card

GPG_ERR_CARD_NOT_PRESENT        Card not present

GPG_ERR_NO_PKCS15_APP           No PKCS15 application

GPG_ERR_NOT_CONFIRMED           Not confirmed

GPG_ERR_CONFIGURATION           Configuration error

GPG_ERR_NO_POLICY_MATCH                 No policy match

GPG_ERR_INV_INDEX               Invalid index

GPG_ERR_INV_ID                  Invalid ID

GPG_ERR_NO_SCDAEMON             No SmartCard daemon

GPG_ERR_SCDAEMON                SmartCard daemon error

GPG_ERR_UNSUPPORTED_PROTOCOL    Unsupported protocol

    GPG:        - An unsupported keyserver protocol.
    GPG_AGENT:  - Invalid shadow_info protocol (not "t1-v1")
    LIBKSBA:    - Unknown OID of the OCSP response bytes
    GPGME:      - GPGME_PROTOCOL_xxx not supported.
    NTBTLS:     - Handshake protocol version not supported.

GPG_ERR_BAD_PIN_METHOD          Bad PIN method

GPG_ERR_CARD_NOT_INITIALIZED    Card not initialized

    SCDAEMON: - A card function is called but the card has not yet
                been initialized.  This may be due to a conflict with
                another card using connection or due to a bug.

GPG_ERR_UNSUPPORTED_OPERATION   Unsupported operation

GPG_ERR_WRONG_KEY_USAGE                 Wrong key usage

    GNUPG: - Key usage not possible with selected algorithm.

GPG_ERR_NOTHING_FOUND           Nothing found

  Indicates that the operation was not possible because nothing has
  been found.  For example an update request for non existent data.

GPG_ERR_WRONG_BLOB_TYPE                 Wrong blob type

    GNUPG: - The keyboxd returns an unexpected blob
             (e.g. OpenPGP was requested but X.509 returned).

GPG_ERR_MISSING_VALUE           Missing value

    GNUPG: - Not enough parameters for a secret key send to gpg-agent.

    GCRYPT: - A required parameter has not been given.


GPG_ERR_HARDWARE                Hardware problem

GPG_ERR_PIN_BLOCKED             PIN blocked

GPG_ERR_USE_CONDITIONS          Conditions of use not satisfied

    GNUPG: - The PIN given to a smartcard is too short or has
             unacceptable characters so that the smartcard does
             not even try to verify it.
           - The smartcard can't do an operation because some
             intermediate command send to a card is missing or the
             card can't use the provided data due to an unsupported
             algorithm.

GPG_ERR_PIN_NOT_SYNCED          PINs are not synced

GPG_ERR_INV_CRL                         Invalid CRL

GPG_ERR_BAD_BER                         BER error

GPG_ERR_INV_BER                         Invalid BER

GPG_ERR_ELEMENT_NOT_FOUND       Element not found

GPG_ERR_IDENTIFIER_NOT_FOUND    Identifier not found

GPG_ERR_INV_TAG                         Invalid tag

GPG_ERR_INV_LENGTH              Invalid length

    GCRYPT: - Bad block length for certain cipher algorithms and
              modes.
            - Bad length of input data; e.g. not a multiple of the
              block length.
            - A length does not match the size of the digest
              algorithm.
            - Length of signature or public key is not as expected
              (e.g. in EdDSA).
    [more to come]
    GNUPG:  - Invalid hash length for a pubkey
    [more to come]

GPG_ERR_INV_KEYINFO             Invalid key info

    KSBA: - Returned if the ASN.1 Keyinfo structure is not valid

GPG_ERR_UNEXPECTED_TAG          Unexpected tag

GPG_ERR_NOT_DER_ENCODED                 Not DER encoded

GPG_ERR_NO_CMS_OBJ              No CMS object

GPG_ERR_INV_CMS_OBJ             Invalid CMS object

GPG_ERR_UNKNOWN_CMS_OBJ                 Unknown CMS object

GPG_ERR_UNSUPPORTED_CMS_OBJ     Unsupported CMS object

GPG_ERR_UNSUPPORTED_ENCODING    Unsupported encoding

    GNUPG: - Returned by Dirmngr if a keyserver returns a HTML document.


GPG_ERR_UNSUPPORTED_CMS_VERSION         Unsupported CMS version


GPG_ERR_UNKNOWN_ALGORITHM       Unknown algorithm

    GCRYPT:  gcry_kdf_proc for an unknown kdf algorithm

GPG_ERR_INV_ENGINE              Invalid crypto engine

    GPGME: Several uses use cases.  For example:
           - Unexpected format of a status line.

GPG_ERR_PUBKEY_NOT_TRUSTED      Public key not trusted
GPG_ERR_DECRYPT_FAILED          Decryption failed
GPG_ERR_KEY_EXPIRED             Key expired
GPG_ERR_SIG_EXPIRED             Signature expired
GPG_ERR_ENCODING_PROBLEM        Encoding problem

GPG_ERR_INV_STATE               Invalid state

    The state (of a protocol) is not possible or not defined at all.

    NTBTLS: - Data received in an unexpected state.
            - A function is called while not being in the right state.

GPG_ERR_DUP_VALUE               Duplicated value

GPG_ERR_MISSING_ACTION          Missing action

    GNUPG: - In G13 the server command "MOUNT" is used without prior
             use of the command "OPEN".

    others: - The libassuan ce-server test program uses this to
              indicate that the client did not connect to the server
              as requested.

GPG_ERR_MODULE_NOT_FOUND        ASN.1 module not found

GPG_ERR_INV_OID_STRING          Invalid OID string

GPG_ERR_INV_TIME                Invalid time

GPG_ERR_INV_CRL_OBJ             Invalid CRL object

GPG_ERR_UNSUPPORTED_CRL_VERSION         Unsupported CRL version


GPG_ERR_INV_CERT_OBJ            Invalid certificate object

    GPGME: - A bad certificate (gpgme_key_t) has been passed to a
             function.  For example it might be incomplete due to a
             missing fingerprint.
    GNUPG: - A certificate has a length of zero.


GPG_ERR_UNKNOWN_NAME            Unknown name

        Used by GPG to indicate an unknown ECC curve name (may also
        indicate missing ECC support).  It is also used to indicate an
        unsuported parameter name in functions which take a name and
        value to update state.  Note that GPG_ERR_UNKNOWN_CURVE is
        used instead by newer code.

GPG_ERR_LOCALE_PROBLEM          A locale function failed

GPG_ERR_NOT_LOCKED              Not locked

GPG_ERR_PROTOCOL_VIOLATION      Protocol violation

    GNUPG: - Used for invalid HTTP responses.


GPG_ERR_INV_MAC                         Invalid MAC

        The length, algo, or other properties of a MAC are not met.
        See also GPG_ERR_BAD_MAC.


GPG_ERR_INV_REQUEST             Invalid request

GPG_ERR_UNKNOWN_EXTN            Unknown extension

GPG_ERR_UNKNOWN_CRIT_EXTN       Unknown critical extension

GPG_ERR_LOCKED                  Locked

GPG_ERR_UNKNOWN_OPTION          Unknown option

GPG_ERR_UNKNOWN_COMMAND         Unknown command

GPG_ERR_NOT_OPERATIONAL         Not operational

GPG_ERR_NO_PASSPHRASE           No passphrase given

GPG_ERR_NO_PIN                  No PIN given

GPG_ERR_NOT_ENABLED             Not enabled

        Similar to GPG_ERR_NOT_SUPPORTED.  In general this error is
        used for disabled features which can be expected to be enabled
        by the user.


GPG_ERR_NO_ENGINE               No crypto engine

GPG_ERR_MISSING_KEY             Missing key

    GNUPG:  - gpg-agent returns this error on import or export if a key
              wrapping transport key has not been specified.
            - It is used when the name "Key" is not found while looking
              up name value pairs of the extended private key format

    GCRYPT: - A key has not been set when calling a symmetric
              encryption function.

GPG_ERR_TOO_MANY                Too many objects

    GPG: - Dirmngr KS_GET called with too many pattern so that the
           maximum Assuan line length would overflow.
         - gpgsm's command export --secret called with too man keys.
    GPGME: - To many patterns in gpgme-tools's KEYLIST command.

GPG_ERR_LIMIT_REACHED           Limit reached

        A programmed limit has been reached.

        GnuPG: gpgtar: Extract directory can't be created because too
        many of directories with a similar name are already existing.

GPG_ERR_NOT_INITIALIZED                 Not initialized

    An operation can't be performed because something has not been
    initialized.  This might be a missing initialization of an entire
    subsystems or a prerequisite for using a function is not
    fulfilled.

GPG_ERR_MISSING_ISSUER_CERT     Missing issuer certificate

GPG_ERR_NO_KEYSERVER            No keyserver available

        No keyserver configured or no keyserver available due to
        missing support for the requested protocol.  Found in Dirmngr.

GPG_ERR_INV_CURVE               Invalid elliptic curve

        The curve parameter is missing or the curve is invalid; for
        example it is not possible to get affine coordinates for the
        public key.

GPG_ERR_UNKNOWN_CURVE           Unknown elliptic curve

        The curve is not known or not supported by the protocol.


GPG_ERR_DUP_KEY                 Duplicated key

        A duplicated key was detected.  For example a unique key in a
        database occurred more than once.  Also used if in a protocol
        an expected key was returned more than once.

GPG_ERR_AMBIGUOUS               Ambiguous search

        A search etc returned an ambigious result.  This usually means
        that the search string was not specific enough.

GPG_ERR_NO_CRYPT_CTX            No crypto context

        A crypto context was expected but not given.  Commonly used by
        Libgcrypt.

GPG_ERR_WRONG_CRYPT_CTX                 Wrong crypto context

        The given crypto context does not match the requirements.  For
        example in Libgcrypt a crypto context has private data
        pertaining to certain algorithms.  This error is for example
        returned if a crypto context initialized for a different
        algorithm is used.

GPG_ERR_BAD_CRYPT_CTX           Bad crypto context

        The is a problem with the crypto context.  For example it has
        not been properly initialized.

GPG_ERR_CRYPT_CTX_CONFLICT      Conflict in the crypto context

        Conflicting use of a crypto context.  For example if a context
        is used with objects that don't match the state of the
        context.

GPG_ERR_BROKEN_PUBKEY           Broken public key

        The public key was mathematically not correctly generated.
        (It would have been nicer if we would have used BAD_PUBKEY for
        this, but that error code is in long time use to describe for
        example policy and encoding problems with a key.  Using
        INV_PUBKEY would have been better for these purposes)

GPG_ERR_BROKEN_SECKEY           Broken secret key

        The secret key was mathematically not correctly generated.

GPG_ERR_MAC_ALGO

    GCRYPT: - MAC key algorithm is not defined or not available.


GPG_ERR_FULLY_CANCELED        Operation fully cancelled

GPG_ERR_UNFINISHED            Operation not yet finished

GPG_ERR_BUFFER_TOO_SHORT      Buffer too short

GPG_ERR_SEXP_INV_LEN_SPEC     Invalid length specifier in S-expression

GPG_ERR_SEXP_STRING_TOO_LONG  String too long in S-expression

GPG_ERR_SEXP_UNMATCHED_PAREN  Unmatched parentheses in S-expression

GPG_ERR_SEXP_NOT_CANONICAL    S-expression not canonical

GPG_ERR_SEXP_BAD_CHARACTER    Bad character in S-expression

GPG_ERR_SEXP_BAD_QUOTATION    Bad quotation in S-expression

GPG_ERR_SEXP_ZERO_PREFIX      Zero prefix in S-expression

GPG_ERR_SEXP_NESTED_DH        Nested display hints in S-expression

GPG_ERR_SEXP_UNMATCHED_DH     Unmatched display hints

GPG_ERR_SEXP_UNEXPECTED_PUNC  Unexpected reserved punctuation in S-expression

GPG_ERR_SEXP_BAD_HEX_CHAR     Bad hexadecimal character in S-expression

GPG_ERR_SEXP_ODD_HEX_NUMBERS  Odd hexadecimal numbers in S-expression

GPG_ERR_SEXP_BAD_OCT_CHAR     Bad octal character in S-expression

GPG_ERR_SUBKEYS_EXP_REV       All subkeys are expired or revoked

GPG_ERR_DB_CORRUPTED          Database is corrupted

GPG_ERR_SERVER_FAILED         Server indicated a failure

GPG_ERR_NO_NAME               No name

    GNUPG: - No component given in gpgconf runs.
           - A field name is missing in an import/export filter.
           - "Domain not found".
           - "Host not found".
           - Host or service name not found (EAI_NONAME).
           - No or erroneous SRV record.

GPG_ERR_NO_KEY          No key

    Some kind of key was not found.

GPG_ERR_LEGACY_KEY        Legacy key

    Used by GnuPG to identify version 2 and 3 OpenPGP key packets.

GPG_ERR_REQUEST_TOO_SHORT       Request too short

    A received request is too short to continue processing.

GPG_ERR_REQUEST_TOO_LONG        Request too long

    A received request is too long to continue processing.  This may
    be due to an internal limitation, a protocol violation, or due to
    the use of a newer version of a protocol.

GPG_ERR_OBJ_TERM_STATE          Object is in termination state

    For cards this is the ISO status word 0x6285 (file is in
    termination state).

GPG_ERR_NO_CERT_CHAIN           No certificate chain

    NTBTLS: - A CA chain has not been set but is required.

GPG_ERR_CERT_TOO_LARGE          Certificate is too large

    NTBTLS: - A certificate is too large to be used by the protocol.

GPG_ERR_INV_RECORD              Invalid record

    NTBTLS: - An invalid record was received

GPG_ERR_BAD_MAC                         The MAC does not verify

    NTBTLS: - MAC verification of the message failed.

GPG_ERR_UNEXPECTED_MSG         Unexpected message

    GNUPG:  - An unexpected WKS message was received.
    NTBTLS: - Unexpected message received.

GPG_ERR_COMPR_FAILED           Compression or decompression failed

    NTBTLS: - As the description says.

GPG_ERR_WOULD_WRAP             A counter would wrap

    NTBTLS: - Too many messages exchanged
    Other:  - A counter would wrap.

GPG_ERR_FATAL_ALERT            Fatal alert message received

    NTBTLS: - Fatal alert message received from the peer.

GPG_ERR_NO_CIPHER              No cipher algorithm

    NTBTLS: - Server and client have no algo in common

GPG_ERR_MISSING_CLIENT_CERT     Missing client certificate

    NTBTLS: - No certificate received from client.

GPG_ERR_CLOSE_NOTIFY            Close notification received

    NTBTLS: - Alert with a close notification received

GPG_ERR_TICKET_EXPIRED          Ticket expired

    NTBTLS: - Session ticket has expired.

GPG_ERR_BAD_TICKET              Bad ticket

    NTBTLS: - Bad new session ticket message.

GPG_ERR_UNKNOWN_IDENTITY        Unknown identity

    NTBTLS: - Unknown PSK identify received

GPG_ERR_BAD_HS_CERT             Bad certificate message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_CERT_REQ         Bad certificate request message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_CERT_VER         Bad certificate verify message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_CHANGE_CIPHER    Bad change cipher message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_CLIENT_HELLO     Bad client hello message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_SERVER_HELLO     Bad server hello message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_SERVER_HELLO_DONE  Bad server hello done message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_FINISHED         Bad finished message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_SERVER_KEX       Bad server key exchange message in handshake

    NTBTLS: - As the description says.

GPG_ERR_BAD_HS_CLIENT_KEX       Bad client key exchange message in handshake

    NTBTLS: - As the description says.


GPG_ERR_BOGUS_STRING            Bogus string

    Used if a protocol sends length prefixed strings which contain a
    Nul byte and further processing would discard the rest of the
    string.  May also be used if a string contains unexpected and
    possible dangerous characters (e.g. control characters in a domain
    name).

GPG_ERR_FORBIDDEN               Forbidden

    The use of a features is not allowed due to insufficient rights.
    Use by gpg-agent as an error codes for restricted commands.

GPG_ERR_KEY_DISABLED            Key disabled

    GNUPG: - The key has been disabled by the user.

GPG_ERR_KEY_ON_CARD             Not possible with a card based key

    GNUPG: - The gpg-agent returns this if a DELETE_KEY commands is
             used for a smartcard based key.

GPG_ERR_INV_LOCK_OBJ            Invalid lock object

    GPGRT: - The provided lock object is not valid.  This indicates an
             internal problem in libgpg-error or more likely a
             programming error.

GPG_ERR_TRUE                    True

    Used to return the boolean value True.  Note that GPG_ERR_NO_ERROR
    (with the value 0) is also often used to indicate the value true.

GPG_ERR_FALSE                   False

    Used to return the boolean value False.


GPG_ERR_ASS_GENERAL            General IPC error

GPG_ERR_ASS_ACCEPT_FAILED      IPC accept call failed

GPG_ERR_ASS_CONNECT_FAILED     IPC connect call failed

GPG_ERR_ASS_INV_RESPONSE       Invalid IPC response

GPG_ERR_ASS_INV_VALUE          Invalid value passed to IPC

GPG_ERR_ASS_INCOMPLETE_LINE    Incomplete line passed to IPC

GPG_ERR_ASS_LINE_TOO_LONG      Line passed to IPC too long

GPG_ERR_ASS_NESTED_COMMANDS    Nested IPC commands

GPG_ERR_ASS_NO_DATA_CB         No data callback in IPC

GPG_ERR_ASS_NO_INQUIRE_CB      No inquire callback in IPC

GPG_ERR_ASS_NOT_A_SERVER       Not an IPC server

GPG_ERR_ASS_NOT_A_CLIENT       Not an IPC client

GPG_ERR_ASS_SERVER_START       Problem starting IPC server

GPG_ERR_ASS_READ_ERROR         IPC read error

GPG_ERR_ASS_WRITE_ERROR         IPC write error

GPG_ERR_ASS_TOO_MUCH_DATA      Too much data for IPC layer

GPG_ERR_ASS_UNEXPECTED_CMD     Unexpected IPC command

GPG_ERR_ASS_UNKNOWN_CMD         Unknown IPC command

GPG_ERR_ASS_SYNTAX             IPC syntax error

GPG_ERR_ASS_CANCELED           IPC call has been cancelled

GPG_ERR_ASS_NO_INPUT           No input source for IPC

GPG_ERR_ASS_NO_OUTPUT          No output source for IPC

GPG_ERR_ASS_PARAMETER          IPC parameter error

GPG_ERR_ASS_UNKNOWN_INQUIRE    Unknown IPC inquire

GPG_ERR_ENGINE_TOO_OLD  Crypto engine too old

GPG_ERR_WINDOW_TOO_SMALL        Screen or window too small

    Pinentry: - The size of the screen is too small.

GPG_ERR_WINDOW_TOO_LARGE        Screen or window too large

GPG_ERR_MISSING_ENVVAR  Required environment variable not set

    Pinentry: - The size of the screen can't be determined.

GPG_ERR_USER_ID_EXISTS  User ID already exists

    GNUPG: - Existing user ID in --quick-gen-key.

GPG_ERR_NAME_EXISTS     Name already exists

GPG_ERR_DUP_NAME                Duplicated name

GPG_ERR_TOO_YOUNG               Objects is too young

    For example used if a file is younger than expected.

GPG_ERR_TOO_OLD                         Objects is too old

    Used if an object is too old to be used.  This is a more generic
    code than GPG_ERR_ENGINE_TOO_OLD or GPG_ERR_CRL_TOO_OLD.

GPG_ERR_UNKNOWN_FLAG            Unknown flag

    The flag is not known.

    GNUPG: - The flag part of the string given to the
             option --default-new-key-algo value is not known.

GPG_ERR_INV_ORDER               Invalid execution order

    GNUPG: - In Dirmngr used for the libdns error code DNS_EORDER.

GPG_ERR_ALREADY_FETCHED         Already fetched

    GNUPG: - In Dirmngr used for the libdns error code DNS_EFETCHED.

GPG_ERR_TRY_LATER               Try again later

    This indicates that a server asked to try again later; thus it is
    different from EAGAIN which is used by the local system.  This
    code is for example used instead of h_error's TRY_AGAIN.

GPG_ERR_WRONG_NAME              Wrong name

    NTBTLS: - Hostname does not match the certificate

GPG_ERR_NO_AUTH                 Not authenticated

    GnuPG: - A smartcard requires authentication

GPG_ERR_BAD_AUTH                Bad authentication

    GnuPG: - A smartcard could not be authenticated.  For example
             a wrong authentication key was used with a PIV card.

GPG_ERR_NO_KEYBOXD              No Keyboxd running

    GnuPG: - The keyboxd component is not running

GPG_ERR_KEYBOXD                 Keyboxd error

    GnuPG: - Malfunction in the keyboxd

GPG_ERR_NO_SERVICE              Service is not running

    A component is not running.  Tnis is a generic version of
    GPG_ERR_NO_AGENT et al.

GPG_ERR_SERVICE                 Service error

    An error occured in a service component.  This is a generic
    version of GPG_ERR_AGENT et al.

GPG_ERR_BAD_PUK                 Bad PUK

    A wrong PIN Unblocking Code was used.

GPG_ERR_NO_RESET_CODE           No reset code

    A reset code has not been configured for the card or token.

GPG_ERR_BAD_RESET_CODE           No reset code

    A wrong reset code was used.

GPG_ERR_SYSTEM_BUG              System bug detected

   The underlying operating system misbehaved.  For example it wrote
   more to a buffer than the told maximum size.


GPG_ERR_DNS_UNKNOWN             Unknown DNS error

   Used by Dirmngr for DNS errors from libdns (DNS_EUNKNOWN);

GPG_ERR_DNS_SECTION             Invalid DNS section

   Used by Dirmngr for DNS errors from libdns (DNS_ESECTION);

GPG_ERR_DNS_ADDRESS             Invalid textual address form

   Used by Dirmngr for DNS errors from libdns (DNS_EADDRESS);

GPG_ERR_DNS_NO_QUERY            Missing DNS query packet

   Used by Dirmngr for DNS errors from libdns (DNS_ENOQUERY);

GPG_ERR_DNS_NO_ANSWER           Missing DNS answer packet

   Used by Dirmngr for DNS errors from libdns (DNS_ENOANSWER);

GPG_ERR_DNS_CLOSED              Connection closed in DNS

   Used by Dirmngr for DNS errors from libdns (DNS_ECONNFIN);

GPG_ERR_DNS_VERIFY              Verification failed in DNS

   Used by Dirmngr for DNS errors from libdns (DNS_EVERIFY);

GPG_ERR_DNS_TIMEOUT             DNS Timeout

   A DNS query timed out

GPG_ERR_LDAP_GENERAL            LDAP General error

   Catch all error for LDAP.  Use when an error code could not be
   mapped to a gpg-error code.

GPG_ERR_LDAP_ATTR_GENERAL       LDAP General attribute error
GPG_ERR_LDAP_NAME_GENERAL       LDAP General name error
GPG_ERR_LDAP_SECURITY_GENERAL   LDAP General security error
GPG_ERR_LDAP_SERVICE_GENERAL    LDAP General service error
GPG_ERR_LDAP_UPDATE_GENERAL     LDAP General update error
GPG_ERR_LDAP_E_GENERAL          LDAP Experimental error code
GPG_ERR_LDAP_X_GENERAL          LDAP Private error code
GPG_ERR_LDAP_OTHER_GENERAL      LDAP Other general error

  The 8 GPG_ERR_LDAP_*_GENERAL error codes may be used to map ranges
  of LDAP errors to one specific code. OpenLDAP uses LDAP_xxx_RANGE(n)
  macros for that mapping.  "Other general error" may be used similar
  to "General error" for mapping of ranges.  Here are macros from
  OpenLDAP for reference:

  #define LDAP_ATTR_ERROR(n)     LDAP_RANGE((n),0x10,0x15) /* 16-21 */
  #define LDAP_NAME_ERROR(n)     LDAP_RANGE((n),0x20,0x24) /* 32-34,36 */
  #define LDAP_SECURITY_ERROR(n) LDAP_RANGE((n),0x2F,0x32) /* 47-50 */
  #define LDAP_SERVICE_ERROR(n)  LDAP_RANGE((n),0x33,0x36) /* 51-54 */
  #define LDAP_UPDATE_ERROR(n)   LDAP_RANGE((n),0x40,0x47) /* 64-69,71 */
  #define LDAP_E_ERROR(n)        LDAP_RANGE((n),0x1000,0x3FFF)
  #define LDAP_X_ERROR(n)        LDAP_RANGE((n),0x4000,0xFFFF)

GPG_ERR_SQL_OK          SQL success

  This code is normally not used because it it mapped to GPG_ERR_NO_ERROR.

GPG_ERR_SQL_ERROR	SQL error

GPG_ERR_SQL_INTERNAL	Internal logic error in SQL library

GPG_ERR_SQL_PERM	Access permission denied (SQL)

GPG_ERR_SQL_ABORT	SQL abort was requested

GPG_ERR_SQL_BUSY	SQL database file is locked

GPG_ERR_SQL_LOCKED	An SQL table in the database is locked

GPG_ERR_SQL_NOMEM	SQL library ran out of core

GPG_ERR_SQL_READONLY	Attempt to write a readonly SQL database

GPG_ERR_SQL_INTERRUPT	SQL operation terminated by interrupt

GPG_ERR_SQL_IOERR	I/O error during SQL operation

GPG_ERR_SQL_CORRUPT	SQL database disk image is malformed

GPG_ERR_SQL_NOTFOUND	Unknown opcode in SQL file control

GPG_ERR_SQL_FULL	Insertion failed because SQL database is full

GPG_ERR_SQL_CANTOPEN	Unable to open the SQL database file

GPG_ERR_SQL_PROTOCOL	SQL database lock protocol error

GPG_ERR_SQL_EMPTY	(internal SQL code: empty)

GPG_ERR_SQL_SCHEMA	SQL database schema changed

GPG_ERR_SQL_TOOBIG	String or blob exceeds size limit (SQL)

GPG_ERR_SQL_CONSTRAINT	SQL abort due to constraint violation

GPG_ERR_SQL_MISMATCH	Data type mismatch (SQL)

GPG_ERR_SQL_MISUSE	SQL library used incorrectly

GPG_ERR_SQL_NOLFS	SQL library uses unsupported OS features

GPG_ERR_SQL_AUTH	Authorization denied (SQL)

GPG_ERR_SQL_FORMAT	(unused SQL code: format)

GPG_ERR_SQL_RANGE	SQL bind parameter out of range

GPG_ERR_SQL_NOTADB	File opened that is not an SQL database file

GPG_ERR_SQL_NOTICE	Notifications from SQL logger

GPG_ERR_SQL_WARNING	Warnings from SQL logger

GPG_ERR_SQL_ROW		SQL has another row ready

GPG_ERR SQL_DONE	SQL has finished executing


# Installed by libgpg-error 1.47-unknown
