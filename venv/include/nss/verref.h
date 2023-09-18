/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/* This header is used inline in a function to ensure that a version string
 * symbol is linked in and not optimized out. A volatile reference is added to
 * the variable identified by NSS_VERSION_VARIABLE.
 *
 * Use this as follows:
 *
 * #define NSS_VERSION_VARIABLE __nss_ssl_version
 * #include "verref.h"
 */

/* Suppress unused variable warnings. */
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4101)
#endif
/* This works for both gcc and clang */
#if defined(__GNUC__) && !defined(NSS_NO_GCC48)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#ifndef NSS_VERSION_VARIABLE
#error NSS_VERSION_VARIABLE must be set before including "verref.h"
#endif
{
    extern const char NSS_VERSION_VARIABLE[];
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((unused))
#endif
    volatile const char _nss_version_c = NSS_VERSION_VARIABLE[0];
}
#undef NSS_VERSION_VARIABLE

#ifdef _MSC_VER
#pragma warning(pop)
#endif
#if defined(__GNUC__) && !defined(NSS_NO_GCC48)
#pragma GCC diagnostic pop
#endif
