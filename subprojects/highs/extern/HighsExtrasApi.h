/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasApi.h
 * @brief API for the HiGHS external deps library for dynamic loading
 */
#ifndef HIGHS_EXTRAS_API_H_
#define HIGHS_EXTRAS_API_H_

#include <utility>

#include "HConfig.h"
#include "HighsExtrasExternalDeps.h"

// ABI note:
// HighsExtras_getApi uses C linkage for symbol lookup only.
// The HighsExtrasApi parameter is a C++ type, so the caller and callee
// must be built with ABI-compatible toolchains. In practice, this means
// using the same compiler/toolset (or a known ABI-compatible equivalent).

// Building / importing / static linking is determined by build macros:
#if defined(_WIN32) || defined(_WIN64)
#define HIGHS_EXTRAS_EXPORT __declspec(dllexport)
#define HIGHS_EXTRAS_IMPORT __declspec(dllimport)
#else
#define HIGHS_EXTRAS_EXPORT __attribute__((visibility("default")))
#define HIGHS_EXTRAS_IMPORT __attribute__((visibility("default")))
#endif

#if defined(HIGHS_EXTRAS_LIBRARY_BUILD)
#define HIGHS_EXTRAS_API HIGHS_EXTRAS_EXPORT
#elif defined(HIGHS_SHARED_EXTRAS_LIBRARY)
#define HIGHS_EXTRAS_API HIGHS_EXTRAS_IMPORT
#else
#define HIGHS_EXTRAS_API
#endif

// C++ API with actual types (outside extern "C" block)
// These use C++ references and HiGHS types directly
#ifdef __cplusplus

// Note: We use extern "C" here to disable C++ name mangling, allowing
// GetProcAddress/dlsym to find the function by its simple name.
//
// While extern "C" is typically for C-compatible interfaces,
// using C++ types (references, classes) works here because:
// 1. Both highspy and highspy-extras are built with the same C++ compiler/ABI
// 2. The ABI version check ensures struct layouts match between builds
// 3. We only need C linkage for symbol lookup, not C type compatibility
extern "C" {

/**
 * Get the version string of the HiGHS extras library, used to verify
 * compatibility with the main HiGHS library.
 *
 * @return Version string (e.g., "1.16.0")
 */
HIGHS_EXTRAS_API const char* HighsExtras_getVersion(void);

// Get metadata for all features
HIGHS_EXTRAS_API size_t HighsExtras_getFeatureCount(void);
HIGHS_EXTRAS_API const char* HighsExtras_getFeatureName(size_t index);
HIGHS_EXTRAS_API const HighsExtrasFeatureInfo* HighsExtras_getFeatureInfo();

// Get the HighsExtrasApi, with appropriate function pointers
HIGHS_EXTRAS_API bool HighsExtras_getApi(HighsExtras::HighsExtrasApi* api);

}  // extern "C"

using get_version_t = decltype(&HighsExtras_getVersion);
using get_feature_info_t = decltype(&HighsExtras_getFeatureInfo);
using get_api_t = decltype(&HighsExtras_getApi);

#endif  // __cplusplus

#endif  // HIGHS_EXTRAS_API_H_
