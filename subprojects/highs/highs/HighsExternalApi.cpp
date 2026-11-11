/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExternalApi.cpp
 * @brief Manages access (dynamic or static) to optional external dependencies
 */

#include "HighsExternalApi.h"

using namespace HighsExtras;

HighsExternalApi& HighsExternalApi::instance() {
  static HighsExternalApi _instance;
  return _instance;
}

void HighsExternalApi::unload() {
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
  HighsExternalApi& inst = instance();

  inst.library_.unload();
  inst.api_ = HighsExtrasApi{};
  inst.extras_feature_info_ = nullptr;
  inst.available_ = false;
  inst.status_.clear();
#endif
}

void HighsExternalApi::reset() {
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
  HighsExternalApi& inst = instance();
  inst.initialized_.reset(new std::once_flag);
#endif
}

#ifdef HIGHS_SHARED_EXTRAS_LIBRARY

#define STRINGFY(s) STRINGFY0(s)
#define STRINGFY0(s) #s

bool HighsExternalApi::tryLoad(const std::string& directory,
                               bool search_if_empty) {
  HighsExternalApi& inst = instance();
  // prevents multiple attempts to load the library
  // ensure thread safety (multiple threads may call tryLoad simultaneously)
  std::call_once(*inst.initialized_, [&]() {
    inst.available_ = false;

    // Load library
#if defined(_WIN32) || defined(_WIN64)
    const char* library_filename = "highs_extras.dll";
#elif defined(__APPLE__)
    const char* library_filename = "libhighs_extras.dylib";
#else
    const char* library_filename = "libhighs_extras.so";
#endif

    bool ok = (search_if_empty || directory.empty() == false) &&
              inst.library_.load(library_filename, directory);

    inst.status_ = "Extras: " + inst.library_.status();

    if (ok) {
      try {
        // Check ABI compatibility
        std::string highs_version = STRINGFY(HIGHS_VERSION_MAJOR) "." STRINGFY(
            HIGHS_VERSION_MINOR) "." STRINGFY(HIGHS_VERSION_PATCH);

        std::string extras_version =
            inst.library_.call<get_version_t>("HighsExtras_getVersion");

        if (extras_version == highs_version) {
          inst.extras_feature_info_ = inst.library_.call<get_feature_info_t>(
              "HighsExtras_getFeatureInfo");

          inst.library_.call<get_api_t>("HighsExtras_getApi", &inst.api_);
        } else {
          inst.status_ = "Extras: ABI version mismatch: expected " +
                         highs_version + ", got " + extras_version + ".";
          ok = false;
        }
      } catch (const std::exception& e) {
        inst.status_ = std::string("Extras: failed to load API: ") + e.what();
        ok = false;
      }

      if (!ok) {
        inst.unload();
      }
    }

    inst.available_ = ok;
  });

  return inst.available_;
}
#else
bool HighsExternalApi::tryLoad(const std::string& directory,
                               bool search_if_empty) {
  (void)directory;
  static std::once_flag flag;

  std::call_once(flag, []() {
    HighsExternalApi& inst = instance();
    HighsExtras_getApi(&inst.api_);
    inst.extras_feature_info_ = HighsExtras_getFeatureInfo();
    inst.status_ = "Extras: Built in";
  });

  return true;
}
#endif
