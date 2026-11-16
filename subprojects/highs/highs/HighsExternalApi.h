/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExternalApi.h
 * @brief Manages access (dynamic or static) to optional external dependencies
 */
#ifndef HIGHS_EXTERNAL_API_H_
#define HIGHS_EXTERNAL_API_H_
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>

#include "HConfig.h"
#include "HighsExternalDeps.h"
#include "HighsExtrasApi.h"
#include "io/HighsIO.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsDynamicLibrary.h"
#include "util/HighsType.h"
#include "util/stringutil.h"

//
// Support nested multiple features, e.g.
// - isAvailable<Feature, ...>
// - isAvailable<require<Feature1, Feature2>, Feature3, ...>
// - isAvailable<require<Feature1, require<Feature2, Feature3>>, Feature4, ...>
//
// Defines support for isAvailable, notice rows, and missing feature names
//

namespace HighsExtras {

// provides the actual feature trait implementation
template <class Feature>
struct feature_ops {
  static bool is_available() {
    const auto* info = Feature::getInfo();
    return info && info->enabled;
  }

  static void append_notice_rows(HighsTextTable<4>& table) {
    if (is_available()) {
      const auto* info = Feature::getInfo();

      HighsTextTable<4>::Row row = {
          {Feature::name(), info->provider, info->version, info->license}};
      table.addRow(row);
    }
  }

  static void append_missing_names(std::set<std::string>& names) {
    if (!is_available()) names.insert(Feature::name());
  }
};

// forward declarations
template <class... Traits>
struct trait_pack_ops;

template <>
struct trait_pack_ops<> {
  static bool all_available() { return true; }
  static void append_notice_rows(HighsTextTable<4>&) {}
  static void append_missing_names(std::set<std::string>&) {}
};

template <class T>
struct feature_traits : feature_ops<T> {};

// handles features or nested feature sets wrapped in require<...>
// and flattens to a list of traits for trait_pack_ops to handle
template <class... Features>
struct feature_traits<require<Features...>> {
  using pack = trait_pack_ops<feature_traits<Features>...>;

  static bool is_available() { return pack::all_available(); }

  static void append_notice_rows(HighsTextTable<4>& table) {
    pack::append_notice_rows(table);
  }

  static void append_missing_names(std::set<std::string>& names) {
    pack::append_missing_names(names);
  }
};

// handles a list of features (not wrapped in require<...>)
// and recursively calls the actual trait implementation for each feature
template <class Trait, class... Rest>
struct trait_pack_ops<Trait, Rest...> {
  static bool all_available() {
    return Trait::is_available() && trait_pack_ops<Rest...>::all_available();
  }

  static void append_notice_rows(HighsTextTable<4>& table) {
    Trait::append_notice_rows(table);
    trait_pack_ops<Rest...>::append_notice_rows(table);
  }

  static void append_missing_names(std::set<std::string>& names) {
    Trait::append_missing_names(names);
    trait_pack_ops<Rest...>::append_missing_names(names);
  }
};

}  // namespace HighsExtras

/**
 * External dependencies for HiGHS that can either be dynamically loaded
 * or linked at compile time, depending on the build configuration.
 *
 * This class handles runtime / static / missing external dependencies
 * those that have a different license from HiGHS.
 */
struct HighsExternalApi {
 public:
  HighsExternalApi() = default;
  ~HighsExternalApi() { unload(); }

  // Prevent copying
  HighsExternalApi(const HighsExternalApi&) = delete;
  HighsExternalApi& operator=(const HighsExternalApi&) = delete;

  static HighsExternalApi& instance();
  HighsExtras::HighsExtrasApi api_;
  const HighsExtrasFeatureInfo* extras_feature_info_ = nullptr;

  static void unload();  // unloads but does not allow reloading without reset()
  static void reset();   // resets the cached initialization attempt for
                         // dynamically loaded extras
  static bool tryLoad(const std::string& directory = "",
                      bool search_if_empty = true);

  static const std::string getLoadStatus() { return instance().status_; }

  // define <feature> functions, e.g. isAvailable<Feature1, Feature2, ...>()
  template <class Feature, class... Features>
  static bool isAvailable() {
    using traits =
        HighsExtras::trait_pack_ops<HighsExtras::feature_traits<Feature>,
                                    HighsExtras::feature_traits<Features>...>;
    tryLoad();
    return traits::all_available();
  }

  static std::string thirdPartyNoticeHeader() {
    return "Includes third-party software components, see "
           "THIRD_PARTY_NOTICES.md for full details";
  }

  template <class Feature, class... Features>
  static std::string getThirdPartyNotice() {
    using traits =
        HighsExtras::trait_pack_ops<HighsExtras::feature_traits<Feature>,
                                    HighsExtras::feature_traits<Features>...>;
    tryLoad();
    HighsTextTable<4>::Row headers{{"key", "name", "version", "license"}};
    HighsTextTable<4> table(headers);
    traits::append_notice_rows(table);

    table.sortRows<0>();  // sort by key for consistent ordering
    return std::string("Third-party components:\n\n") + table.render();
  }

  template <class Feature, class... Features>
  static std::string getMissingFeatures() {
    using traits =
        HighsExtras::trait_pack_ops<HighsExtras::feature_traits<Feature>,
                                    HighsExtras::feature_traits<Features>...>;
    tryLoad();
    std::set<std::string> names;
    traits::append_missing_names(names);

    return joinString(names.begin(), names.end(), std::string(", "));
  }

  template <class Feature, class... Features>
  static void logUnavailable(const HighsLogOptions& log_options,
                             const HighsLogType type, const char* format = "",
                             ...) {
    if (isAvailable<Feature, Features...>()) return;

    va_list argptr;
    va_start(argptr, format);
    std::array<char, kIoBufferSize> msgbuffer = {};
    const int len =
        vsnprintf(msgbuffer.data(), msgbuffer.size(), format, argptr);
    va_end(argptr);

    const std::string load_failure = tryLoad() ? "" : ". " + getLoadStatus();
    highsLogUser(log_options, type,
                 "%s%sThe following features are unavailable: %s%s\n",
                 msgbuffer.data(), (len > 0) ? " " : "",
                 getMissingFeatures<Feature, Features...>().c_str(),
                 load_failure.c_str());
  }

 private:
  std::string status_;

#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
  HighsDynamicLibrary library_;
  std::unique_ptr<std::once_flag> initialized_{new std::once_flag};
  bool available_ = false;
#endif
};

namespace HighsExtras {

inline HighsStatus initialize(const std::string& path,
                              bool search_default_dirs = true) {
  return HighsExternalApi::tryLoad(path, search_default_dirs)
             ? HighsStatus::kOk
             : HighsStatus::kError;
}

// allows static access to feature APIs via
// HighsExtras::<feature>::<method>(...)
template <class Methods>
feature_api<Methods>& wrapper_storage<extras_family>::getApi() {
  return HighsExternalApi::instance().api_.template as<Methods>();
}

inline const HighsExtrasFeatureInfo* wrapper_storage<extras_family>::getInfo() {
  return HighsExternalApi::instance().extras_feature_info_;
}

}  // namespace HighsExtras

#endif  // HIGHS_EXTERNAL_API_H_
