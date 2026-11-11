/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsDynamicLibrary.h
 * @brief Lightweight cross-platform runtime shared library loader
 */
#ifndef HIGHS_UTIL_HIGHSDYNAMICLIBRARY_H_
#define HIGHS_UTIL_HIGHSDYNAMICLIBRARY_H_

#include <stdexcept>
#include <string>
#include <type_traits>

class HighsDynamicLibrary {
 public:
  HighsDynamicLibrary() = default;
  ~HighsDynamicLibrary();

  HighsDynamicLibrary(const HighsDynamicLibrary&) = delete;
  HighsDynamicLibrary& operator=(const HighsDynamicLibrary&) = delete;

  bool load(const std::string& filename, const std::string& path = "");
  void unload();

  bool isLoaded() const { return handle_ != nullptr; }

  const std::string& status() const { return status_; }

  // resolve a symbol to a function pointer of the given type
  template <typename FuncType>
  bool resolve(FuncType& target, const char* name) const {
    target = nullptr;
    if (!handle_) return false;

    target = reinterpret_cast<FuncType>(resolveRaw(name));
    return target != nullptr;
  }

  // resolve and call a function directly (throws on failure)
  template <typename FuncType, typename... Args>
  auto call(const char* name, Args... args) const
      -> decltype(std::declval<FuncType>()(std::declval<Args>()...)) {
    FuncType func;
    if (!resolve(func, name)) {
      throw std::runtime_error(std::string("Failed to resolve symbol: ") +
                               name);
    }
    return func(args...);
  }

 private:
  void* resolveRaw(const char* name) const;

  void* handle_ = nullptr;
  std::string status_;
};

#endif  // HIGHS_UTIL_HIGHSDYNAMICLIBRARY_H_