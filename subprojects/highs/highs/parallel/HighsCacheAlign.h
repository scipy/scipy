/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_CACHE_ALIGN_H_
#define HIGHS_CACHE_ALIGN_H_

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>

#include "HighsSchedulerConstants.h"

namespace highs {

struct cache_aligned {
  static constexpr std::size_t alignment() {
    return HighsSchedulerConstants::kCacheLineSize;
  }

  static void* alloc(std::size_t size) {
    using std::uintptr_t;

    uintptr_t ptr =
        reinterpret_cast<uintptr_t>(::operator new(size + alignment()));

    char* aligned_ptr = reinterpret_cast<char*>((ptr | (alignment() - 1)) + 1);
    std::memcpy(aligned_ptr - sizeof(uintptr_t), &ptr, sizeof(uintptr_t));
    return reinterpret_cast<void*>(aligned_ptr);
  }

  static void free(void* aligned_ptr) {
    using std::uintptr_t;

    if (aligned_ptr) {
      void* freeptr;
      std::memcpy(&freeptr,
                  reinterpret_cast<char*>(aligned_ptr) - sizeof(uintptr_t),
                  sizeof(uintptr_t));
      ::operator delete(freeptr);
    }
  }

  template <typename T>
  struct Deleter {
    void operator()(T* ptr) const {
      ptr->~T();
      free(ptr);
    }
  };

  template <typename T>
  struct Deleter<T[]> {
    void operator()(T* ptr) const { free(ptr); }
  };

  template <typename T>
  using unique_ptr = std::unique_ptr<T, Deleter<T>>;

  template <typename T>
  using shared_ptr = std::shared_ptr<T>;

  template <typename T, typename... Args>
  static shared_ptr<T> make_shared(Args&&... args) {
    return shared_ptr<T>(new (alloc(sizeof(T))) T(std::forward<Args>(args)...),
                         Deleter<T>());
  }

  template <typename T, typename... Args>
  static unique_ptr<T> make_unique(Args&&... args) {
    return unique_ptr<T>(new (alloc(sizeof(T))) T(std::forward<Args>(args)...));
  }

  template <typename T, typename... Args>
  static unique_ptr<T[]> make_unique_array(std::size_t N) {
    return unique_ptr<T[]>(static_cast<T*>(alloc(sizeof(T) * N)));
  }
};

}  // namespace highs

#endif
