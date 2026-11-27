/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "parallel/HighsTaskExecutor.h"

#include "parallel/HighsParallel.h"

#if defined(__linux__)
#include <sched.h>
#include <unistd.h>

#include <fstream>
#include <set>
#include <string>
#elif defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOGDI
#define NOGDI
#endif
#include <windows.h>

#include <algorithm>
#include <vector>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif

using namespace highs;

// Fallback: assume 2 vthreads per core
static unsigned int fallback_core_count() {
  return (std::thread::hardware_concurrency() + 1) / 2;
}

// Returns the number of physical cores available to this process,
// respecting the OS affinity mask.
unsigned int highs::parallel::available_core_count() {
#if defined(__linux__)
  // Use dynamically-allocated cpu_set to support systems with >1024 CPUs
  int num_cpus = sysconf(_SC_NPROCESSORS_CONF);
  if (num_cpus < 1) num_cpus = CPU_SETSIZE;
  size_t setsize = CPU_ALLOC_SIZE(num_cpus);
  cpu_set_t* affinity = CPU_ALLOC(num_cpus);
  if (affinity == nullptr) return fallback_core_count();

  if (sched_getaffinity(0, setsize, affinity) == 0) {
    // Count unique (package, core) pairs across the affinity set
    std::set<std::pair<int, int>> physical_cores;
    for (int cpu = 0; cpu < num_cpus; ++cpu) {
      if (!CPU_ISSET_S(cpu, setsize, affinity)) continue;
      std::string prefix =
          "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/";
      std::ifstream pkg_file(prefix + "physical_package_id");
      std::ifstream core_file(prefix + "core_id");
      if (!pkg_file.is_open() || !core_file.is_open()) {
        // Topology files unavailable; fall back
        CPU_FREE(affinity);
        return fallback_core_count();
      }
      int package_id, core_id;
      pkg_file >> package_id;
      core_file >> core_id;
      physical_cores.insert({package_id, core_id});
    }
    CPU_FREE(affinity);
    if (!physical_cores.empty())
      return static_cast<unsigned int>(physical_cores.size());
  } else {
    CPU_FREE(affinity);
  }
#elif defined(_WIN32) || defined(_WIN64)
  // Query which processor groups this process spans
  USHORT process_group_count = 0;
  std::vector<USHORT> process_groups;
  GetProcessGroupAffinity(GetCurrentProcess(), &process_group_count, nullptr);
  if (process_group_count <= 0) return fallback_core_count();
  process_groups.resize(process_group_count);
  if (!GetProcessGroupAffinity(GetCurrentProcess(), &process_group_count,
                               process_groups.data()))
    return fallback_core_count();

  // For single-group processes, get the affinity mask for precise filtering;
  // for multi-group, all bits set means no filtering.
  DWORD_PTR process_mask = ~static_cast<DWORD_PTR>(0);
  if (process_group_count == 1) {
    DWORD_PTR system_mask;
    GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask);
  }

  // Query topology across all groups via GetLogicalProcessorInformationEx
  DWORD length = 0;
  GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &length);
  if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) return fallback_core_count();

  // Variable-size entries packed in a byte buffer
  std::vector<char> buffer(length);
  auto entry_at = [&](DWORD offset) {
    return reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(
        buffer.data() + offset);
  };
  if (!GetLogicalProcessorInformationEx(RelationProcessorCore, entry_at(0),
                                        &length))
    return fallback_core_count();

  unsigned int physical_cores = 0;
  DWORD offset = 0;
  while (offset < length) {
    auto* entry = entry_at(offset);
    if (entry->Relationship == RelationProcessorCore) {
      for (WORD i = 0; i < entry->Processor.GroupCount; i++) {
        WORD group = entry->Processor.GroupMask[i].Group;
        KAFFINITY mask = entry->Processor.GroupMask[i].Mask;
        // Only count cores that match the affinity mask and process groups
        if ((mask & process_mask) &&
            std::find(process_groups.begin(), process_groups.end(), group) !=
                process_groups.end())
          physical_cores++;
      }
    }
    offset += entry->Size;
  }
  if (physical_cores > 0) return physical_cores;
#elif defined(__APPLE__)
  // macOS does not support affinity; return physical core count
  int n;
  size_t size = sizeof(n);
  if (sysctlbyname("hw.physicalcpu", &n, &size, nullptr, 0) == 0 && n > 0)
    return static_cast<unsigned int>(n);
#endif
  return fallback_core_count();
}

#ifdef _WIN32
static thread_local HighsSplitDeque* threadLocalWorkerDequePtr{nullptr};
HighsSplitDeque*& HighsTaskExecutor::threadLocalWorkerDeque() {
  return threadLocalWorkerDequePtr;
}

static thread_local HighsTaskExecutor::ExecutorHandle globalExecutorHandle{};

HighsTaskExecutor::ExecutorHandle&
HighsTaskExecutor::threadLocalExecutorHandle() {
  return globalExecutorHandle;
}
#else
thread_local HighsSplitDeque* HighsTaskExecutor::threadLocalWorkerDequePtr{
    nullptr};
thread_local HighsTaskExecutor::ExecutorHandle
    HighsTaskExecutor::globalExecutorHandle{};
#endif

void HighsTaskExecutor::ExecutorHandle::dispose() {
  if (ptr == nullptr) return;
  if (isMain) {
    ptr->stopWorkerThreads(false);
  }

  // check to see if we are the last handle and if so, delete the executor
  if (--ptr->referenceCount == 0) {
    cache_aligned::Deleter<HighsTaskExecutor>()(ptr);
  }

  ptr = nullptr;
}
