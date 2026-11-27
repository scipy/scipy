/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsAvailableMemory.cpp
 * @brief Memory available to the process
 */
#include "util/HighsAvailableMemory.h"

#include <algorithm>
#include <atomic>
#include <fstream>
#include <limits>
#include <string>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#endif

#if defined(__linux__)
// First number in a file, or negative if there is none, as for "max"
static double readNumber(const std::string& path) {
  std::ifstream file(path);
  double value;
  return file >> value ? value : -1;
}

// Value of a "key value" line in a file, or zero if there is none
static double readKeyValue(const std::string& path, const std::string& key) {
  std::ifstream file(path);
  std::string name;
  double value;
  while (file >> name >> value) {
    if (name == key) return value;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return 0;
}

// Memory left under the cgroup limits on this process, which are what a
// container or a systemd unit is held to. Limits apply at every level
static double cgroupAvailableBytes() {
  std::string v1_path, v2_path;
  bool has_v1 = false, has_v2 = false;
  std::ifstream self("/proc/self/cgroup");
  std::string line;
  while (std::getline(self, line)) {
    const size_t first = line.find(':');
    const size_t second = line.find(':', first + 1);
    if (first == std::string::npos || second == std::string::npos) continue;
    const std::string controllers = line.substr(first + 1, second - first - 1);
    if (controllers.empty()) {
      has_v2 = true;
      v2_path = line.substr(second + 1);
    } else if (("," + controllers + ",").find(",memory,") !=
               std::string::npos) {
      has_v1 = true;
      v1_path = line.substr(second + 1);
    }
  }

  double available = std::numeric_limits<double>::infinity();
  auto walk = [&](std::string path, const std::string& root,
                  const std::string& limit_file, const std::string& usage_file,
                  const std::string& inactive_key) {
    // Inside a container the path may not exist under root, and root is then
    // the container's own cgroup, so reach it however the walk goes
    while (true) {
      const std::string dir = root + path + "/";
      const double limit = readNumber(dir + limit_file);
      const double usage = readNumber(dir + usage_file);
      if (limit >= 0 && usage >= 0) {
        // Inactive page cache is reclaimed before the limit is enforced
        const double reclaimable =
            readKeyValue(dir + "memory.stat", inactive_key);
        available =
            std::min(available, std::max(0.0, limit - (usage - reclaimable)));
      }
      const size_t slash = path.rfind('/');
      if (slash == std::string::npos) break;
      path.erase(slash);
    }
  };
  if (has_v2)
    walk(v2_path, "/sys/fs/cgroup", "memory.max", "memory.current",
         "inactive_file");
  if (has_v1)
    walk(v1_path, "/sys/fs/cgroup/memory", "memory.limit_in_bytes",
         "memory.usage_in_bytes", "total_inactive_file");
  return available;
}
#endif

// Atomic, since HiPO may check memory on other threads
static std::atomic<double> available_memory_bytes_for_testing{-1};

void setAvailableMemoryBytesForTesting(const double bytes) {
  available_memory_bytes_for_testing = bytes;
}

bool availableMemoryBytesSetForTesting() {
  return available_memory_bytes_for_testing >= 0;
}

double availableMemoryBytes() {
  if (available_memory_bytes_for_testing >= 0)
    return available_memory_bytes_for_testing;
#if defined(_WIN32)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  return GlobalMemoryStatusEx(&status)
             ? static_cast<double>(status.ullAvailPhys)
             : -1;
#elif defined(__APPLE__)
  vm_statistics64_data_t stats;
  mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
  if (host_statistics64(mach_host_self(), HOST_VM_INFO64,
                        reinterpret_cast<host_info64_t>(&stats),
                        &count) != KERN_SUCCESS)
    return -1;
  return static_cast<double>(stats.free_count + stats.inactive_count) *
         vm_kernel_page_size;
#elif defined(__linux__)
  const double available =
      readKeyValue("/proc/meminfo", "MemAvailable:") * 1024;
  if (available <= 0) return -1;
  return std::min(available, cgroupAvailableBytes());
#else
  return -1;
#endif
}
