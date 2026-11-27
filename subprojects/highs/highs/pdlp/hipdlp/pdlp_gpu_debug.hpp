#ifndef PDLP_GPU_DEBUG_HPP__
#define PDLP_GPU_DEBUG_HPP__

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#if PDLP_DEBUG_LOG

bool vecDiff(const std::vector<double>& x1, const std::vector<double>& x2,
             double tol, const std::string& func_name = "vecDiff") {
  if (x1.size() != x2.size()) {
    printf("[DEBUG %s] FAILED: Vector sizes do not match! CPU=%zu, GPU=%zu\n",
           func_name.c_str(), x1.size(), x2.size());
    return false;
  }
  if (x1.empty()) {
    return true;
  }

  HighsInt num_diffs = 0;
  double max_diff = 0.0;
  HighsInt max_diff_idx = -1;

  for (size_t i = 0; i < x1.size(); ++i) {
    double diff = std::abs(x1[i] - x2[i]);
    if (diff > max_diff) {
      max_diff = diff;
      max_diff_idx = (HighsInt)i;
    }
    if (diff > tol) {
      num_diffs++;
      if (num_diffs <= 10) {
        // THIS IS THE PRINTF YOU ARE SEEING
        printf("[DEBUG %s] Difference at index %zu: %.4f vs %.4f\n",
               func_name.c_str(), i, x1[i], x2[i]);
      }
    }
  }

  if (num_diffs > 0) {
    if (num_diffs > 10) {
      printf("[DEBUG %s] ...and %d more differences found.\n",
             func_name.c_str(), (num_diffs - 10));
    }
    printf(
        "[DEBUG %s] SUMMARY: FAILED. %d/%zu elements differ. Max diff=%.8e at "
        "idx %d\n",
        func_name.c_str(), num_diffs, x1.size(), max_diff, max_diff_idx);
    return false;
  }
  // This will print a success message if it passes
  // printf("[DEBUG %s] SUCCESS: All %zu elements within tol=%.1e. Max
  // diff=%.8e\n",
  //        func_name.c_str(), x1.size(), tol, max_diff);
  return true;
}

#endif

#endif  // PDLP_GPU_DEBUG_HPP__
