#include <algorithm>
#include <vector>

#include "catch.hpp"
#include "lp_data/HConst.h"
#include "util/HighsRandom.h"
#include "util/HighsSort.h"

// No commas in test case name.
TEST_CASE("HiGHS_sort", "[highs_data]") {
  int num_values = 10;
  std::vector<int> indices;
  std::vector<int> int_values;
  std::vector<double> double_values;
  std::vector<double> original_double_values;
  indices.resize(1 + num_values);
  int_values.resize(num_values);
  double_values.resize(1 + num_values);
  original_double_values.resize(1 + num_values);

  // Set up a vector of random number and their corresponding indices
  HighsRandom random;
  for (int ix = 0; ix < num_values; ix++) {
    double_values[1 + ix] = random.fraction();
    original_double_values[1 + ix] = double_values[1 + ix];
    indices[1 + ix] = ix;
  }
  // Sort the vector of random number and their corresponding indices
  maxheapsort(&double_values[0], &indices[0], num_values);

  // Check that the random numbers are ascending and that the indices
  // point from the original values to their new positions
  bool error0 = false;
  bool error1 = false;
  double previous_double = -1e200;
  for (int ix = 0; ix < num_values; ix++) {
    //    printf("%2d: %2d %12g %12g\n", ix, indices[1+ix], double_values[1+ix],
    //    original_double_values[1+ix]);
    error0 = error0 || double_values[1 + ix] < previous_double;
    previous_double = double_values[1 + ix];
    error1 = error1 ||
             double_values[1 + ix] == original_double_values[indices[1 + ix]];
  }

  REQUIRE(error0 == false);
  REQUIRE(error1 == false);

  // Use the indices of the previous sort as a vector of integers to sort
  for (int ix = 0; ix < num_values; ix++) {
    double_values[ix] = double_values[ix + 1];
    int_values[ix] = indices[1 + ix];
  }
  std::make_heap(int_values.begin(), int_values.end());
  std::sort_heap(int_values.begin(), int_values.end());
  //  maxheapsort(&int_values[0], num_values);

  bool ok;
  // Check that the values in the vector of doubles are ascending - can do
  // strict test
  ok = increasing_set_ok(&double_values[0], num_values, 0, 1, true);
  REQUIRE(ok == true);

  // Check that the values in the vector of integers are ascending - maybe can't
  // do strict test
  ok = increasing_set_ok(&int_values[0], num_values, 0, num_values, false);
  REQUIRE(ok == true);

  num_values = 14;
  std::vector<int> set;
  std::vector<double> lb;
  std::vector<double> ub;
  set.resize(num_values);
  lb.resize(num_values);
  ub.resize(num_values);

  set[0] = 6;
  lb[0] = 60;
  ub[0] = 60;
  set[1] = 7;
  lb[1] = 6;
  ub[1] = 6;
  set[2] = 8;
  lb[2] = 60;
  ub[2] = 60;
  set[3] = 13;
  lb[3] = 600;
  ub[3] = 1200;
  set[4] = 4;
  lb[4] = 70;
  ub[4] = 70;
  set[5] = 5;
  lb[5] = 16;
  ub[5] = 16;
  set[6] = 2;
  lb[6] = 70;
  ub[6] = 70;
  set[7] = 3;
  lb[7] = 7;
  ub[7] = 7;
  set[8] = 11;
  lb[8] = 200;
  ub[8] = 1400;
  set[9] = 0;
  lb[9] = 75;
  ub[9] = 75;
  set[10] = 1;
  lb[10] = 12;
  ub[10] = 12;
  set[11] = 14;
  lb[11] = 0;
  ub[11] = 1400;
  set[12] = 9;
  lb[12] = 6;
  ub[12] = 6;
  set[13] = 15;
  lb[13] = 600;
  ub[13] = 1200;

  std::vector<int> sorted_set;
  std::vector<double> sorted_lb;
  std::vector<double> sorted_ub;
  sorted_set.resize(num_values);
  sorted_lb.resize(num_values);
  sorted_ub.resize(num_values);

  sortSetData(num_values, &set[0], &lb[0], &ub[0], NULL, &sorted_set[0],
              &sorted_lb[0], &sorted_ub[0], NULL);

  int prev_ix = -HIGHS_CONST_I_INF;
  for (int k0 = 0; k0 < num_values; k0++) {
    int ix = sorted_set[k0];
    REQUIRE(ix >= prev_ix);
    int k1 = -HIGHS_CONST_I_INF;
    for (int k_1 = 0; k_1 < num_values; k_1++) {
      if (set[k_1] == ix) {
        k1 = k_1;
        break;
      }
    }
    REQUIRE(k1 > -HIGHS_CONST_I_INF);
    REQUIRE(sorted_lb[k0] == lb[k1]);
    REQUIRE(sorted_ub[k0] == ub[k1]);
  }
}
