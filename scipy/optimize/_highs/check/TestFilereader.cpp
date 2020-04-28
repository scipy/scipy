#include <cstdio>

#include "Highs.h"
#include "catch.hpp"
#include "io/FilereaderEms.h"
#include "io/HMPSIO.h"
#include "io/HMpsFF.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"

TEST_CASE("filereader-edge-cases", "[highs_filereader]") {
  std::string model = "";
  std::string model_file;
  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;
  HighsOptions options;

  // Several tests don't pass, but should, so possibly skip them
  const bool test_garbage_mps = false;
  const bool test_garbage_ems = true;
  const bool test_garbage_lp = false;

  Highs highs(options);
  const HighsInfo& info = highs.getHighsInfo();

  // Try to run HiGHS with default options. No model loaded so fails
  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::Error);

  // Set model_file to non-existent file and try to run HiGHS
  model = "";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::Error);

  // Set model_file to non-supported file type and try to run HiGHS
  model = "model";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".xyz";
  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::Error);

  // Set model_file to existing MPS file and run HiGHS
  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::OK);
  REQUIRE(info.simplex_iteration_count == 86);

  model = "garbage";
  if (test_garbage_mps) {
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    return_status = highs.setHighsOptionValue("model_file", model_file);
    REQUIRE(return_status == HighsStatus::OK);

    read_status = highs.readModel(model_file);
    REQUIRE(read_status == HighsStatus::Error);
  }

  if (test_garbage_ems) {
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".ems";
    return_status = highs.setHighsOptionValue("model_file", model_file);
    REQUIRE(return_status == HighsStatus::OK);

    read_status = highs.readModel(model_file);
    REQUIRE(read_status == HighsStatus::Error);
  }

  if (test_garbage_lp) {
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
    return_status = highs.setHighsOptionValue("model_file", model_file);
    REQUIRE(return_status == HighsStatus::OK);

    read_status = highs.readModel(model_file);
    REQUIRE(read_status == HighsStatus::Error);
  }
}
TEST_CASE("filereader-free-format-parser", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsStatus status;

  // Read mps
  HighsOptions options;

  Highs highs(options);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_free = highs.getLp();

  status = highs.setHighsOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::OK);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_fixed = highs.getLp();

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

// No commas in test case name.
TEST_CASE("filereader-read-mps-ems-lp", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsStatus status;

  // Read mps
  HighsOptions options;

  Highs highs(options);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);
  HighsLp lp_mps = highs.getLp();

  // Write lp
  std::string filename_lp = "adlittle.lp";
  status = highs.writeModel(filename_lp);
  REQUIRE(status == HighsStatus::OK);

  /*
  bool are_the_same;
  // Write ems
  std::string filename_ems = "adlittle.ems";
  status = highs.writeModel(filename_ems);
  REQUIRE(status == HighsStatus::OK);

  // Read ems and compare with mps
  std::cout << "Reading " << filename_ems << std::endl;
  status = highs.readModel(filename_ems);
  REQUIRE(status == HighsStatus::OK);

  std::cout << "Compare LP from .ems and .mps" << std::endl;
  are_the_same = lp_mps == highs.getLp();
  REQUIRE(are_the_same);

  std::remove(filename_ems.c_str());
  */

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  const HighsInfo& info = highs.getHighsInfo();
  double mps_objective_function_value = info.objective_function_value;

  // Read lp and compare objective with mps
  std::cout << "Reading " << filename_lp << std::endl;
  status = highs.readModel(filename_lp);
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  REQUIRE(mps_objective_function_value == info.objective_function_value);

  std::remove(filename_lp.c_str());
}

TEST_CASE("filereader-integrality-constraints", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/small_mip.mps";

  // integer variables are COL03,COL04 so x[2], x[3].
  const std::vector<int> kIntegers{0, 0, 1, 1, 0, 0, 0, 0};

  HighsStatus status;
  HighsOptions options;

  Highs highs(options);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_free = highs.getLp();

  REQUIRE(lp_free.integrality_.size() == lp_free.numCol_);
  REQUIRE(lp_free.integrality_ == kIntegers);

  // Read mps with fixed format parser.
  status = highs.setHighsOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::OK);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_fixed = highs.getLp();

  REQUIRE(lp_fixed.integrality_.size() == lp_fixed.numCol_);
  REQUIRE(lp_fixed.integrality_ == kIntegers);

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

TEST_CASE("filereader-dualize", "[highs_data]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  // Read mps.
  HighsOptions options;
  options.model_file = filename;

  HighsLp lp;
  HMpsFF parser{};
  FreeFormatParserReturnCode result = parser.loadProblem(stdout, filename, lp);
  REQUIRE(result == FreeFormatParserReturnCode::SUCCESS);

  HighsLp primal;
  HighsStatus status;
  status = transformIntoEqualityProblem(lp, primal);
  REQUIRE(status == HighsStatus::OK);

  Highs highs_lp;
  HighsModelStatus model_status;
  status = highs_lp.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  status = highs_lp.run();
  model_status = highs_lp.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  Highs highs_primal;
  status = highs_primal.passModel(primal);
  REQUIRE(status == HighsStatus::OK);
  status = highs_primal.run();
  model_status = highs_lp.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  double lp_objective;
  highs_lp.getHighsInfoValue("objective_function_value", lp_objective);
  double primal_objective;
  highs_lp.getHighsInfoValue("objective_function_value", primal_objective);

  double diff_equality = lp_objective - primal_objective;
  REQUIRE(diff_equality < 0.00000001);

  HighsLp dual;
  status = dualizeEqualityProblem(primal, dual);
  REQUIRE(status == HighsStatus::OK);
  Highs highs_dual;
  status = assessLp(dual, options);
  REQUIRE(status == HighsStatus::OK);
  status = highs_dual.passModel(dual);
  REQUIRE(status == HighsStatus::OK);
  status = highs_dual.run();
  model_status = highs_lp.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  double dual_objective;
  highs_dual.getHighsInfoValue("objective_function_value", dual_objective);

  double diff_dual = primal_objective + dual_objective;
  REQUIRE(diff_dual < 0.00000001);
}
