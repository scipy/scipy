#include <cstdio>

#include "Highs.h"
#include "catch.hpp"
#include "io/FilereaderEms.h"
#include "io/HMPSIO.h"
#include "io/HMpsFF.h"
#include "io/HighsIO.h"
#include "io/LoadProblem.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"

TEST_CASE("free-format-parser", "[highs_filereader]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  // Read mps.
  HighsLp lp_free_format, lp_fixed_format;
  bool are_the_same = false;

  std::vector<int> integerColumn;
  FilereaderRetcode status = readMPS(
      stdout, filename.c_str(), -1, -1, lp_fixed_format.numRow_,
      lp_fixed_format.numCol_, lp_fixed_format.numInt_, lp_fixed_format.sense_,
      lp_fixed_format.offset_, lp_fixed_format.Astart_, lp_fixed_format.Aindex_,
      lp_fixed_format.Avalue_, lp_fixed_format.colCost_,
      lp_fixed_format.colLower_, lp_fixed_format.colUpper_,
      lp_fixed_format.rowLower_, lp_fixed_format.rowUpper_, integerColumn,
      lp_fixed_format.col_names_, lp_fixed_format.row_names_);
  lp_fixed_format.nnz_ = lp_fixed_format.Avalue_.size();
  if (status == FilereaderRetcode::OK) {
    HMpsFF parser{};
    FreeFormatParserReturnCode result =
        parser.loadProblem(stdout, filename, lp_free_format);
    if (result != FreeFormatParserReturnCode::SUCCESS)
      status = FilereaderRetcode::PARSERERROR;
    if (status == FilereaderRetcode::OK)
      are_the_same = lp_free_format == lp_fixed_format;
  }

  // In case you want to compare.
  // FilereaderEms ems;
  // ems.writeModelToFile(options, "fixed.ems", lp_fixed_format);
  // ems.writeModelToFile(options, "free.ems", lp_free_format);

  REQUIRE(are_the_same);
}

// No commas in test case name.
TEST_CASE("read-mps-ems", "[highs_filereader]") {
  HighsOptions options;

  std::cout << std::string(HIGHS_DIR) << std::endl;

  options.model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  // Read mps.
  HighsLp lp_mps;
  HighsStatus read_status = loadLpFromFile(options, lp_mps);
  REQUIRE(read_status == HighsStatus::OK);

  // Write ems.
  FilereaderEms ems;
  ems.writeModelToFile(options, "adlittle.ems", lp_mps);

  // Read ems and compare.
  options.model_file = "adlittle.ems";  // todo: check how to specify path

  HighsLp lp_ems;
  HighsStatus ems_read_status = loadLpFromFile(options, lp_ems);
  REQUIRE(ems_read_status == HighsStatus::OK);

  bool are_the_same = lp_mps == lp_ems;
  REQUIRE(are_the_same);

  std::remove(options.model_file.c_str());
}

TEST_CASE("integrality-constraints", "[highs_filereader]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/small_mip.mps";

  HighsOptions options;
  options.model_file = filename;
  // integer variables are COL03,COL04 so x[2], x[3].
  const std::vector<int> kIntegers{0, 0, 1, 1, 0, 0, 0, 0};

  // Read mps with fixed format parser.
  HighsLp lp_fixed;
  options.mps_parser_type_free = false;

  HighsStatus read_status = loadLpFromFile(options, lp_fixed);
  REQUIRE(read_status == HighsStatus::OK);
  REQUIRE(lp_fixed.integrality_.size() == lp_fixed.numCol_);
  REQUIRE(lp_fixed.integrality_ == kIntegers);

  // Read mps with free format parser.
  HighsLp lp_free;
  options.mps_parser_type_free = true;

  read_status = loadLpFromFile(options, lp_free);
  REQUIRE(read_status == HighsStatus::OK);
  REQUIRE(lp_free.integrality_.size() == lp_free.numCol_);
  REQUIRE(lp_free.integrality_ == kIntegers);
}

TEST_CASE("dualize", "[highs_data]") {
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
