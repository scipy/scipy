/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file TestOsi.cpp
 * @brief Osi/HiGHS unit test
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include <iostream>

#include "CoinHelperFunctions.hpp"
#include "CoinPragma.hpp"
#include "HCheckConfig.h"
#include "OsiHiGHSSolverInterface.hpp"
#include "OsiUnitTests.hpp"

using namespace OsiUnitTest;

int main(int argc, const char* argv[]) {
#ifndef COINSAMPLEFOUND
  std::cerr << "Path to Data/Sample not known. Cannot run tests without sample "
               "MPS files."
            << std::endl;
  return EXIT_FAILURE;
#endif

  // Synchronise C++ stream i/o with C stdio.
  std::ios::sync_with_stdio();
  setbuf(stderr, 0);
  setbuf(stdout, 0);

  // Suppress an popup window that Windows shows in response to a crash.
  WindowsErrorPopupBlocker();

  OsiUnitTest::verbosity = 10;
  // OsiUnitTest::haltonerror = 2;

  int nerrors;
  int nerrors_expected;

  std::string mpsDir = COINSAMPLEDIR;
  mpsDir.push_back(CoinFindDirSeparator());

  // Do common solverInterface testing by calling the base class testing method.
  {
    OsiHiGHSSolverInterface highsSi;
    OSIUNITTEST_CATCH_ERROR(
        OsiSolverInterfaceCommonUnitTest(&highsSi, mpsDir, ""), {}, highsSi,
        "osi common unittest");
  }

  // Test Osi{Row,Col}Cut routines
  {
    OsiHiGHSSolverInterface highsSi;
    testingMessage("Testing OsiRowCut with OsiHiGHSSolverInterface\n");
    OSIUNITTEST_CATCH_ERROR(OsiRowCutUnitTest(&highsSi, mpsDir), {}, highsSi,
                            "rowcut unittest");
  }
  {
    OsiHiGHSSolverInterface highsSi;
    testingMessage("Testing OsiColCut with OsiHiGHSSolverInterface\n");
    OSIUNITTEST_CATCH_ERROR(OsiColCutUnitTest(&highsSi, mpsDir), {}, highsSi,
                            "colcut unittest");
  }
  {
    OsiHiGHSSolverInterface highsSi;
    testingMessage("Testing OsiRowCutDebugger with OsiHiGHSSolverInterface\n");
    OSIUNITTEST_CATCH_ERROR(OsiRowCutDebuggerUnitTest(&highsSi, mpsDir), {},
                            highsSi, "rowcut debugger unittest");
  }

#ifdef COINNETLIBFOUND
  //  We have run the fast unit tests.
  //  If there were no errors, then also run the Netlib problems.
  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);
  if (nerrors <= nerrors_expected) {
    std::string netlibDir = COINNETLIBDIR;
    netlibDir.push_back(CoinFindDirSeparator());

    // Create vector of solver interfaces
    std::vector<OsiSolverInterface*> vecSi(1, new OsiHiGHSSolverInterface);

    testingMessage("Testing OsiSolverInterface on Netlib problems.\n");
    OSIUNITTEST_CATCH_ERROR(OsiSolverInterfaceMpsUnitTest(vecSi, netlibDir), {},
                            "highs", "netlib unittest");

    delete vecSi[0];
  } else {
    testingMessage(
        "Skip testing OsiSolverInterface on Netlib problems as there were "
        "unexpected errors in previous runs.\n");
  }
#else
  testingMessage(
      "Skip testing OsiSolverInterface on Netlib problems as path to "
      "Data/Netlib not known.\n");
#endif

  // We're done. Report on the results.
  std::cout.flush();
  outcomes.print();

  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);
  if (nerrors > nerrors_expected)
    std::cerr << "Tests completed with " << nerrors - nerrors_expected
              << " unexpected errors." << std::endl;
  else
    std::cerr << "All tests completed successfully\n";

  return nerrors - nerrors_expected;
}
