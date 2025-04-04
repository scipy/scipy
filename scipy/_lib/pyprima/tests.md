**Code must be battle-tested before becoming software.**
The development of PRIMA is driven by [intensive and extensive tests](https://github.com/orgs/libprima/discussions/39)
automated by [GitHub Actions](https://docs.github.com/en/actions).

Since each GitHub Team account can only run at most 60 GitHub Actions workflows concurrently, I have
to distribute this large amount of tests to multiple Team accounts as follows.

- [Tests](https://github.com/libprima/prima/actions) at [libprima/prima](https://github.com/libprima/prima)

    - [![CMake build](https://github.com/libprima/prima/actions/workflows/cmake.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake.yml)
    - [![CMake build, nagfor](https://github.com/libprima/prima/actions/workflows/cmake_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake_nagfor.yml)
    - [![CMake build, macOS ARM64](https://github.com/libprima/prima/actions/workflows/cmake_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake_mac.yml)
    - [![CMake build, macOS ARM64, nagfor](https://github.com/libprima/prima/actions/workflows/cmake_mac_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake_mac_nagfor.yml)
    - [![Test nagfor, macOS ARM64](https://github.com/libprima/prima/actions/workflows/test_nagfor_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/test_nagfor_mac.yml)
    - [![Test matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/test_matlab_mac.yml)
    - [![Stress test, matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/stress_test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/stress_test_matlab_mac.yml)
    - [![Parallel test, matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/parallel_test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/parallel_test_matlab_mac.yml)
    - [![Recursive test, matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/recursive_test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/recursive_test_matlab_mac.yml)
    - [![Test nagfor](https://github.com/libprima/prima/actions/workflows/test_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/test_nagfor.yml)
    - [![Build Python wheels](https://github.com/libprima/prima/actions/workflows/build_python.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/build_python.yml)
    - [![Compile MEX](https://github.com/libprima/prima/actions/workflows/compile_mex.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/compile_mex.yml)
    - [![Lint the Fortran code and the MEX gateways with nagfor](https://github.com/libprima/prima/actions/workflows/lint_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/lint_nagfor.yml)

- <a name="verification"></a> [Tests](https://github.com/libsprima/prima/actions) at [libsprima/prima](https://github.com/libsprima/prima)
    - [![Verification, small](https://github.com/libsprima/prima/actions/workflows/verify_small.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/verify_small.yml)
    - [![Verification, big](https://github.com/libsprima/prima/actions/workflows/verify_big.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/verify_big.yml)
    - [![Verification, large](https://github.com/libsprima/prima/actions/workflows/verify_large.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/verify_large.yml)
    - [![Plot performance profiles for lincoa](https://github.com/libsprima/prima/actions/workflows/profile_lincoa_small.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/profile_lincoa_small.yml)

- <a name="profiling"></a> [Tests](https://github.com/primapack/prima/actions) at [primapack/prima](https://github.com/primapack/prima)

    - [![Plot performance profiles for cobyla, small](https://github.com/primapack/prima/actions/workflows/profile_cobyla_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_cobyla_small.yml)
    - [![Plot performance profiles for uobyqa, small](https://github.com/primapack/prima/actions/workflows/profile_uobyqa_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_uobyqa_small.yml)
    - [![Plot performance profiles for newuoa, small](https://github.com/primapack/prima/actions/workflows/profile_newuoa_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_newuoa_small.yml)
    - [![Plot performance profiles for bobyqa, small](https://github.com/primapack/prima/actions/workflows/profile_bobyqa_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_bobyqa_small.yml)
    - [![Plot performance profiles for PRIMA, small](https://github.com/primapack/prima/actions/workflows/profile_prima_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_prima_small.yml)

- [Tests](https://github.com/fortlab/prima/actions) at [fortlab/prima](https://github.com/fortlab/prima)

    - [![Plot performance profiles for all problems, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_all_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_all_sq.yml)
    - [![Plot performance profiles for cobyla, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_cobyla_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_cobyla_small_sq.yml)
    - [![Plot performance profiles for uobyqa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_uobyqa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_uobyqa_small_sq.yml)
    - [![Plot performance profiles for newuoa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_newuoa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_newuoa_small_sq.yml)
    - [![Plot performance profiles for bobyqa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_bobyqa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_bobyqa_small_sq.yml)
    - [![Plot performance profiles for lincoa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_lincoa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_lincoa_small_sq.yml)

- [Tests](https://github.com/primalib/prima/actions) at [primalib/prima](https://github.com/primalib/prima)

    - [![Verification, archiva](https://github.com/primalib/prima/actions/workflows/verify_archiva.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/verify_archiva.yml)
    - [![Plot performance profiles for all problems](https://github.com/primalib/prima/actions/workflows/profile_all.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_all.yml)
    - [![Plot performance profiles, single](https://github.com/primalib/prima/actions/workflows/profile_single.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_single.yml)
    - [![Plot performance profiles, quadruple](https://github.com/primalib/prima/actions/workflows/profile_quadruple.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_quadruple.yml)
    - [![Plot performance profiles, int16 and int64](https://github.com/primalib/prima/actions/workflows/profile_integer_kind.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_integer_kind.yml)
    - [![Plot performance profiles, infnan](https://github.com/primalib/prima/actions/workflows/profile_infnan.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_infnan.yml)
    - [![Plot performance profiles, various compiler options](https://github.com/primalib/prima/actions/workflows/profile_compiler_options.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_compiler_options.yml)
    - [![Plot performance profiles, various npt](https://github.com/primalib/prima/actions/workflows/profile_npt.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_npt.yml)

- [Tests](https://github.com/opt4ai/prima/actions) at [opt4ai/prima](https://github.com/opt4ai/prima)

    - [![Test MATLAB](https://github.com/opt4ai/prima/actions/workflows/test_matlab.yml/badge.svg)](https://github.com/opt4ai/prima/actions/workflows/test_matlab.yml)
    - [![Test MATLAB, Linux](https://github.com/opt4ai/prima/actions/workflows/test_matlab_linux.yml/badge.svg)](https://github.com/opt4ai/prima/actions/workflows/test_matlab_linux.yml)

- [Tests](https://github.com/opt2ai/prima/actions) at [opt2ai/prima](https://github.com/opt2ai/prima)

    - [![Test MATLAB, macOS Intel](https://github.com/opt2ai/prima/actions/workflows/test_matlab_mac_intel.yml/badge.svg)](https://github.com/opt2ai/prima/actions/workflows/test_matlab_mac_intel.yml)
    - [![Test MATLAB, Windows](https://github.com/opt2ai/prima/actions/workflows/test_matlab_windows.yml/badge.svg)](https://github.com/opt2ai/prima/actions/workflows/test_matlab_windows.yml)

- [Tests](https://github.com/sprimalib/prima/actions) at [sprimalib/prima](https://github.com/sprimalib/prima)<a name="stress-tests"></a>

    - <a name="stress"></a>[![Stress test on large problems, MATLAB](https://github.com/sprimalib/prima/actions/workflows/stress_test_matlab.yml/badge.svg)](https://github.com/sprimalib/prima/actions/workflows/stress_test_matlab.yml)
    - [![Stress test on large problems, Fortran](https://github.com/sprimalib/prima/actions/workflows/stress_test_fortran.yml/badge.svg)](https://github.com/sprimalib/prima/actions/workflows/stress_test_fortran.yml)

- [Tests](https://github.com/zequipe/prima/actions) at [zequipe/prima](https://github.com/zequipe/prima)
    - [![Parallel test, MATLAB](https://github.com/zequipe/prima/actions/workflows/parallel_test_matlab.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/parallel_test_matlab.yml)
    - [![Recursive test, MATLAB](https://github.com/zequipe/prima/actions/workflows/recursive_test_matlab.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/recursive_test_matlab.yml)
    - [![Test Flang](https://github.com/zequipe/prima/actions/workflows/test_flang.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_flang.yml)
    - [![Test Flang in AMD AOCC](https://github.com/zequipe/prima/actions/workflows/test_aflang.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_aflang.yml)
    - [![Test nvfortran](https://github.com/zequipe/prima/actions/workflows/test_nvfortran.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_nvfortran.yml)
    - [![Test Oracle sunf95](https://github.com/zequipe/prima/actions/workflows/test_sunf95.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_sunf95.yml)
    - [![Test g95](https://github.com/zequipe/prima/actions/workflows/test_g95.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_g95.yml)
    - [![Test RESCUE and IDZ, classical](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_classical.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_classical.yml)
    - [![Test RESCUE and IDZ, modernized](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_modernized.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_modernized.yml)

- [Tests](https://github.com/equipez/prima/actions) at [equipez/prima](https://github.com/equipez/prima)
    - [![Test gfortran on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_gfortran_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_gfortran_kunpeng.yml)
    - [![Test Flang on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_flang_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_flang_kunpeng.yml)
    - [![Test nvfortran on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_nvfortran_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_nvfortran_kunpeng.yml)
    - [![Test armflang on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_armflang_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_armflang_kunpeng.yml)
    - [![CMake build on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/cmake_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/cmake_pi.yml)
    - [![Test gfortran on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_gfortran_pi64.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_gfortran_pi64.yml)
    - [![Test Flang on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_flang_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_flang_pi.yml)
    - [![Test armflang on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_armflang_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_armflang_pi.yml)
    - [![Test nvfortran on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_nvfortran_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_nvfortran_pi.yml)

- [Tests](https://github.com/s-prima/prima/actions) at [s-prima/prima](https://github.com/s-prima/prima)

    - [![Lint the Fortran code and the MEX gateways on GitHub hosted runners](https://github.com/s-prima/prima/actions/workflows/lint_hosted.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/lint_hosted.yml)
    - [![Test gfortran](https://github.com/s-prima/prima/actions/workflows/test_gfortran.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/test_gfortran.yml)
    - [![Test ifort](https://github.com/s-prima/prima/actions/workflows/test_ifort.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/test_ifort.yml)
    - [![Test ifx](https://github.com/s-prima/prima/actions/workflows/test_ifx.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/test_ifx.yml)
