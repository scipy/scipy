include(CTest)

# prepare Catch library
set(CATCH_INCLUDE_DIR ${HIGHS_SOURCE_DIR}/src/external/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/lp_data)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/ipx/src)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/basiclu/src)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/ipx/include)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/basiclu/include)

# Make test executable
set(TEST_SOURCES
    TestMain.cpp
    TestOptions.cpp
    TestIO.cpp
    TestSort.cpp
    TestSetup.cpp
    TestFilereader.cpp
    TestInfo.cpp
    TestBasisSolves.cpp
    TestLpValidation.cpp
    TestLpModification.cpp
    Avgas.cpp)

if (IPX_ON)
    set (TEST_SOURCES ${TEST_SOURCES} TestIpx.cpp)
endif()

add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests libhighs Catch)

if (OSI_FOUND AND BUILD_TESTING)
    pkg_check_modules(OSITEST osi-unittests)
    if (OSITEST_FOUND)
        include_directories(${HIGHS_SOURCE_DIR}/src)
	add_executable(osi_unit_tests TestOsi.cpp)
	target_link_libraries(osi_unit_tests OsiHighs Catch ${OSITEST_LIBRARIES} CoinUtils)
	target_include_directories(osi_unit_tests PUBLIC ${OSITEST_INCLUDE_DIRS} ${HIGHS_SOURCE_DIR}/src/interfaces)
	target_compile_options(osi_unit_tests PUBLIC ${OSITEST_CFLAGS_OTHER})
    endif (OSITEST_FOUND)
endif()

if(FORTRAN_FOUND)
   set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
   include_directories(${HIGHS_SOURCE_DIR}/src)
   add_executable(fortrantest TestFortranAPI.f90)
   target_link_libraries(fortrantest libhighs FortranHighs)
   target_include_directories(fortrantest PUBLIC ${HIGHS_SOURCE_DIR}/src/interfaces)
 else()
 endif(FORTRAN_FOUND)

# check the C API
add_executable(capi_unit_tests TestCAPI.c)
target_link_libraries(capi_unit_tests libhighs)
add_test(NAME capi_unit_tests COMMAND capi_unit_tests)

# Check whether test executable builds OK.
add_test(NAME unit-test-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target unit_tests
         )


# Avoid that several build jobs try to concurretly build.
set_tests_properties(unit-test-build
                     PROPERTIES
                     RESOURCE_LOCK unittestbin)

# create a binary running all the tests in the executable
add_test(NAME unit_tests_all COMMAND unit_tests --success)
set_tests_properties(unit_tests_all
                    PROPERTIES
                    DEPENDS unit-test-build)

if (OSITEST_FOUND)

add_test(NAME osi-unit-test-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target osi_unit_tests
         )

# Avoid that several build jobs try to concurretly build.
set_tests_properties(osi-unit-test-build
PROPERTIES
RESOURCE_LOCK osiunittestbin)

pkg_search_module(COINSAMPLE coindatasample)
if (COINSAMPLE_FOUND)
   pkg_get_variable(COINSAMPLEDIR coindatasample datadir)
endif ()

pkg_search_module(COINNETLIB coindatanetlib)
if (COINNETLIB_FOUND)
   pkg_get_variable(COINNETLIBDIR coindatanetlib datadir)
endif ()

configure_file(${HIGHS_SOURCE_DIR}/check/HCheckConfig.h.in ${HIGHS_BINARY_DIR}/HCheckConfig.h)

# create a binary running all the tests in the executable
add_test(NAME osi_unit_tests_all COMMAND osi_unit_tests)
set_tests_properties(osi_unit_tests_all
PROPERTIES
DEPENDS osi-unit-test-build)

endif()

if (GAMS_FOUND)
add_test(NAME rungams COMMAND ${HIGHS_SOURCE_DIR}/check/rungams.sh ${GAMS_ROOT})
endif (GAMS_FOUND)

# An individual test can be added with the command below but the approach
# above with a single add_test for all the unit tests automatically detects all
# TEST_CASEs in the source files specified in TEST_SOURCES. Do not define any
# tests in TestMain.cpp and do not define CATCH_CONFIG_MAIN anywhere else.
# add_test(NAME correct-print-test COMMAND unit_tests correct-print)

# --------------------------------------
# Another way of adding the tests. Needs a script from github repo and a
# Catch2 installation. So add tests manually if there is no build issues.
# catch_discover_tests(unit_test)

# --------------------------------------
# Run instance tests.
#
# define the set of feasible instances
set(successInstances
    "25fv47\;2888\; 5.501846\;"
    "80bau3b\;3760\; 9.872242\;"
    "adlittle\;81\; 2.254950\;"
    "afiro\;21\;-4.647531\;"
    "etamacro\;507\;-7.557152\;"
    "greenbea\;5249\;-7.255525\;"
    "shell\;624\; 1.208825\;"
    "stair\;555\;-2.512670\;"
    "standata\;74\; 1.257699\;"
    "standgub\;67\; 1.257699\;"
    "standmps\;215\; 1.406017\;"
    )

set(infeasibleInstances
    "bgetam\;        infeasible"
    "box1\;          infeasible"
    "ex72a\;         infeasible"
    "forest6\;       infeasible"
    "galenet\;       infeasible"
    "gams10am\;      infeasible"
    "klein1\;        infeasible"
#    "refinery\;      infeasible"
    "woodinfe\;      infeasible"
    )

set(unboundedInstances
     "gas11\;         unbounded"
    )

set(failInstances
    )

# define settings
set(settings
    "--presolve=off"
    "--presolve=on"
#    "--parallel=on"
    )

# define a macro to add tests
#
# add_instancetests takes an instance group and a status
# that the solver should report as arguments
macro(add_instancetests instances solutionstatus)
# loop over the instances
foreach(instance ${${instances}})
    # add default tests
    # treat the instance as a tuple (list) of two values
    list(GET instance 0 name)
    list(GET instance 1 iter)

    if(${solutionstatus} STREQUAL "Optimal")
        list(GET instance 2 optval)
    endif()

    # specify the instance and the settings load command
    set(inst "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps")

    # loop over all settings
    foreach(setting ${settings})
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting}
              ${inst})

        set_tests_properties (${name}${setting} PROPERTIES
                DEPENDS unit_tests_all)
        set_tests_properties (${name}${setting} PROPERTIES
                PASS_REGULAR_EXPRESSION
                "Model   status      : ${solutionstatus}")

        if(${solutionstatus} STREQUAL "Optimal")
            if(${setting} STREQUAL "--presolve=off")
                set_tests_properties (${name}${setting} PROPERTIES
                        PASS_REGULAR_EXPRESSION
                        "Simplex   iterations: ${iter}\nObjective value     : ${optval}")
            else()
                set_tests_properties (${name}${setting} PROPERTIES
                        PASS_REGULAR_EXPRESSION
                        "Objective value     : ${optval}")
            endif()
        endif()
    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

# add tests for success and fail instances
add_instancetests(successInstances "Optimal")
add_instancetests(failInstances "Fail")
add_instancetests(infeasibleInstances "Infeasible")
add_instancetests(unboundedInstances "Unbounded")

# todo: test objective value

# todo: test iteration count

# manually add presolve tests until we fix presolve
set(setting "--presolve=on")

set(inst "${HIGHS_SOURCE_DIR}/check/instances/adlittle.mps")
add_test(NAME "presolve-adlittle" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-adlittle" PROPERTIES PASS_REGULAR_EXPRESSION "Model   status      : Optimal")
set_tests_properties ("presolve-adlittle" PROPERTIES DEPENDS unit_tests_all)

set(inst "${HIGHS_SOURCE_DIR}/check/instances/scrs8.mps")
add_test(NAME "presolve-scrs8" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-scrs8" PROPERTIES PASS_REGULAR_EXPRESSION "Model   status      : Optimal")
set_tests_properties ("presolve-scrs8" PROPERTIES DEPENDS unit_tests_all)


set(inst "${HIGHS_SOURCE_DIR}/check/instances/woodinfe.mps")
add_test(NAME "presolve-woodinfe" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-woodinfe" PROPERTIES PASS_REGULAR_EXPRESSION "Model   status      : Infeasible")
set_tests_properties ("presolve-woodinfe" PROPERTIES DEPENDS unit_tests_all)


set(inst "${HIGHS_SOURCE_DIR}/check/instances/gas11.mps")
add_test(NAME "presolve-gas11" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-gas11" PROPERTIES PASS_REGULAR_EXPRESSION "Model   status      : Unbounded")
set_tests_properties ("presolve-gas11" PROPERTIES DEPENDS unit_tests_all)
