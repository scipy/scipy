# create highs binary using library without pic
add_executable(highs RunHighs.cpp ${HIGHS_SOURCE_DIR}/src/lp_data/HighsRuntimeOptions.h)
target_link_libraries(highs libhighs)

# install the binary
install(TARGETS highs EXPORT highs-targets
        RUNTIME DESTINATION bin)
