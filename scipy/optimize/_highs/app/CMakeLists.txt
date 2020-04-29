# create highs binary using library without pic
add_executable(highs)
target_sources(highs PRIVATE RunHighs.cpp)
target_include_directories(highs PRIVATE ${HIGHS_SOURCE_DIR}/app)
target_link_libraries(highs libhighs)

# install the binary
install(TARGETS highs EXPORT highs-targets
        RUNTIME DESTINATION bin)
