cmake_minimum_required(VERSION 3.20)

project(cgkmc_cpp LANGUAGES CXX)

add_library(cgkmc_core
    src/config.cpp
    src/dump_writer.cpp
    src/simulation.cpp
)
target_include_directories(cgkmc_core PUBLIC include)
target_compile_features(cgkmc_core PUBLIC cxx_std_20)

add_executable(cgkmc_run src/main.cpp)
target_link_libraries(cgkmc_run PRIVATE cgkmc_core)

option(CGKMC_BUILD_TESTS "Build cgkmc C++ tests" ON)
if(CGKMC_BUILD_TESTS)
    enable_testing()
    add_executable(cgkmc_core_tests tests/test_core.cpp)
    target_link_libraries(cgkmc_core_tests PRIVATE cgkmc_core)
    add_test(NAME cgkmc_core_tests COMMAND cgkmc_core_tests)
endif()

option(CGKMC_BUILD_BENCHMARKS "Build cgkmc C++ smoke benchmarks" ON)
if(CGKMC_BUILD_BENCHMARKS)
    add_executable(cgkmc_petn_smoke benchmarks/petn_smoke.cpp)
    target_link_libraries(cgkmc_petn_smoke PRIVATE cgkmc_core)
endif()
