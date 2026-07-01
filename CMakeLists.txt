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
