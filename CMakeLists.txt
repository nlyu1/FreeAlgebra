cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(FreeAlgebraProject)

# Specify C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID MATCHES "Clang")
  # Compiler options
  add_compile_options(-Wall -Wextra -pedantic)
  add_compile_options(-Wshadow -Wnon-virtual-dtor -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wmisleading-indentation -Wnull-dereference -Wdouble-promotion -Wformat=2)
  # add_compile_options(-O3)
  # add_compile_options(-fsanitize=address)
  # add_link_options(-fsanitize=address)
endif()


# Debug and testing features
enable_testing()
set(CMAKE_BUILD_TYPE Debug)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()


# Set the path to the extracted LibTorch or use find_package if installed globally
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/pytorch/2.2.0_3")
include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories("/opt/homebrew/Cellar/fmt/10.2.1_1/include")
include_directories("/opt/homebrew/Cellar/gperftools/2.15")

# Dependencies
find_package(Torch REQUIRED)
set(TORCH_LIBRARIES "torch"
                    "torch_library"
                    "/opt/homebrew/Cellar/pytorch/2.2.0_3/lib/libc10.dylib"
)
# message(STATUS "TORCH_LIBRARIES: ${TORCH_LIBRARIES}")
find_package(fmt REQUIRED)

# Define the INTERFACE library
add_library(FreeAlgebraLib INTERFACE)
# Specify include directories
target_include_directories(FreeAlgebraLib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(FreeAlgebraLib INTERFACE ${TORCH_LIBRARIES} fmt::fmt)

if(COMMAND target_precompile_headers)
  target_precompile_headers(FreeAlgebraLib INTERFACE
    <vector>
    <string>
    <sstream>
    <cassert>
    <cstring>
    <unordered_map>
    <Eigen/Dense>
    <iomanip>
    <iostream>
    <cmath>
    <complex>
    <Eigen/Dense>
    <torch/torch.h>
    <fmt/core.h>
  )
endif()

find_library(PROFILER_LIB NAMES profiler PATHS /usr/local/lib /opt/homebrew/lib)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE FreeAlgebraLib)
target_link_libraries(main PRIVATE ${PROFILER_LIB})



# add_executable(test_finitealgebras ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_finitealgebras.cpp)
# target_link_libraries(test_finitealgebras PRIVATE FreeAlgebraLib fmt::fmt ${TORCH_LIBRARIES})
# add_test(NAME test_finitealgebras COMMAND test_finitealgebras)

# add_executable(test_fermions ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_fermions.cpp)
# target_link_libraries(test_fermions PRIVATE FreeAlgebraLib fmt::fmt ${TORCH_LIBRARIES})
# add_test(NAME test_fermions COMMAND test_fermions)