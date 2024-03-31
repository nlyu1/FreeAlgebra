cmake_minimum_required(VERSION 3.12)
project(PrettyComplexProject)

# Specify C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Libraries
include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories("/opt/homebrew/Cellar/fmt/10.2.1_1/include")

# Fetch pybind11
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.6.1  # Optionally, specify a version tag
)
FetchContent_MakeAvailable(pybind11)

# Define our module
pybind11_add_module(freealgebra MODULE freealgebra.cpp)

find_package(fmt REQUIRED) # Find fmt package
target_link_libraries(freealgebra PRIVATE fmt::fmt) # Link fmt library