cmake_minimum_required(VERSION 3.12)
project(FreeAlgebraProject)

# Specify C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID MATCHES "Clang")
  # Add the most pedantic warnings and treat all warnings as errors.
  add_compile_options(-Wall -Wextra -pedantic)
  
  # Additional warnings that are not included in -Wall or -Wextra
  add_compile_options(-Wshadow -Wnon-virtual-dtor -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wmisleading-indentation -Wnull-dereference -Wdouble-promotion -Wformat=2)
  
  # Enable AddressSanitizer (a fast memory error detector)
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()



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


# Add an executable for tests
add_executable(test_algebra test.cpp) # Add this line, assuming test.cpp is your test file
target_link_libraries(test_algebra PRIVATE fmt::fmt)