cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(FreeAlgebraProject)

# Specify C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID MATCHES "Clang")
  # Compiler options
  add_compile_options(-Wall -Wextra -pedantic)
  add_compile_options(-Wshadow -Wnon-virtual-dtor -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wmisleading-indentation -Wnull-dereference -Wdouble-promotion -Wformat=2)
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

# Set the path to the extracted LibTorch or use find_package if installed globally
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/libtorch")
include_directories("/opt/homebrew/Cellar/fmt/10.2.1_1/include")

# Required libraries 
find_package(Torch REQUIRED)
find_package(fmt REQUIRED)


# Define the INTERFACE library
add_library(FreeAlgebraLib INTERFACE)
# Specify include directories
target_include_directories(FreeAlgebraLib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Dependencies
find_package(Torch REQUIRED)
find_package(fmt REQUIRED)
target_link_libraries(FreeAlgebraLib INTERFACE ${TORCH_LIBRARIES} fmt::fmt)


add_executable(main main.cpp)
target_link_libraries(main PRIVATE FreeAlgebraLib)

# set_target_properties(main PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
# )