cmake_minimum_required(VERSION 3.10)
project(Simulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable MPI
find_package(MPI REQUIRED)

# Include directories for MPI
include_directories(${MPI_INCLUDE_PATH})

# Specify the source directory
add_subdirectory(src)

# Output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

