# src/CMakeLists.txt

# Enable MPI support based on user option
option(USE_MPI "Use MPI for parallel compilation" ON)

# Find MPI package if USE_MPI is ON
if(USE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "Using MPI: ${MPI_C_COMPILER}")
else()
    message(STATUS "Building without MPI")
endif()


# Add the executable
add_executable(simulation
    main.cpp
    Particle.cpp
    Field.cpp
    Equilibrium.cpp
    ../inih/INIReader.cpp 
    ../inih/ini.c
)

# Include headers (including MPI if enabled)
target_include_directories(simulation PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_INCLUDE_PATH}  # Only add MPI include path if MPI is enabled
)

# Output the executable in the bin directory
set_target_properties(simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Link MPI libraries if USE_MPI is ON
if(USE_MPI)
    target_link_libraries(simulation MPI::MPI_CXX)  # Using modern CMake MPI target
else()
    target_link_libraries(simulation)  # No MPI linkage for serial build
endif()