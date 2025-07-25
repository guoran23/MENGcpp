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

# # to debug, only serial run
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fsanitize=address -fno-omit-frame-pointer")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")


# Add the executable
add_executable(simulation
    main.cpp
    Particle.cpp
    Field.cpp
    Equilibrium.cpp
    ../inih/INIReader.cpp 
    ../inih/ini.c
    util_math.cpp
    util_io.cpp
    SplineCubicNd.cpp
    SplineCubic2d1f.cpp
    FieldExt.cpp
    ParticleExtCls.cpp
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