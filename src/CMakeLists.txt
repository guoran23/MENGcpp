# src/CMakeLists.txt

# Add the executable
add_executable(simulation
    main.cpp
    Particle.cpp
    Field.cpp
    Equilibrium.cpp
    ../inih/INIReader.cpp 
    ../inih/ini.c
)

# Include headers
target_include_directories(simulation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Output the executable in the bin directory
set_target_properties(simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Link MPI libraries
target_link_libraries(simulation ${MPI_LIBRARIES})
