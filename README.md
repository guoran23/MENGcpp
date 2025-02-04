# MENG++ Code

## Multiscale Exploration of Nonlinear Gyrokinetics in C++

### Build and Compile

1. Load the necessary modules:
   ```bash
   module load cmake gcc openmpi
   ```
2. Configure the build using CMake:
   ```bash
   cmake --preset mpi
   ```
3. Compile the code:
   ```bash
   cmake --build ./build_MPI/
   ```
4. The executable file will be located at:
   ```plaintext
   ./build_MPI/bin/simulation
   ```
   and the input files are located at:
   ```plaintext
   ./inputfiles/
   ```

