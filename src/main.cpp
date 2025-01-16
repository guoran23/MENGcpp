#include <mpi.h>
#include <iostream>
#include "Field.h"
#include "Particle.h"
#include "Equilibrium.h"

class MPIManager {
public:
    MPIManager(int argc, char** argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }

    ~MPIManager() { MPI_Finalize(); }

    int getRank() const { return rank; }
    int getSize() const { return size; }

private:
    int rank, size;
};

class Simulation {
public:
    Simulation(int argc, char** argv) : mpiManager(argc, argv), rank(mpiManager.getRank()), size(mpiManager.getSize()) {}

    void run() {
        std::cout << "Hello, World! I am process " << rank << " of " << size << ".\n";
        
        // Create a magnetic field
        Field field(10.0);

        // Create a particle species
        ParticleSpecies electrons(2, 1.0, -1.0);
        
        // Create an equilibrium state
        Equilibrium eq(101325.0, 300.0);
        eq.readInput("input.ini");
        
        eq.printState();
        std::cout << "qbar(0) = " << eq.get_qbar(0.0) << std::endl;

        // Create particle object for simulation
        Particle particle(eq, rank, size);

        double timeStep = 0.1;
        for (int i = 0; i < 10; ++i) {
            // Simulate particle movement (commented out for now)
            // particle.move(field, timeStep);
        }

        std::cout << "\nFinal State:\n";
        // particle.displayState(); // Uncomment if needed
        eq.displayEquilibrium();
    }

private:
    MPIManager mpiManager;
    int rank, size;
};

int main(int argc, char** argv) {
    Simulation sim(argc, argv);
    sim.run();
    return 0;
}
