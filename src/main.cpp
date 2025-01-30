#include "Equilibrium.h"
#include "Field.h"
#include "MPIManager.h"
#include "Particle.h"
#include <fstream>
#include <iostream>
#include <mpi.h>
// #include <petsc.h>

class Simulation {
public:
  Simulation(int argc, char **argv)
      : mpiManager(MPIManager::getInstance(argc, argv)),
        rank(mpiManager.getRank()), size(mpiManager.getSize()) {}

  void run() {
    std::cout << "Hello, World! I am process " << rank << " of " << size
              << ".\n";

    // Create a magnetic field
    Field field(10.0);

    // Create a particle species
    ParticleSpecies electrons(2, 1.0, -1.0);

    // Create an equilibrium state
    Equilibrium eq(101325.0, 300.0);
    eq.readInput("input.ini");

    eq.printState();
    std::cout << "qbar(0) = " << eq.get_qbar(0.0) << std::endl;
    std::cout << "eq.rmaxis = " << eq.rmaxis << std::endl;
    std::cout << "cos theta = " << eq.getcos_straight_adhoc(1.0, 0.0, eq.rmaxis)
              << std::endl;
    std::cout << "R(1,0) = " << eq.getR(1.0, 0.0) << std::endl;

    // Create particle object for simulation
    Particle particle(eq, rank, size);

    double timeStep = 0.1;
    for (int i = 0; i < 10; ++i) {
      // Simulate particle movement
      // particle.move(field, timeStep);
    }

    std::cout << "\nFinal State:\n";
    // particle.displayState(); // Uncomment if needed
    eq.displayEquilibrium();
  }

  void testParticle() {

    Equilibrium equ;
    // Create a particle species
    int nptot_all = 10;
    double mass = 1.0;
    double zcharge = -1.0;
    ParticleSpecies pt(nptot_all, mass, zcharge);
    pt.setSpId(0);

    int fic;
    double t1, t2;

    // INPUT from file
    int icase, irk, nrun, iset_track;
    double dt_o_Ttr;

    // Start Main program
    // Set Run Parameters
    std::ifstream input_file("input.particle");
    if (input_file.is_open()) {
      std::string line;
      while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '=')) {
          if (key == "icase") {
            iss >> icase;
          } else if (key == "irk") {
            iss >> irk;
          } else if (key == "nrun") {
            iss >> nrun;
          } else if (key == "dt_o_Ttr") {
            iss >> dt_o_Ttr;
          } else if (key == "iset_track") {
            iss >> iset_track;
          }
        }
      }
      input_file.close();
    } else {
      if (rank == 0)
        std::cerr << "Error opening input file!" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
      std::cout << "icase: " << icase << ", irk: " << irk << ", nrun: " << nrun
                << ", dt_o_Ttr: " << dt_o_Ttr << ", iset_track: " << iset_track
                << std::endl;
    }

    if (rank == 0)
      std::cout << "====Test Particle Starts====" << std::endl;
    t1 = MPI_Wtime();

    // Equilibrium part
    if (rank == 0)
      std::cout << "========Init Equilibrium in Testparticle========"
                << std::endl;
    // Create an equilibrium state
    equ.readInput("input.ini"); // instance the equ
    equ.printState();
    equ.writeWallRZ();

    // Field part (commented out)
    // if (rank == 0)
    //   std::cout << "========Init,Test Field in Testparticle========"
    //             << std::endl;
    // fd.field_cls_init(equ);

    // Particle part
    if (rank == 0)
      std::cout << "========Init,Test Particle in Testparticle========"
                << std::endl;
    // pt.particle_ext_cls_init(equ);

    if (icase == 0) {
      // Test Particle Check
      if (rank == 0) {
        std::cout << "========TEST Particle Class w/o field========"
                  << std::endl;
      }
      pt.particle_cls_test(equ, irk, nrun, dt_o_Ttr, iset_track);
    } else if (icase == 1) {
      // Particle Test with Field (commented out)
      // if (rank == 0)
      //     std::cout << "========TEST Particle Class with field========"
      //               << std::endl;
      // pt.particle_ext_cls_test(equ, fd, irk, nrun, dt_o_Ttr, iset_track);
    }

    t2 = MPI_Wtime();
    if (rank == 0) {
      std::cout << "---- cputime: " << t2 - t1 << " ----" << std::endl;
    }

    if (rank == 0)
      std::cout << "====Test Particle Ends====" << std::endl;
  }

private:
  MPIManager &mpiManager;
  int rank, size;
};

int main(int argc, char **argv) {
  Simulation sim(argc, argv);
  // sim.run();
  sim.testParticle();
  return 0;
}
