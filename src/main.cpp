/*
 * Author: Guo Meng
 * Email: guo.meng@ipp.mpg.de
 * Created Date: 2024-12-10
 * Last Modified: 2025-03-10
 * License: MIT License
 *
 * Description:
 * This file is part of MENG++ project.
 */
#include <chrono>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <stdexcept>

#include "Equilibrium.h"
#include "Field.h"
#include "FieldExt.h"
#include "GKEM2D1FCls.h"
#include "MPIManager.h"
#include "Particle.h"
#include "ParticleExtCls.h"
#include "SplineCubicNd.h"
#include "util_io.h"

class TestCls {
private:
  MPIManager &mpiManager;
  int rank, size;

public:
  TestCls(int argc, char **argv)
      : mpiManager(MPIManager::getInstance(argc, argv)),
        rank(mpiManager.getRank()), size(mpiManager.getSize()) {}

  void run() {
    std::cout << "Hello, World! I am process " << rank << " of " << size
              << ".\n";

    // Create a particle species
    ParticleSpecies electrons(2, 1.0, -1.0);

    // Create an equilibrium state
    Equilibrium eq;
    eq.readInput("input.ini");
    eq.printState();
    std::cout << "qbar(0) = " << eq.get_qbar(0.0) << std::endl;
    std::cout << "eq.rmaxis = " << eq.rmaxis << std::endl;
    std::cout << "cos theta = " << eq.getcos_straight_adhoc(1.0, 0.0, eq.rmaxis)
              << std::endl;
    std::cout << "R(1,0) = " << eq.getR(1.0, 0.0) << std::endl;

    // Create particle object for test
    Particle particle(eq, rank, size);

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
    if (rank == 0)
      std::cout << "========Init,Test Field in Testparticle========"
                << std::endl;
    FieldCls fd;
    fd.field_cls_init(equ);

    // Particle part
    if (rank == 0)
      std::cout << "========Init,Test Particle in Testparticle========"
                << std::endl;
    ParticleExtCls pt_ext(equ, pt, rank, size);

    if (icase == 0) {
      // Test Particle Check
      if (rank == 0) {
        std::cout << "========TEST Particle Class w/o field========"
                  << std::endl;
      }
      pt.particle_cls_test(equ, irk, nrun, dt_o_Ttr, iset_track);
    } else if (icase == 1) {
      // Particle Test with Field (commented out)
      if (rank == 0)
        std::cout << "========TEST Particle Class with field========"
                  << std::endl;
      std::cout << "========Only test whether it can compile and run.========"
                << std::endl;
      std::cout << "========Init,Test FieldExt========" << std::endl;
      pt_ext.particle_ext_cls_test(equ, fd, irk, nrun, dt_o_Ttr, iset_track);
    }

    t2 = MPI_Wtime();
    if (rank == 0) {
      std::cout << "---- cputime: " << t2 - t1 << " ----" << std::endl;
    }

    if (rank == 0) {
      std::cout << "====Test Particle Ends====" << std::endl;
    }
  }

  void testSplineNd() {

    if (rank == 0)
      std::cout << "========== Test SplineNd Starts ==========" << std::endl;
    SplineCubicNdCls spc1d, f_bsp;

    int isource = 1;

    int fic, fjc, fkc;
    int idiffx = 0, idiffy = 0, idiffz = 0;
    int nqx = 40, nqy = 41, nqz = 42;
    std::vector<double> xq1d(nqx), yq1d(nqy), zq1d(nqz), fq1d;
    double t0, t1, t2;

    // Initialize the spline (Assume SplineCubicNd has testilar methods)
    spc1d.spc_cls_init(isource);
    spc1d.spc_cls_test();

    int ndim = spc1d.getDimension();
    if (ndim >= 1) {
      for (int i = 0; i < nqx; i++)
        xq1d[i] = spc1d.getZMin(0) - spc1d.getDz(0) +
                  (spc1d.getZWidth(0) + 2 * spc1d.getDz(0)) * i / (nqx - 1);
    }

    if (ndim >= 2) {
      for (int i = 0; i < nqy; i++)
        yq1d[i] = spc1d.getZMin(1) - spc1d.getDz(1) +
                  (spc1d.getZWidth(1) + 2 * spc1d.getDz(1)) * i / (nqy - 1);
    }

    if (ndim >= 3) {
      for (int i = 0; i < nqz; i++)
        zq1d[i] = spc1d.getZMin(2) - spc1d.getDz(2) +
                  (spc1d.getZWidth(2) + 2 * spc1d.getDz(2)) * i / (nqz - 1);
    }

    // Allocate for function values
    if (ndim == 1)
      fq1d.resize(nqx);
    else if (ndim == 2)
      fq1d.resize(nqx * nqy);
    else if (ndim == 3)
      fq1d.resize(nqx * nqy * nqz);

    // !calc. interpolated values
    // Perform interpolation (timing included)
    auto start = std::chrono::high_resolution_clock::now();
    t1 = MPI_Wtime();
    if (ndim == 1) {
      // !=========================1D interpolation=======================
      if (rank == 0) {
        std::cout
            << "--------f_bsp 1d test starts: initialize, evaluate--------"
            << std::endl;
      }
      // call f_bsp%initialize1d(spc1d%z1d_arr(1)%v,
      // spc1d%fval(1:spc1d%nnode_arr(1)), spc1d%bc_arr(1), spc1d%ext_arr(1))
      f_bsp.initialize1d(spc1d.getZ1d(0), spc1d.getFval(), spc1d.getBC(0),
                         spc1d.getExt(0));
      for (int i = 0; i < nqx; i++)
        f_bsp.evaluate1d(xq1d[i], idiffx, fq1d[i]);
    } else if (ndim == 2) {
      // !=========================3D interpolation=======================
      if (rank == 0) {
        std::cout
            << "--------f_bsp 2d test starts: initialize, evaluate--------"
            << std::endl;
      }
      f_bsp.initialize2d(spc1d.getZ1d(0), spc1d.getZ1d(1), spc1d.getFval(),
                         spc1d.getBC(0), spc1d.getBC(1), spc1d.getExt(0),
                         spc1d.getExt(1));
      std::cout << "Finished initialize 2d" << std::endl;
      for (int i = 0; i < nqx; i++)
        for (int j = 0; j < nqy; j++)
          f_bsp.evaluate2d(xq1d[i], yq1d[j], idiffx, idiffy, fq1d[i + j * nqx]);
    } else if (ndim == 3) {
      // !=========================3D interpolation=======================
      if (rank == 0) {
        std::cout
            << "--------f_bsp 3d test starts: initialize, evaluate--------"
            << std::endl;
      }
      f_bsp.initialize3d(spc1d.getZ1d(0), spc1d.getZ1d(1), spc1d.getZ1d(2),
                         spc1d.getFval(), spc1d.getBC(0), spc1d.getBC(1),
                         spc1d.getBC(2), spc1d.getExt(0), spc1d.getExt(1),
                         spc1d.getExt(2));
      for (int i = 0; i < nqx; i++)
        for (int j = 0; j < nqy; j++)
          for (int k = 0; k < nqz; k++)
            f_bsp.evaluate3d(xq1d[i], yq1d[j], zq1d[k], idiffx, idiffy, idiffz,
                             fq1d[i + j * nqx + k * nqx * nqy]);
    }
    auto end = std::chrono::high_resolution_clock::now();
    t2 = MPI_Wtime();

    if (rank == 0) {
      std::chrono::duration<double> elapsed = end - start;
      std::cout << "Interpolation Time: " << elapsed.count() << " seconds"
                << std::endl;
      std::cout << "---- cputime: " << t2 - t1 << " ----" << std::endl;
    }

    if (rank == 0) {
      UtilIO::write_arr1d_r8(fq1d, "data_fval_intp_dense.txt", true);
      // UtilIO::write_arr1d_r8(f_bsp.fspl, "data_fspl_oop.txt");
      // UtilIO::write_arr1d_r8(f_bsp.fval, "data_fval_intp.txt");
    }

    if (rank == 0)
      std::cout << "========== Test SplineNd Ends ==========" << std::endl;
  }

  int testField() {
    if (rank == 0) {
      std::cout << "========== Test Field Starts ==========" << std::endl;
    }

    Equilibrium equ;
    equ.readInput("input.ini");
    // FieldCls field;
    // field.field_cls_init(equ);
    std::vector<ParticleSpecies> particleList; // Create the actual vector
    // std::vector<ParticleSpecies> &pt = particleList; // Reference to it

    // Create a particle species
    int nptot_all = 10;
    double mass = 1.0;
    double zcharge = -1.0;
    ParticleSpecies pt(nptot_all, mass, zcharge);
    pt.setSpId(0);

    FieldExtCls fd_ext;
    FieldCls fd;
    std::cout << "========Init,Test FieldExt========" << std::endl;
    fd.field_cls_init(equ);
    // Create particle object for test
    Particle particle(equ, rank, size);

    fd_ext.init(equ, particle);
    fd_ext.field_ext_cls_test(equ, pt);
    std::vector<double> WWW;
    WWW.assign(fd_ext.ntor1d.size(), 0.0); //

    fd_ext.field_ext_cls_calc_W(WWW, equ, particle, fd_ext.phik, fd_ext.ntor1d);
    // Print the calculated WWW
    std::cout << "Calculated WWW: [";
    for (size_t i = 0; i < WWW.size(); ++i) {
      std::cout << WWW[i];
      if (i != WWW.size() - 1)
        std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    // try {
    //     field.readInput("input.ini");
    // } catch (const std::exception &e) {
    //     std::cerr << e.what() << std::endl;
    //     return EXIT_FAILURE;  // Return failure status
    // }

    return EXIT_SUCCESS; // Return success status
  }
};

int main(int argc, char **argv) {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  // TestCls test(argc, argv);
  // // test.run();
  // // test.testParticle();
  // // test.testSplineNd();
  // test.testField();
  // return 0;

  GKEM2D1FCls gkem2d1f(argc, argv);
  gkem2d1f.test(); // test with itest = 0 1 2
  
  auto end_time = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds = end_time - start_time;
  std::cout << "Execution time: " << elapsed_seconds.count() << " seconds" << std::endl;
}
