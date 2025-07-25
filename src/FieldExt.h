#ifndef FIELDEXTCLS_H
#define FIELDEXTCLS_H

#include "Equilibrium.h"
#include "Field.h" // Base class
#include "MPIManager.h"
#include "Particle.h"
#include <complex>
#include <mpi.h>
#include <vector>

class FieldExtCls : public FieldCls {
private:
  int rank, size;

public:
  // int istart_shm, iend_shm;
  // MPI_Comm comm_shm;
  // Fields & moments arrays
  int ntotdof2d1f;
  std::vector<int> idxdof2d1f;

  std::vector<std::complex<double>> denskMkj, phik, jparkMkj, apark;
  std::vector<std::complex<double>>
      amplitude_arr; // Amplitude array for perturbations
  std::vector<double> omega0;
  // int imixvar = 0; //0, no mixvar, 1 mixing variable

  // Constructor & Destructor
  FieldExtCls() {
    std::cout << "FieldExtCls default constructor called" << std::endl;
  }
  
  ~FieldExtCls() = default;

  // Methods
  void init(const Equilibrium &equ, Particle &pt);
  void initializePerturbations(const Equilibrium &equ);
  void restart(int flag){};
  void field_ext_cls_test(const Equilibrium &equ, ParticleSpecies &pt);
  void
  field_ext_cls_calc_W(std::vector<double> &WWW, const Equilibrium &equ,
                       Particle &pt,
                       const std::vector<std::complex<double>> &phik_c,
                       const std::vector<int> &ntor1d);
                       
  void field_ext_cls_record1d(int istep) {
    // timer.timer_cls_tic(4);

    double Etot_phi = 0.0;

    for (int fic = 0; fic <= ntotfem2d1f; ++fic) {
      Etot_phi -= std::real(phik[fic] * denskMkj[fic]);
    }

    if (rank == 0) {
      std::string sfile = "data_energy.txt";
      std::ofstream file;

      if (istep == 0 && !irestart) {
        file.open(sfile, std::ios::out); // overwrite
      } else {
        file.open(sfile, std::ios::app); // append
      }

      if (file.is_open()) {
        file << std::scientific << std::setprecision(5) << Etot_phi << "\n";
        file.close();
      }
    }

    // timer.timer_cls_toc(4);
  }
};
#endif // FIELDEXTCLS_H
