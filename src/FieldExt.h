#ifndef FIELDEXTCLS_H
#define FIELDEXTCLS_H

#include "Equilibrium.h"
#include "Field.h" // Base class
#include "MPIManager.h"
#include "Particle.h"
#include <complex>
#include <vector>

class FieldExtCls : public FieldCls {
private:
  int rank, size;

public:
  // Fields & moments arrays
  int ntotfem2d1f, ntotdof2d1f;
  std::vector<int> idxdof2d1f;

  std::vector<std::complex<double>> denskMkj, phik, jparkMkj, apark;
  // std::vector<std::complex<double>> aparsk, aparhk;
  // int imixvar = 0; //0, no mixvar, 1 mixing variable

  // Constructor & Destructor
  FieldExtCls() {
    std::cout << "FieldExtCls default constructor called" << std::endl;
  };
  ~FieldExtCls() = default;

  // Methods
  void init(const Equilibrium &equ, Particle &pt);
  void initializePerturbations(const Equilibrium &equ);
  void restart(int flag){};
  void field_ext_cls_test(const Equilibrium &equ, ParticleSpecies &pt);
  void field_ext_cls_calc_W(std::vector<double> &WWW, const Equilibrium &equ,
                            Particle &pt,
                            const std::vector<std::complex<double>> &phik_c,
                            const std::vector<int> &ntor1d);
  void field_ext_cls_calc_T_oneSpecies(
      std::vector<std::complex<double>> &TTT, const Equilibrium &equ,
      ParticleSpecies &species, const std::vector<double> &partrad0,
      const std::vector<double> &parttheta0,
      const std::vector<double> &partphitor0,
      const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
      const std::vector<double> &partw0, const std::vector<double> &partfog0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      std::vector<double> &dwdt, int icase,
      const std::vector<std::complex<double>> &phik_c,
      const std::vector<int> &ntor1d);
};
#endif // FIELDEXTCLS_H

