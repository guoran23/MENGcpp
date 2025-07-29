#ifndef PARTICLEEXTCLS_H
#define PARTICLEEXTCLS_H

#include "Equilibrium.h"
#include "Field.h"
#include "MPIManager.h"
#include "Particle.h"
#include "util_math.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// =====================================================
// Class ParticleExtCls extends ParticleCls
// Fortran: type, extends(particle_cls)::particle_ext_cls
// =====================================================
class ParticleExtCls : public Particle {
public:
  bool batch_g2p = false;

  // 默认构造函数（如果需要）
  ParticleExtCls() : Particle() {
    std::cout << "ParticleExtCls default constructor called" << std::endl;
  }

  // 自定义构造函数
  ParticleExtCls(Equilibrium &equ_in, int rank_in, int mpisize_in)
      : Particle(equ_in, rank_in, mpisize_in) // 初始化基类
  {
    std::cout << "ParticleExtCls custom constructor called" << std::endl;
  }
  ParticleExtCls(Equilibrium &equ_in, ParticleSpecies &sp_in, int rank_in,
                 int mpisize_in)
      : Particle(equ_in, sp_in, rank_in, mpisize_in) // 初始化基类
  {
    std::cout << "ParticleExtCls custom constructor called" << std::endl;
  }

  void add_vectors(std::vector<std::complex<double>> &a,
                   const std::vector<std::complex<double>> &b) {
    if (a.size() != b.size()) {
      throw std::runtime_error("add_vectors: size mismatch between vectors");
    }

    for (size_t i = 0; i < a.size(); ++i) {
      a[i] += b[i];
    }
  }

  void print_complex_vector(const std::string &name,
                         const std::vector<std::complex<double>> &omega_1_tmp) {
    
      std::cout << name << " = ";
      for (const auto &val : omega_1_tmp) {
        std::cout << std::scientific << std::setprecision(15) << "("
                  << val.real() << "," << val.imag() << " i ) ";
      }
      std::cout << std::endl;
    
  }

  // --------------------
  void particle_ext_cls_test(Equilibrium &equ, FieldCls &fd, int irk, int nrun,
                             double dt_o_Ttr, int iset_track);

  void particle_ext_cls_init(const Equilibrium &equ,
                             std::optional<int> sp_id = std::nullopt);
  void particle_ext_cls_init(const Equilibrium &equ, int rank, int mpisize_in);
  //
  void particle_ext_cls_dxvpardt123EM2d1f_2sp(
      const Equilibrium &equ, const FieldCls &fd,
      const std::vector<std::complex<double>> &phik,
      const std::vector<std::complex<double>> &apark,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amp,
      std::vector<ParticleCoords> &dxvdt,
      std::vector<std::complex<double>> &TTT_allsp);

  // ---- Time-stepping (Electromagnetic) ----
  void particle_ext_cls_dxvpardt123EM2d1f(
      const int speciesIndex, const Equilibrium &equ, const FieldCls &fd,
      const std::vector<std::complex<double>> &phik,
      const std::vector<std::complex<double>> &apark,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amp,
      std::vector<double> &partrad0,
      std::vector<double> &parttheta0,
      std::vector<double> &partphitor0,
      std::vector<double> &partvpar0, std::vector<double> &partmu0,
      std::vector<double> &partw0, std::vector<double> &partfog0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      std::vector<double> &dwdt, std::vector<std::complex<double>> &TTT_onesp);
  void particle_ext_cls_dxvpardt123EMgeneral(
      const int speciesIndex, const Equilibrium &equ, const FieldCls &fd,
      std::vector<double> &partrad0,
      std::vector<double> &parttheta0,
      std::vector<double> &partphitor0,
      std::vector<double> &partvpar0, std::vector<double> &partmu0,
      std::vector<double> &partw0, std::vector<double> &partfog0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      std::vector<double> &dwdt, int icase,
      const std::vector<std::complex<double>> &phik,
      const std::vector<std::complex<double>> &apark,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amp,
      std::vector<std::complex<double>> &TTT_onesp);
  std::vector<std::complex<double>>
  calc_T_onePar(const Equilibrium &equ, const FieldCls &fd,
                const double &partrad0, const double &parttheta0,
                const double &partphitor0, const double &partvpar0,
                const double &partw0, const double &vd_rad,
                const double &vd_the, const double &vd_phi,
                const std::vector<std::complex<double>> &phik_c,
                const std::vector<int> &ntor1d,
                const std::vector<std::complex<double>> &amp);
};

#endif // PARTICLEEXTCLS_H
