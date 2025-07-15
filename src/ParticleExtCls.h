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

  // --------------------
  void particle_ext_cls_test(Equilibrium &equ, FieldCls &fd, int irk, int nrun,
                             double dt_o_Ttr, int iset_track);

  void particle_ext_cls_init(const Equilibrium &equ,
                             std::optional<int> sp_id = std::nullopt);
  void particle_ext_cls_init(const Equilibrium &equ, int rank, int mpisize_in);
  //
  void particle_ext_cls_dxvpardt123EM2d1f_2sp(
      Equilibrium &equ, FieldCls &fd,
      const std::vector<std::complex<double>> &phik,
      const std::vector<std::complex<double>> &apark,
      const std::vector<int> &ntor1d, 
      const std::vector<std::complex<double>> &amp,
      std::vector<ParticleCoords> &dxvdt);

  // ---- Time-stepping (Electromagnetic) ----
  void particle_ext_cls_dxvpardt123EM2d1f(
      int speciesIndex, Equilibrium &equ, FieldCls &fd,
      const std::vector<std::complex<double>> &phik,
      const std::vector<std::complex<double>> &apark,
      const std::vector<int> &ntor1d, 
      const std::vector<std::complex<double>> &amp,
      const std::vector<double> &partrad0,
      const std::vector<double> &parttheta0,
      const std::vector<double> &partphitor0,
      const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
      const std::vector<double> &partw0, const std::vector<double> &partfog0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      std::vector<double> &dwdt);
  void particle_ext_cls_dxvpardt123EMgeneral(
      int speciesIndex, Equilibrium &equ, FieldCls &fd,
      const std::vector<double> &partrad0,
      const std::vector<double> &parttheta0,
      const std::vector<double> &partphitor0,
      const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
      const std::vector<double> &partw0, const std::vector<double> &partfog0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      std::vector<double> &dwdt, int icase,
      const std::vector<std::complex<double>> &phik,
      const std::vector<std::complex<double>> &apark,
      const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp);
};

#endif // PARTICLEEXTCLS_H
