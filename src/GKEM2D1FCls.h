#ifndef GKEM2D1FCLS_H
#define GKEM2D1FCLS_H

#include "Equilibrium.h"
#include "FieldExt.h"
#include "MPIManager.h"
#include "Particle.h"
#include "ParticleExtCls.h"
#include <memory>

class GKEM2D1FCls {
private:
  MPIManager &mpiManager;

public:
  // 基本控制参数
  int *nstart = nullptr;
  int nsp, imixvar, ipullback;
  bool nl_pt_perturbed, nl_fd_solve;
  int nrun, ischeme;
  double dtoTA, dtoTN;
  int niter_fd, niter_pt;
  double time_all, time_save;
  int rank, size;
  
  // 物理对象
  Equilibrium equ;
  FieldExtCls fd;
  std::unique_ptr<ParticleExtCls> pt;

  // 工作空间
  std::vector<double> WWW;
  std::vector<std::complex<double>> TTT;
  std::vector<std::complex<double>> omega;
  std::vector<double> omega_0;
  std::vector<std::complex<double>> omega_1;
  std::vector<ParticleCoords> xv0, dxv_sum, dxvdt;
  std::vector<std::complex<double>> denskMkj_tot, jparkMkj_tot;
  std::vector<std::complex<double>> aparsk0, daparsk_sum, daparskdt;

  // 计时
  double t1 = 0.0, t2 = 0.0;

  // 构造与析构
  GKEM2D1FCls();
  GKEM2D1FCls(int argc, char **argv)
      : mpiManager(MPIManager::getInstance(argc, argv)),
        rank(mpiManager.getRank()), size(mpiManager.getSize()) {}

  // 成员函数
  void initialize();
  void finalize();
  void record(int istep){};
  void record2h5(int istep){};
  void onestep_rk4();
  void onestep_heun();
  void calc_T_allSpecies();
  void test();
};
#endif // GKEM2D1FCLS_H
void GKEM2D1FCls::initialize() {

  equ.readInput("input.ini");
  pt = std::make_unique<ParticleExtCls>(equ, rank, size); // 初始化粒子
  nsp = pt->getNsp(); // 获取粒子种类数
  double mass_ion = pt->mass_ion;

  nstart = new int(0); // 指针分配
  // 打印配置信息
  if (rank == 0) {
    printf(">> nsp=%d, imixvar=%d, ipullback=%d, nrun=%d\n", nsp, imixvar,
           ipullback, nrun);
    printf("   ischeme=%d, niter_fd=%d, niter_pt=%d\n", ischeme, niter_fd,
           niter_pt);
  }

  // 分配粒子
  // pt.resize(nsp);
  xv0.resize(nsp);
  dxvdt.resize(nsp);
  dxv_sum.resize(nsp);

  for (int fsc = 0; fsc < nsp; ++fsc) {
    ParticleSpecies &species = pt->group.getSpecies(fsc);
    int nptot = species.getNptot();
    xv0[fsc].initialize(nptot);
    dxvdt[fsc].initialize(nptot);
    dxv_sum[fsc].initialize(nptot);

    if (rank == 0) {
      printf("nptot * mpisize = %d\n", nptot * size);
    }
  }

  fd.init(equ, *pt);
  int lenntor = fd.lenntor;
  // for perturbations
  WWW.resize(lenntor, 0.0);
  std::complex<double> zero_c(0.0, 0.0);
  TTT.resize(lenntor, zero_c);
  omega.resize(lenntor, zero_c);
  omega_0.resize(lenntor, 0.0);
  omega_1.resize(lenntor, zero_c);
  // Calculate the WWW vector
  fd.field_ext_cls_calc_W(WWW, equ, *pt, fd.phik, fd.ntor1d);
  // Print the calculated WWW
  std::cout << "Calculated WWW: [";
  for (size_t i = 0; i < WWW.size(); ++i) {
    std::cout << WWW[i];
    if (i != WWW.size() - 1)
      std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  denskMkj_tot.resize(fd.ntotfem2d1f, {0.0, 0.0});
  jparkMkj_tot.resize(fd.ntotfem2d1f, {0.0, 0.0});

  if (imixvar == 1) {
    aparsk0.resize(fd.ntotfem2d1f, {0.0, 0.0});
    daparsk_sum.resize(fd.ntotfem2d1f, {0.0, 0.0});
    daparskdt.resize(fd.ntotfem2d1f, {0.0, 0.0});
  }

  if (rank == 0) {
    std::cout << ">> nsp = " << nsp << std::endl;
    std::cout << ">> mass_ion = " << mass_ion << std::endl;
    std::cout << ">> ntotfem2d1f = " << fd.getNtotfem2d1f() << std::endl;

    printf("----solve initial NJAP----\n");
  }

  // 计算初始的N, J, A, P
  record(0);

  // 检查是否是restart
  bool ptirestart = false;
  for (int fsc = 0; fsc < nsp; ++fsc) {
    // ptirestart |= pt[fsc].irestart;
  }

  if (fd.irestart || ptirestart) {
    record2h5(0);
  } else {
    *nstart = 0;
  }

  if (rank == 0) {
    printf("Already run step #: %d\n", *nstart);
    printf("time_all=%e, time_save=%e\n", time_all, time_save);
  }
}
