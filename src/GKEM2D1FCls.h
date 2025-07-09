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
  int lenntor; // ntor1d.size()

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
  void print_omega1(const std::vector<std::complex<double>> &omega_1_tmp,
                    int rank) {
    if (rank == 0) {
      std::cout << "omega_1: ";
      for (const auto &val : omega_1_tmp) {
        std::cout << std::scientific << std::setprecision(15) << "("
                  << val.real() << "," << val.imag() << " i ) ";
      }
      std::cout << std::endl;
    }
  };
  void gkem_cls_solve_delta_omega(
      std::vector<double> &omega_0_in,
      std::vector<std::complex<double>> &omega_1_out,
      //  std::vector<ParticleCoords> &xv0,
      const std::vector<std::vector<double>> &partfog0_allsp,
      const std::vector<ParticleCoords> &dxvdt) {
    // 计算delta omega

    std::vector<std::complex<double>> TTT_tmp;
    TTT_tmp.resize(lenntor);
    for (int fsc = 0; fsc < nsp; ++fsc) {
      std::vector<double> partfog0 = partfog0_allsp[fsc];
      ParticleSpecies &species = this->pt->group.getSpecies(fsc);
      ParticleCoords &coords = species.getCoords();
      int nptot = species.getNptot();

      fd.field_ext_cls_calc_T_oneSpecies(
          TTT_tmp, equ, species, coords.partrad, coords.parttheta,
          coords.partphitor, coords.partvpar, coords.partw, partfog0,
          dxvdt[fsc].partrad, dxvdt[fsc].parttheta, dxvdt[fsc].partphitor,
          dxvdt[fsc].partvpar, dxvdt[fsc].partw,
          /*icase*/ 0, fd.phik, fd.ntor1d);
      // 将TTT_tmp的值累加到TTT中
      for (int i = 0; i < lenntor; ++i) {
        TTT[i] += TTT_tmp[i];
      }
    } // nsp

    // MPI 通信
    MPI_Allreduce(MPI_IN_PLACE, TTT.data(), lenntor, MPI_DOUBLE_COMPLEX,
                  MPI_SUM, MPI_COMM_WORLD);
    // 计算omega
    for (int i = 0; i < lenntor; ++i) {
      omega_1_out[i] =
          std::complex<double>(0.0, 1.0) * TTT[i] / 2.0 / WWW[i]; // 计算omega_1
                                                                  // 更新omega
    }
  }
  void onestep_rk4(std::vector<ParticleCoords> &xv0,
                   std::vector<ParticleCoords> &dxv_sum,
                   std::vector<ParticleCoords> &dxvdt, double dt) {
    double dthalf = dt / 2.0;
    double dt1o6 = dt / 6.0;
    double dt2o6 = dt * 2.0 / 6.0;
    std::complex<double> zero_c(0.0, 0.0);

    for (int fsc = 0; fsc < nsp; ++fsc) {
      ParticleSpecies &species = pt->group.getSpecies(fsc);
      int nptot = species.getNptot();
      dxv_sum[fsc].initialize(nptot); // dxv_sum initialize to 0.0
      xv0[fsc] = species.getCoords(); // xv0= pt
    }

    // Step 1
    std::vector<std::vector<double>> partfog0_allsp, partmu0_allsp;
    partfog0_allsp.resize(nsp);
    partmu0_allsp.resize(nsp);
    for (int fsc = 0; fsc < nsp; ++fsc) {
      ParticleSpecies &species = this->pt->group.getSpecies(fsc);
      partfog0_allsp[fsc] = species.partfog;
      partmu0_allsp[fsc] = species.partmu;
    }
    std::vector<std::complex<double>> omega_1_tmp(lenntor, zero_c);
    gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, partfog0_allsp, dxvdt);
    print_omega1(omega_1_tmp, rank);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(equ, fd, fd.phik, fd.apark,
                                               fd.ntor1d, xv0, partmu0_allsp,
                                               partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    // Step 2
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);

    gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, partfog0_allsp, dxvdt);
    print_omega1(omega_1_tmp, rank);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(equ, fd, fd.phik, fd.apark,
                                               fd.ntor1d, xv0, partmu0_allsp,
                                               partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    // Step 3
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);
    gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, partfog0_allsp, dxvdt);
    print_omega1(omega_1_tmp, rank);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(equ, fd, fd.phik, fd.apark,
                                               fd.ntor1d, xv0, partmu0_allsp,
                                               partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    // // Step 4
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dt);

    gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, partfog0_allsp, dxvdt);
    print_omega1(omega_1_tmp, rank);
    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(equ, fd, fd.phik, fd.apark,
                                               fd.ntor1d, xv0, partmu0_allsp,
                                               partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    // Final update of aparsk and particle coordinates

    pt->particle_setvalue2sp_from_coords(*pt, xv0);
    pt->particle_add_coords2sp(*pt, dxv_sum, 1.0);
  };

  void onestep_euler();
  void onestep_heun();
  void calc_T_allSpecies();
  void test(int nrun_in) {
    nrun = nrun_in;
    initialize();
    // 测试函数
    if (rank == 0) {
      std::cout << "Testing GKEM2D1FCls with nrun = " << nrun << std::endl;
    }
    for (int i = 0; i < nrun; ++i) {
      if (rank == 0) {
        std::cout << "Running step " << i + 1 << " of " << nrun << std::endl;
      }
      onestep_rk4(xv0, dxv_sum, dxvdt, dtoTA);
    }
    if (rank == 0) {
      std::cout << "Test completed." << std::endl;
    }
  };
};
#endif // GKEM2D1FCLS_H
void GKEM2D1FCls::initialize() {

  // read input
  dtoTA = 1.0;
  nrun = 2;

  //
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
  lenntor = fd.lenntor;
  // for perturbations
  WWW.resize(lenntor, 0.0);
  std::complex<double> zero_c(0.0, 0.0);
  TTT.resize(lenntor, zero_c);
  omega.resize(lenntor, zero_c);
  omega_0.resize(lenntor, 1.0);
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

  denskMkj_tot.resize(fd.ntotfem2d1f, zero_c);
  jparkMkj_tot.resize(fd.ntotfem2d1f, zero_c);

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
