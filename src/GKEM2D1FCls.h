#ifndef GKEM2D1FCLS_H
#define GKEM2D1FCLS_H

#include "Equilibrium.h"
#include "FieldExt.h"
#include "MPIManager.h"
#include "Particle.h"
#include "ParticleExtCls.h"
#include <memory>

#define CHECK_NAN(x, label)                                                    \
  if (std::isnan(std::real(x)) || std::isnan(std::imag(x)))                    \
    std::cout << "Rank " << rank << ": NaN in " << label << " at step "        \
              << step << std::endl;

class GKEM2D1FCls {
private:
  MPIManager &mpiManager;

public:
  // 基本控制参数
  int *nstart = nullptr;
  int nsp, imixvar, ipullback;
  bool nl_pt_perturbed, nl_fd_solve;
  bool irestart = false;
  int nrun, ischeme;
  double dtoTA, dtoTN;
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
  std::vector<std::complex<double>> amp0, damp_sum, dampdt;
  std::vector<std::complex<double>> denskMkj_tot, jparkMkj_tot;
  std::vector<std::complex<double>> aparsk0, daparsk_sum, daparskdt;

  // 计时
  double t1 = 0.0, t2 = 0.0;

  //
  int ntrack = 10; // 轨道数

  // 构造与析构
  GKEM2D1FCls();
  GKEM2D1FCls(int argc, char **argv)
      : mpiManager(MPIManager::getInstance(argc, argv)),
        rank(mpiManager.getRank()), size(mpiManager.getSize()) {}
  ~GKEM2D1FCls() {
    if (nstart) {
      delete nstart; // 释放nstart指针
      nstart = nullptr;
    }
    if (pt) {
      pt.reset(); // 释放pt指针
    }
    std::cout << "GKEM2D1FCls destructor called" << std::endl;
  }

  // 成员函数
  void gkem_cls_initialize();
  void gkem_cls_finalize();
  void gkem_cls_readInput(const std::string &inputFile);
  void gkem_cls_record(int irun) {
    int irun_rec = irun + *nstart;

    // 写amplitude数据
    if (rank == 0) {
      write_amplitude(irun);
      // 写粒子轨迹数据（假设变量都已经准备好）
      for (int fsc = 0; fsc < nsp; ++fsc) {
        ParticleSpecies &species = this->pt->group.getSpecies(fsc);
        int ntrack_in = ntrack;
        ntrack_in = std::max(ntrack_in, species.getNptot());
        write_particle(irun_rec, species, ntrack_in);
      }
    }
  };
  void write_amplitude(int irun) {
    if (rank != 0)
      return;

    std::string sfile = "amplitude.txt";
    std::ofstream outfile;

    if (irun == 0 && !irestart)
      outfile.open(sfile, std::ios::out);
    else
      outfile.open(sfile, std::ios::app);

    for (const auto &amp : fd.amplitude_arr) {
      outfile << std::scientific << std::setprecision(15) << amp << " ";
    }
    outfile << std::endl;
    outfile.close();
  };
  void write_particle(int irun, ParticleSpecies &species, int ntrack) {
    // C++文件操作
    std::ofstream outfile;

    // 启动计时
    // timer.clsTic(4);

    if (rank == 0) {
      int nfile = 102;
      // 使用安全的 std::ostringstream 代替 sprintf
      std::ostringstream oss;
      // oss << std::setw(10) << sp_id; // 格式化 sp_id，宽度为 10
      oss << species.getSpId();
      std::string cfile = oss.str(); // 将结果存储为字符串
      std::string sfile = "data_track" + cfile + ".txt"; // 构造最终的文件名

      // 根据irun和irestart来选择文件打开方式
      if (irun == 0 && !irestart) {
        outfile.open(sfile, std::ios::out); // 写模式
      } else {
        outfile.open(sfile, std::ios::app); // 追加模式
      }

      // 循环处理所有轨道
      const double Bref = equ.getBref();
      const double rhoN = equ.rhoN;
      double mass = species.getMass();
      double zcharge = species.getCharge();
      ParticleCoords &coords = species.getCoords();
      std::vector<double> partmu = species.partmu;
      for (int fic = 0; fic < ntrack; ++fic) {
        int fidx = fic; // itrack[fic];
        double ptB = equ.getB(coords.partrad[fidx], coords.parttheta[fidx]);

        double ptpsi_can = equ.getpsi(coords.partrad[fidx]) +
                           coords.partvpar[fidx] * rhoN *
                               equ.getfpol_r(coords.partrad[fidx]) * mass *
                               Bref / (zcharge * ptB);

        double ptenergy =
            mass * (coords.partvpar[fidx] * coords.partvpar[fidx] / 2 +
                    partmu[fidx] * ptB);

        // 输出数据
        outfile << std::scientific << std::setprecision(16)
                << coords.partrad[fidx] << " " << coords.parttheta[fidx] << " "
                << coords.partphitor[fidx] << " " << coords.partvpar[fidx]
                << " " << partmu[fidx] << " " << coords.partw[fidx] << " "
                << equ.getR(coords.partrad[fidx], coords.parttheta[fidx]) << " "
                << equ.getZ(coords.partrad[fidx], coords.parttheta[fidx]) << " "
                << ptenergy << " " << ptpsi_can << std::endl;
      }

      // 关闭文件
      outfile.close();
    }
  };
  void gkem_cls_record2h5(int istep){};
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
      const std::vector<std::complex<double>> &amp,
      //  std::vector<ParticleCoords> &xv0,
      const std::vector<std::vector<double>> &partfog0_allsp,
      const std::vector<ParticleCoords> &dxvdt) {
    // 计算delta omega

    std::vector<std::complex<double>> TTT_tmp;
    TTT_tmp.resize(lenntor, std::complex<double>(0.0, 0.0));
    TTT.resize(lenntor, std::complex<double>(0.0, 0.0));
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
          /*icase*/ 0, fd.phik, fd.ntor1d, amp);

      for (int i = 0; i < lenntor; ++i) {
        if (!std::isfinite(TTT_tmp[i].real()) ||
            !std::isfinite(TTT_tmp[i].imag())) {
          std::cerr << "Rank " << rank << ": TTT_tmp[" << i
                    << "] = " << TTT_tmp[i] << " is NaN at species = " << fsc
                    << std::endl;
        }
      }

      // 将TTT_tmp的值累加到TTT中
      for (int i = 0; i < lenntor; ++i) {
        TTT[i] += TTT_tmp[i];
      }
    } // nsp

    for (int i = 0; i < lenntor; ++i) {
      if (!std::isfinite(TTT[i].real()) || !std::isfinite(TTT[i].imag())) {
        std::cout << "Rank " << rank << " | TTT[" << i << "] = " << TTT[i]
                  << " (invalid) at step " << std::endl;
      }

      if (std::abs(WWW[i]) < 1e-12 || !std::isfinite(WWW[i])) {
        std::cout << "Rank " << rank << " | WWW[" << i << "] = " << WWW[i]
                  << " (invalid) at step " << std::endl;
      }
    }

    for (int i = 0; i < lenntor; ++i) {
      {
        std::cout << "before MPI sum calculated T: " << std::endl;
        std::cout << "Rank " << rank << " | TTT[" << i << "] = " << TTT[i]
                  << std::endl;
      }
    }

    // MPI 通信
    MPI_Allreduce(MPI_IN_PLACE, TTT.data(), lenntor, MPI_DOUBLE_COMPLEX,
                  MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < lenntor; ++i) {
      {
        std::cout << "after MPI sum calculated T: " << std::endl;
        std::cout << "Rank " << rank << " | TTT[" << i << "] = " << TTT[i]
                  << std::endl;
      }
    }

    // 计算omega
    for (int i = 0; i < lenntor; ++i) {
      omega_1_out[i] =
          std::complex<double>(0.0, 1.0) * TTT[i] / 2.0 / WWW[i]; // 计算omega_1
                                                                  // 更新omega
    }
  }

  void add_dAdt(std::vector<std::complex<double>> &amp_tmp,
                const std::vector<std::complex<double>> &omega_1_tmp,
                const double dt) {
    //
    std::vector<std::complex<double>> dA;
    dA.resize(lenntor, std::complex<double>(0.0, 0.0));
    std::complex<double> i_c(0.0, 1.0);
    for (size_t itor = 0; itor < amp_tmp.size(); ++itor) {
      // dA/dt = - i * omega_1 * A
      dA[itor] = -i_c * omega_1_tmp[itor] * amp_tmp[itor]; // 计算dA
      amp_tmp[itor] += dA[itor] * dt;
    }
  };
  void onestep_rk4(std::vector<ParticleCoords> &xv0,
                   std::vector<ParticleCoords> &dxv_sum,
                   std::vector<ParticleCoords> &dxvdt, double dt,
                   double time_now) {
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
    std::vector<std::complex<double>> amp_tmp;
    amp_tmp = fd.amplitude_arr; // 复制amplitude_arr到amp_tmp
    amp0 = fd.amplitude_arr;

    std::vector<double> fast_ossilation_phase;
    std::vector<std::complex<double>> amp0_with_phase;
    amp0_with_phase.resize(lenntor, zero_c);
    fast_ossilation_phase.resize(lenntor, 0.0);
    for (int itor = 0; itor < lenntor; ++itor) {
      fast_ossilation_phase[itor] = std::cos(omega_0[itor] * time_now);
      amp0_with_phase[itor] =
          fd.amplitude_arr[itor] *
          std::complex<double>(fast_ossilation_phase[itor], 0.0);
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
    gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, amp0, partfog0_allsp,
                               dxvdt);
    print_omega1(omega_1_tmp, rank);
    // add_dAdt(amp_tmp, omega_1_tmp, dthalf);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp0_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    // Step 2
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);

    // gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, amp_tmp,
    // partfog0_allsp,
    //                            dxvdt);
    // print_omega1(omega_1_tmp, rank);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp0_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    // Step 3
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);
    // gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, amp_tmp,
    // partfog0_allsp,
    //                            dxvdt);
    // print_omega1(omega_1_tmp, rank);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp0_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    // // Step 4
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dt);

    // gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, amp_tmp,
    // partfog0_allsp,
    //                            dxvdt);
    // print_omega1(omega_1_tmp, rank);
    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp0_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    // Final update of aparsk and particle coordinates

    pt->particle_setvalue2sp_from_coords(*pt, xv0);
    pt->particle_add_coords2sp(*pt, dxv_sum, 1.0);

    for (size_t itor = 0; itor < omega_1_tmp.size(); ++itor) {
      amp_tmp[itor] *=
          std::exp(std::complex<double>(0.0, -1.0) * dt * omega_1_tmp[itor]);
      // amp_tmp[itor] *=
      //     std::exp(std::complex<double>(0.0, -1.0) * omega_0[itor] * dt);
    }
    fd.amplitude_arr = amp_tmp; // 更新fd.amplitude_arr
  };
  void onestep_rk4_testParticle(std::vector<ParticleCoords> &xv0,
                                std::vector<ParticleCoords> &dxv_sum,
                                std::vector<ParticleCoords> &dxvdt, double dt,
                                double time_now) {
    double dthalf = dt / 2.0;
    double dt1o6 = dt / 6.0;
    double dt2o6 = dt * 2.0 / 6.0;
    std::complex<double> zero_c(0.0, 0.0);
    std::vector<double> fast_ossilation_phase;
    std::vector<std::complex<double>> amp_with_phase;
    amp_with_phase.resize(lenntor, zero_c);
    fast_ossilation_phase.resize(lenntor, 0.0);
    for (int itor = 0; itor < lenntor; ++itor) {
      fast_ossilation_phase[itor] = std::cos(omega_0[itor] * time_now);
      amp_with_phase[itor] =
          fd.amplitude_arr[itor] *
          std::complex<double>(fast_ossilation_phase[itor], 0.0);
    }

    for (int fsc = 0; fsc < nsp; ++fsc) {
      ParticleSpecies &species = pt->group.getSpecies(fsc);
      int nptot = species.getNptot();
      dxv_sum[fsc].initialize(nptot); // dxv_sum initialize to 0.0
      xv0[fsc] = species.getCoords(); // xv0= pt
    }
    std::vector<std::complex<double>> amp_tmp;
    amp_tmp = fd.amplitude_arr; // 复制amplitude_arr到amp_tmp
    amp0 = fd.amplitude_arr;

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
    gkem_cls_solve_delta_omega(omega_0, omega_1_tmp, amp0, partfog0_allsp,
                               dxvdt);
    print_omega1(omega_1_tmp, rank);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    // Step 2
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    // Step 3
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    // // Step 4
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dt);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase, xv0,
        partmu0_allsp, partfog0_allsp, dxvdt);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    // Final update of aparsk and particle coordinates

    pt->particle_setvalue2sp_from_coords(*pt, xv0);
    pt->particle_add_coords2sp(*pt, dxv_sum, 1.0);
  };

  void onestep_euler();
  void onestep_heun();
  void test() {
    gkem_cls_initialize();
    // 测试函数
    if (rank == 0) {
      std::cout << "Testing GKEM2D1FCls with nrun = " << nrun << std::endl;
      if (rank == 0) {
        std::cout << "Initial amplitude---." << std::endl;
        std::cout << "amplitude_arr: ";
        for (const auto &amp : fd.amplitude_arr) {
          std::cout << std::scientific << std::setprecision(15) << amp << " ";
        }
        std::cout << std::endl;
      }
    }
    double time_now = 0.0; // 初始化时间
    for (int i = 0; i < nrun; ++i) {
      if (rank == 0) {
        std::cout << "Running step " << i + 1 << " of " << nrun << std::endl;
      }
      onestep_rk4(xv0, dxv_sum, dxvdt, dtoTA, time_now);
      if (rank == 0) {
        std::cout << "Step " << i + 1 << " completed." << std::endl;
        std::cout << std::endl;
      }
      gkem_cls_record(i + 1);
      time_now += dtoTA; // 更新时间
    }
    if (rank == 0) {
      std::cout << "Test completed." << std::endl;
    }
  };
  void testParticle() {
    gkem_cls_initialize();
    // 测试函数
    if (rank == 0) {
      std::cout << "Testing GKEM2D1FCls with nrun = " << nrun << std::endl;
      if (rank == 0) {
        std::cout << "Initial amplitude---." << std::endl;
        std::cout << "amplitude_arr: ";
        for (const auto &amp : fd.amplitude_arr) {
          std::cout << std::scientific << std::setprecision(15) << amp << " ";
        }
        std::cout << std::endl;
      }
    }
    double time_now = 0.0; // 初始化时间
    for (int i = 0; i < nrun; ++i) {
      if (rank == 0) {
        std::cout << "Running step " << i + 1 << " of " << nrun << std::endl;
      }
      onestep_rk4_testParticle(xv0, dxv_sum, dxvdt, dtoTA, time_now);
      if (rank == 0) {
        std::cout << "Step " << i + 1 << " completed." << std::endl;
        std::cout << std::endl;
      }
      gkem_cls_record(i + 1);
    }
    if (rank == 0) {
      std::cout << "Test completed." << std::endl;
    }
    time_now += dtoTA; // 更新时间
  };
};
#endif // GKEM2D1FCLS_H
void GKEM2D1FCls::gkem_cls_readInput(const std::string &inputFile) {
  INIReader reader(inputFile);

  // Check if reading was successful
  if (reader.ParseError() < 0) {
    throw std::runtime_error("Error: Unable to open or parse input file '" +
                             inputFile + "'");
  }
  std::cout << "GKEM2D1FCls readInput: " << inputFile << std::endl;
  // Read key-value pairs
  nrun = reader.GetInteger("MENG", "nrun", 2);
  dtoTA = reader.GetReal("MENG", "dtoTA", 0.01);
}
void GKEM2D1FCls::gkem_cls_initialize() {

  // read input
  gkem_cls_readInput("input.ini");

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
    printf("   ischeme=%d\n", ischeme);
  }

  // 分配粒子
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
  this->omega_0 = fd.omega0; // 直接从fd获取omega_0
  // for perturbations
  WWW.resize(lenntor, 0.0);
  std::complex<double> zero_c(0.0, 0.0);
  TTT.resize(lenntor, zero_c);
  omega.resize(lenntor, zero_c);
  omega_0.resize(lenntor, 1.0);
  omega_1.resize(lenntor, zero_c);

  amp0.resize(lenntor, std::complex<double>(0.0, 0.0));
  damp_sum.resize(lenntor, std::complex<double>(0.0, 0.0));
  dampdt.resize(lenntor, std::complex<double>(0.0, 0.0));

  // Calculate the WWW vector
  fd.field_ext_cls_calc_W(WWW, equ, *pt, fd.phik, fd.ntor1d, fd.amplitude_arr);
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
  gkem_cls_record(0);

  // 检查是否是restart
  bool ptirestart = false;
  for (int fsc = 0; fsc < nsp; ++fsc) {
    // ptirestart |= pt[fsc].irestart;
  }

  if (fd.irestart || ptirestart) {
    gkem_cls_record2h5(0);
  } else {
    *nstart = 0;
  }

  if (rank == 0) {
    printf("Already run step #: %d\n", *nstart);
    printf("time_all=%e, time_save=%e\n", time_all, time_save);
  }
}
