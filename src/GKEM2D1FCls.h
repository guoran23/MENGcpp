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
  int nsp;
  int imixvar = 0;
  int ipullback = 0;
  bool nl_pt_perturbed, nl_fd_solve;
  bool irestart = false;
  int nrun, ischeme;
  double dtoTA, dtoTN;
  double time_all, time_save;
  int rank, size;
  int lenntor;   // ntor1d.size()
  int itest = 0; // 0: test particle motion
                 // 1: euler, 2: rk4
  bool idiag_dxvdt;
  int idiag_dxvdt_step;

  // 物理对象
  Equilibrium equ;
  FieldExtCls fd;
  std::unique_ptr<ParticleExtCls> pt;

  // 工作空间
  std::vector<double> WWW;
  std::vector<std::complex<double>> TTT;
  std::vector<std::complex<double>> TTT_global;
  std::vector<std::complex<double>> omega;
  std::vector<double> omega_0;
  std::vector<std::complex<double>> omega_1;
  std::vector<std::complex<double>> amp_with_phase_diag;
  // for RK methods
  std::vector<ParticleCoords> xv0, dxv_sum, dxvdt;
  std::vector<std::complex<double>> amp0, damp_sum, dampdt;

  std::vector<std::complex<double>> denskMkj_tot, jparkMkj_tot;
  double particle_tot_Energy = 0.0;
  double particle_tot_pert_Energy = 0.0;

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
  void onestep_heun();
  void gkem_cls_record2h5(int istep);

  void gkem_cls_record(int irun) {
    int irun_rec = irun + *nstart;
    calc_particle_tot_Energy(); // 计算粒子总能量
    // 写amplitude数据
    if (rank == 0) {
      write_amplitude(irun);
      write_amp_with_phase(irun);
      write_omega1(irun);
      write_particle_tot_Energy(irun);

      if (idiag_dxvdt && irun % idiag_dxvdt_step == 0) {
        write_species_particles("data_xv0_", xv0, irun);
        write_species_particles("data_dxvdt_", dxvdt, irun);
      }
      // 写粒子轨迹数据（假设变量都已经准备好）
      for (int fsc = 0; fsc < nsp; ++fsc) {
        ParticleSpecies &species = this->pt->group.getSpecies(fsc);
        int ntrack_in = ntrack;
        ntrack_in = std::min(ntrack_in, species.getNptot());
        write_particle(irun_rec, species, ntrack_in);
      }
    }
  }

  template <typename Container>
  void write_data(const std::string &filename, const Container &data,
                  int irun) {
    if (rank != 0)
      return;

    std::ofstream outfile;
    if (irun == 0 && !irestart)
      outfile.open(filename, std::ios::out);
    else
      outfile.open(filename, std::ios::app);

    for (const auto &val : data) {
      double re = std::real(val);
      double im = std::imag(val);

      if (std::isnan(re) || std::isnan(im)) {
        outfile << "Nan Nan ";
      } else {
        outfile << std::scientific << std::setprecision(15) << re << " " << im
                << " ";
      }
    }
    outfile << std::endl;
    outfile.close();
  }
  void write_amplitude(int irun) {
    write_data("data_Amplitude.txt", fd.amplitude_arr, irun);
  }

  void write_omega1(int irun) {
    write_data("data_omega1.txt", this->omega_1, irun);
  }

  void write_amp_with_phase(int irun) {
    write_data("data_amp_with_phase.txt", this->amp_with_phase_diag, irun);
  }

  void calc_particle_tot_Energy() {
    // 计算粒子总能量
    particle_tot_Energy = 0.0;
    particle_tot_pert_Energy = 0.0;
    for (int fsc = 0; fsc < nsp; ++fsc) {
      ParticleSpecies &species = this->pt->group.getSpecies(fsc);
      ParticleCoords &coords = species.getCoords();
      int nptot = species.getNptot();

      for (int fic = 0; fic < nptot; ++fic) {
        double ptB = equ.getB(coords.partrad[fic], coords.parttheta[fic]);
        double ptenergy = species.getMass() *
                          (coords.partvpar[fic] * coords.partvpar[fic] / 2 +
                           species.partmu[fic] * ptB);
        particle_tot_pert_Energy += ptenergy * coords.partw[fic];
        particle_tot_Energy += ptenergy * species.partfog[fic];
      }
    }
    // MPI
    MPI_Allreduce(MPI_IN_PLACE, &particle_tot_Energy, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &particle_tot_pert_Energy, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
  }

  void write_particle_tot_Energy(int irun) {
    if (rank != 0)
      return;

    std::string sfile = "particle_tot_energy.txt";
    std::ofstream outfile;

    if (irun == 0 && !irestart)
      outfile.open(sfile, std::ios::out);
    else
      outfile.open(sfile, std::ios::app);

    outfile << std::scientific << std::setprecision(15) << particle_tot_Energy
            << " " << particle_tot_pert_Energy << std::endl;
    outfile.close();
  }

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
  }
  void write_species_particles(const std::string &prefix,
                               const std::vector<ParticleCoords> &species_data,
                               int irun) {
    for (size_t sp = 0; sp < species_data.size(); ++sp) {
      const auto &p = species_data[sp];

      // 构建文件名，例如 data_xv0_0.txt
      std::string filename = prefix + "SP_" + std::to_string(sp) + ".txt";
      std::ofstream outfile;

      if (irun == 0)
        outfile.open(filename, std::ios::out);
      else
        outfile.open(filename, std::ios::app);

      size_t np = p.partrad.size(); // 假设所有分量长度相等
      for (size_t i = 0; i < np; ++i) {
        outfile << std::scientific << std::setprecision(15) << p.partrad[i]
                << " " << p.parttheta[i] << " " << p.partphitor[i] << " "
                << p.partvpar[i] << " " << p.partw[i] << "\n";
      }

      outfile.close();
    }
  }

  void print_complex_vec(const std::string &name,
                         const std::vector<std::complex<double>> &omega_1_tmp,
                         int rank) {
    if (rank == 0) {
      std::cout << name << " = ";
      for (const auto &val : omega_1_tmp) {
        std::cout << std::scientific << std::setprecision(15) << "("
                  << val.real() << "," << val.imag() << " i ) ";
      }
      std::cout << std::endl;
    }
  }

  void
  gkem_cls_solve_delta_omega(std::vector<std::complex<double>> &omega_1_out,
                             std::vector<std::complex<double>> &amp) {
    // MPI 通信
    constexpr std::complex<double> zero_c(0.0, 0.0);
    std::fill(TTT_global.begin(), TTT_global.end(), zero_c);
    MPI_Allreduce(TTT.data(), TTT_global.data(), lenntor, MPI_DOUBLE_COMPLEX,
                  MPI_SUM, MPI_COMM_WORLD);

    // 计算omega
    double rhoN = equ.rhoN;
    for (int i = 0; i < lenntor; ++i) {
      omega_1_out[i] = std::complex<double>(0.0, 1.0) * TTT_global[i] / 2.0 /
                       WWW[i] / std::norm(amp[i]); // 计算omega
      omega_1_out[i] =  - omega_1_out[i] / 6.28;
    }
    // std::cout << "-----rank = " << rank << std::endl;
    // print_complex_vec("TTT", TTT, 0);
    // print_complex_vec("TTT_global", TTT_global, 0);
    // print_complex_vec("omega1=", omega_1_out,0);
    // std::cout << "-----rank =" << rank << " end-----" << std::endl;
  }

  std::vector<std::complex<double>>
  add_real_and_complex(const std::vector<double> &real_vec,
                       const std::vector<std::complex<double>> &complex_vec) {

    size_t n = real_vec.size();
    std::vector<std::complex<double>> result(n);
    for (size_t i = 0; i < n; ++i) {
      result[i] = std::complex<double>(real_vec[i], 0.0) + complex_vec[i];
    }
    return result;
  }

  void calc_amp_with_phase(std::vector<std::complex<double>> &amp_with_phase,
                           const std::vector<std::complex<double>> &amp_in,
                           const double time) {
    constexpr std::complex<double> i_c(0.0, 1.0);
    constexpr std::complex<double> zero_c(0.0, 0.0);
    std::vector<std::complex<double>> fast_ossilation_phase(lenntor, zero_c);
    amp_with_phase.resize(lenntor, zero_c);
    for (int itor = 0; itor < lenntor; ++itor) {
      fast_ossilation_phase[itor] =
          std::exp(-i_c * (this->omega_0[itor] * time));
      amp_with_phase[itor] =
          amp_in[itor] * fast_ossilation_phase[itor]; // 计算amplitude
    }
  }

  void add_dAdt(std::vector<std::complex<double>> &A_out,
                const std::vector<std::complex<double>> &dAdt,
                const double dt) {
    for (size_t itor = 0; itor < A_out.size(); ++itor) {
     
      if (this->itest == 0) {
        A_out[itor] +=  0.0;
      } else {
         A_out[itor] += dAdt[itor] * dt;
      }
    }

  }

  void calc_dAdt(std::vector<std::complex<double>> &dAdt,
                 const std::vector<std::complex<double>> &amp_tmp,
                 const std::vector<std::complex<double>> &omega_1_tmp) {
    //
    dAdt.resize(lenntor, std::complex<double>(0.0, 0.0));
    constexpr std::complex<double> i_c = {0.0, 1.0};
    for (size_t itor = 0; itor < amp_tmp.size(); ++itor) {
      // dA/dt = - i * omega_1 * A
      dAdt[itor] = -i_c * omega_1_tmp[itor] * amp_tmp[itor]; // 计算dA
    }
  }

  void onestep_rk4(std::vector<ParticleCoords> &xv0,
                   std::vector<ParticleCoords> &dxv_sum,
                   std::vector<ParticleCoords> &dxvdt, double dt,
                   double time_now) {
    double dthalf = dt / 2.0;
    double dt1o6 = dt / 6.0;
    double dt2o6 = dt * 2.0 / 6.0;
    constexpr std::complex<double> zero_c(0.0, 0.0);

    for (int fsc = 0; fsc < nsp; ++fsc) {
      ParticleSpecies &species = pt->group.getSpecies(fsc);
      int nptot = species.getNptot();
      dxv_sum[fsc].initialize(nptot); // dxv_sum initialize to 0.0
      xv0[fsc] = species.getCoords(); // xv0= pt
    }
    std::fill(damp_sum.begin(), damp_sum.end(), zero_c);

    std::vector<std::complex<double>> amp_tmp;
    amp_tmp = fd.amplitude_arr; // 复制amplitude_arr到amp_tmp
    amp0 = fd.amplitude_arr;

    std::vector<std::complex<double>> amp0_with_phase(lenntor, zero_c);
    std::vector<std::complex<double>> amp_with_phase_tmp(lenntor, zero_c);
    calc_amp_with_phase(amp0_with_phase, amp0, time_now);
    amp_with_phase_tmp = amp0_with_phase;
    // print_complex_vec("amp0_with_phase", amp0_with_phase, 0);
    this->amp_with_phase_diag = amp0_with_phase;
    std::vector<std::complex<double>> omega_1_tmp(lenntor, zero_c);

    // Step 1, at t = 0
    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp0_with_phase, dxvdt, TTT);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    gkem_cls_solve_delta_omega(omega_1_tmp, amp0);
    calc_dAdt(dampdt, amp0, omega_1_tmp);
    add_dAdt(damp_sum, dampdt, dt1o6);

    // Step 2
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);
    add_dAdt(amp_tmp, dampdt, dthalf);
    calc_amp_with_phase(amp_with_phase_tmp, amp_tmp, time_now + dthalf);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase_tmp, dxvdt, TTT);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    gkem_cls_solve_delta_omega(omega_1_tmp, amp_tmp);
    calc_dAdt(dampdt, amp_tmp, omega_1_tmp);
    add_dAdt(damp_sum, dampdt, dt2o6);

    // Step 3
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dthalf);

    amp_tmp = amp0;
    add_dAdt(amp_tmp, dampdt, dthalf);
    calc_amp_with_phase(amp_with_phase_tmp, amp_tmp, time_now + dthalf);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase_tmp, dxvdt, TTT);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt2o6);

    gkem_cls_solve_delta_omega(omega_1_tmp, amp_tmp);
    calc_dAdt(dampdt, amp_tmp, omega_1_tmp);
    add_dAdt(damp_sum, dampdt, dt2o6);

    // // Step 4
    pt->particle_setvalue2sp_from_coords(*pt, xv0); // reset to initial
    pt->particle_add_coords2sp(*pt, dxvdt, dt);

    amp_tmp = amp0;
    add_dAdt(amp_tmp, dampdt, dt);
    calc_amp_with_phase(amp_with_phase_tmp, amp_tmp, time_now + dt);

    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp_with_phase_tmp, dxvdt, TTT);
    pt->particle_coords_cls_axpy2sp(dxv_sum, dxvdt, dt1o6);

    gkem_cls_solve_delta_omega(omega_1_tmp, amp_tmp);
    calc_dAdt(dampdt, amp_tmp, omega_1_tmp);
    add_dAdt(damp_sum, dampdt, dt1o6);

    // Final update of Amp and particle coordinates
    pt->particle_setvalue2sp_from_coords(*pt, xv0);
    pt->particle_add_coords2sp(*pt, dxv_sum, 1.0);
    add_dAdt(fd.amplitude_arr, damp_sum, 1.0);
    this->omega_1 = omega_1_tmp; // update omega_1 for recording
  }

  void onestep_euler(std::vector<ParticleCoords> &xv0,
                     std::vector<ParticleCoords> &dxvdt, double dt,
                     double time_now) {

    constexpr std::complex<double> zero_c(0.0, 0.0);

    for (int fsc = 0; fsc < nsp; ++fsc) {
      ParticleSpecies &species = pt->group.getSpecies(fsc);
      int nptot = species.getNptot();
      xv0[fsc] = species.getCoords(); // xv0= pt
    }
    std::vector<std::complex<double>> amp_tmp(lenntor, zero_c);
    amp_tmp = fd.amplitude_arr; // 复制amplitude_arr到amp_tmp
    amp0 = fd.amplitude_arr;

    std::vector<std::complex<double>> amp0_with_phase(lenntor, zero_c);
    calc_amp_with_phase(amp0_with_phase, amp0, time_now);
    this->amp_with_phase_diag = amp0_with_phase;

    // Step 1
    pt->particle_ext_cls_dxvpardt123EM2d1f_2sp(
        equ, fd, fd.phik, fd.apark, fd.ntor1d, amp0_with_phase, dxvdt, TTT);
    std::vector<std::complex<double>> omega_1_tmp(lenntor, zero_c);
    gkem_cls_solve_delta_omega(omega_1_tmp, amp0);
    calc_dAdt(dampdt, amp0, omega_1_tmp);
    pt->particle_add_coords2sp(*pt, dxvdt, dt);
    add_dAdt(amp_tmp, dampdt, dt);

    fd.amplitude_arr = amp_tmp;  // 更新fd.amplitude_arr
    this->omega_1 = omega_1_tmp; // 更新omega_1
  }

  void test() {
    gkem_cls_initialize();
    // 测试函数
    if (rank == 0) {
      std::cout << "Testing GKEM2D1FCls with nrun = " << nrun << std::endl;
      std::cout << "dtoTN = " << dtoTN << std::endl;
      std::cout << "itest = " << itest << std::endl;
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
      if (itest == 0 || itest == 2) {
        onestep_rk4(xv0, dxv_sum, dxvdt, dtoTN, time_now);
      } else if (itest == 1) {
        onestep_euler(xv0, dxvdt, dtoTN, time_now);
      } else {
        throw std::runtime_error("Invalid itest number: " +
                                 std::to_string(itest));
      }

      gkem_cls_record(i + 1);
      time_now += dtoTN; // 更新时间
    }
    if (rank == 0) {
      std::cout << "Test completed." << std::endl;
    }
  }
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
  this->dtoTN = reader.GetReal("MENG", "dtoTN", 0.01);
  itest = reader.GetInteger("MENG", "itest", 0);
  idiag_dxvdt = reader.GetBoolean("MENG", "idiag_dxvdt", false);
  idiag_dxvdt_step = reader.GetInteger("Meng", "idiag_dxvdt_step", 100);
}

void GKEM2D1FCls::gkem_cls_initialize() {

  // read input
  gkem_cls_readInput("input.ini");

  //
  equ.readInput("input.ini");
  pt = std::make_unique<ParticleExtCls>(equ, rank, size); // 初始化粒子
  nsp = pt->getNsp(); // 获取粒子种类数
  double mass_bk = pt->mass_bk;
  std::vector<double> vts_vec(nsp, 0.0);
  std::vector<double> paux_T_transit_vec(nsp, 0.0);
  std::vector<double> twoPi_o_T_transit_vec(nsp, 0.0);
  for (int isp = 0; isp < nsp; ++isp) {
    ParticleSpecies &species = this->pt->group.getSpecies(isp);
    vts_vec[isp] = species.vts;
    paux_T_transit_vec[isp] = species.paux_T_transit;
    twoPi_o_T_transit_vec[isp] = 2.0 * M_PI / paux_T_transit_vec[isp];
  }

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
  // for perturbations
  WWW.resize(lenntor, 0.0);
  constexpr std::complex<double> zero_c(0.0, 0.0);
  omega_0.resize(lenntor, 0.0); // omega_0 is real vector
  omega.resize(lenntor, zero_c);
  omega_1.resize(lenntor, zero_c);
  TTT.resize(lenntor, zero_c);
  TTT_global.resize(lenntor, zero_c);

  this->omega_0 = fd.omega0;
  omega = add_real_and_complex(omega_0, omega_1);

  //
  amp0.resize(lenntor, std::complex<double>(0.0, 0.0));
  damp_sum.resize(lenntor, std::complex<double>(0.0, 0.0));
  dampdt.resize(lenntor, std::complex<double>(0.0, 0.0));
  amp_with_phase_diag.resize(lenntor, std::complex<double>(0.0, 0.0));

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

  // denskMkj_tot.resize(fd.ntotfem2d1f, zero_c);
  // jparkMkj_tot.resize(fd.ntotfem2d1f, zero_c);

  if (rank == 0) {
    std::cout << ">> nsp = " << nsp << std::endl;
    std::cout << ">> Background mass_bk = " << mass_bk
              << ", ion_densprof = " << pt->ion_densprof << std::endl;
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
    // gkem_cls_record2h5(0);
  } else {
    *nstart = 0;
  }

  if (rank == 0) {
    print_vector("omega_0", omega_0);
    std::cout << "nsp = " << nsp << std::endl;
    print_vector("vts_vec", vts_vec);
    print_vector("paux_T_transit_vec", paux_T_transit_vec);
    print_vector("twoPi_o_T_transit_vec", twoPi_o_T_transit_vec);
    printf("Already run step #: %d\n", *nstart);
    printf("time_all=%e, time_save=%e\n", time_all, time_save);
  }
}
