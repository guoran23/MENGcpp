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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "../inih/INIReaderExt.h"
#include "Equilibrium.h"
#include "Field.h"
#include "MPIManager.h"
#include "util_math.h"

#include <chrono> // for time-based seed
#include <cmath>
#include <functional> // std::plus, std::multiplies
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <random>
#include <string>
#include <vector>

class ParticleCoords {
private:
  int nptot;

public:
  std::vector<double> partrad, parttheta, partphitor;
  std::vector<double> partvpar;
  std::vector<double> partw;
  // Constructor
  ParticleCoords() = default;
  ParticleCoords(int total_particles)
      : nptot(total_particles), partrad(total_particles),
        parttheta(total_particles), partphitor(total_particles),
        partvpar(total_particles), partw(total_particles) {}

  void initialize(int total_particles) {
    nptot = total_particles;
    partrad.resize(nptot);
    parttheta.resize(nptot);
    partphitor.resize(nptot);
    partvpar.resize(nptot);
    partw.resize(nptot);
    setZero(); // Initialize all coordinates to zero
  }

  void setValue(const std::vector<double> &partrad_input,
                const std::vector<double> &parttheta_input,
                const std::vector<double> &partphitor_input,
                const std::vector<double> &partvpar_input,
                const std::vector<double> &partw_input) {
    std::copy(partrad_input.begin(), partrad_input.end(), partrad.begin());
    std::copy(parttheta_input.begin(), parttheta_input.end(),
              parttheta.begin());
    std::copy(partphitor_input.begin(), partphitor_input.end(),
              partphitor.begin());
    std::copy(partvpar_input.begin(), partvpar_input.end(), partvpar.begin());
    std::copy(partw_input.begin(), partw_input.end(), partw.begin());
  }

  void setZero() {
    std::fill(partrad.begin(), partrad.end(), 0.0);
    std::fill(parttheta.begin(), parttheta.end(), 0.0);
    std::fill(partphitor.begin(), partphitor.end(), 0.0);
    std::fill(partvpar.begin(), partvpar.end(), 0.0);
    std::fill(partw.begin(), partw.end(), 0.0);
  }

  // Optimized axpy operation using std::transform
  // y=a*x+y, i.e., a*x plus y (PETSC convension)
  // partrad = partrad + dt * dcoords.partrad
  void axpy(const ParticleCoords &dcoords, double dt) {
    std::transform(partrad.begin(), partrad.end(), dcoords.partrad.begin(),
                   partrad.begin(),
                   [dt](double y, double x) { return y + dt * x; });
    std::transform(parttheta.begin(), parttheta.end(),
                   dcoords.parttheta.begin(), parttheta.begin(),
                   [dt](double y, double x) { return y + dt * x; });
    std::transform(partphitor.begin(), partphitor.end(),
                   dcoords.partphitor.begin(), partphitor.begin(),
                   [dt](double y, double x) { return y + dt * x; });
    std::transform(partvpar.begin(), partvpar.end(), dcoords.partvpar.begin(),
                   partvpar.begin(),
                   [dt](double y, double x) { return y + dt * x; });
    std::transform(partw.begin(), partw.end(), dcoords.partw.begin(),
                   partw.begin(),
                   [dt](double y, double x) { return y + dt * x; });
  }

  // partcoords = dt * partcoords
  void ax2a(double dt) {
    std::transform(partrad.begin(), partrad.end(), partrad.begin(),
                   [dt](double x) { return x * dt; });
    std::transform(parttheta.begin(), parttheta.end(), parttheta.begin(),
                   [dt](double x) { return x * dt; });
    std::transform(partphitor.begin(), partphitor.end(), partphitor.begin(),
                   [dt](double x) { return x * dt; });
    std::transform(partvpar.begin(), partvpar.end(), partvpar.begin(),
                   [dt](double x) { return x * dt; });
    std::transform(partw.begin(), partw.end(), partw.begin(),
                   [dt](double x) { return x * dt; });
  }

  void print() const {
    std::cout << "Positions (radial, theta, phi):\n";
    for (int i = 0; i < nptot; ++i) {
      std::cout << "Particle " << i << " -> "
                << "radial: " << partrad[i] << ", "
                << "theta: " << parttheta[i] << ", "
                << "phi: " << partphitor[i] << "\n";
    }
    std::cout << "\nVelocities (vpar, w):\n";
    for (int i = 0; i < nptot; ++i) {
      std::cout << "Particle " << i << " -> "
                << "vpar: " << partvpar[i] << ", "
                << "w: " << partw[i] << "\n";
    }
    std::cout << std::endl;
  }
};

class ParticleSpecies {
private:
  int nptot;
  long long nptot_all = 1000000; // particle # in all PEs
  std::string species_name;
  double mass = 1.0;
  double zcharge = -1.0;
  ParticleCoords coords;
  double Tem = 1.0;
  double nsonN = 1.0;
  int sp_id = 0; // start from 0
  double Cp2g;
  double Cp2g2d1f;
  double Cp2g1d;
  int ideltaf = 1; // 0123

  int ngyro = 1; // number of gyros per particle, default is 1
  int ischeme_motion;
  int v_par0 = 1, v_d = 1, v_mirror = 1, v_ExB = 1, v_Epar = 1, idwdt = 1;

protected:
  int rank = 0;
  int size = 0;

public:
  int idensprof = 1; // density profile type
  std::vector<double> dens_coef1d = {0.49123, 0.298228, 0.198739,
                                     0.0037567601019};
  int iTemprof = 1; // temperature profile type
  std::vector<double> Tem_coef1d = {0.49123, 0.298228, 0.198739,
                                    0.0037567601019};
  std::vector<double> partmu;
  std::vector<double> partfog, partg0;
  double rmin = 0.2, rmax = 0.8;
  // derived parameters
  double rwid = rmax - rmin;
  double rc = (rmax + rmin) * 0.5;
  double vts = std::sqrt(Tem / mass);
  // derived parameters using equ
  double paux_T_transit = 0.0;
  //
  int ntrack = 10;
  std::vector<int> itrack;
  bool irestart = false;

  // Constructor
  ParticleSpecies(size_t nparticles, double m, double q,
                  const std::string &name = "unknown")
      : nptot(nparticles), coords(nparticles), mass(m), zcharge(q),
        species_name(name) {

    MPIManager &mpiManager = MPIManager::getInstance(); // 获取MPIManager的实例
    rank = mpiManager.getRank(); // 获取当前进程的rank
    size = mpiManager.getSize(); // 获取总的进程数

    // Initialize partmu, partfog, partg0
    partmu.resize(nptot);
    partfog.resize(nptot);
    partg0.resize(nptot);
  }

  // Accessors
  int getRank() const { return rank; }
  int getNptot() const { return nptot; }
  int getNgyro() const { return ngyro; }
  long long getNptotAll() const { return nptot_all; }
  double getMass() const { return mass; }
  double getCharge() const { return zcharge; }
  int getSpId() const { return sp_id; }
  ParticleCoords &getCoords() { return coords; }
  const ParticleCoords &getCoords() const { return coords; }
  void setCoords(const ParticleCoords &other) {
    coords.setValue(other.partrad, other.parttheta, other.partphitor,
                    other.partvpar, other.partw);
  }
  double getTem() const { return Tem; }
  double getNsonN() const { return nsonN; }
  double getCp2g() const { return Cp2g; }
  double getCp2g2d1f() const { return Cp2g2d1f; }
  double getCp2g1d() const { return Cp2g1d; }
  void setNptotAll(long long np) { nptot_all = np; }
  void setTem(double t) { Tem = t; }
  void setNsonN(double n) { nsonN = n; }
  void setSpId(int id) { sp_id = id; }
  void setCp2g(double cp2g) { Cp2g = cp2g; }
  void setCp2g2d1f(double cp2g2d1f) { Cp2g2d1f = cp2g2d1f; }
  void setCp2g1d(double cp2g1d) { Cp2g1d = cp2g1d; }
  double getideltaf() const { return ideltaf; }
  void setideltaf(int id) { ideltaf = id; }
  int getIschemeMotion() const { return ischeme_motion; }
  void setIschemeMotion(int scheme) {
    if (scheme >= 1 && scheme <= 2) {
      ischeme_motion = scheme;
    } else {
      std::cerr << "Error: ischeme_motion must be between 0 and 3."
                << std::endl;
      std::abort();
    }
  }
  void setNgyro(int n) {
    if (n > 0) {
      ngyro = n;
    } else {
      std::cerr << "Error: ngyro must be positive." << std::endl;
      std::abort();
    }
  }
  int getvPar0() const { return v_par0; }
  void setvPar0(int v) {
    if (v == 0 || v == 1) {
      v_par0 = v;
    } else {
      std::cerr << "Error: v_par0 must be 0 or 1." << std::endl;
      std::abort();
    }
  }
  int getvD() const { return v_d; }
  void setvD(int v) {
    if (v == 0 || v == 1) {
      v_d = v;
    } else {
      std::cerr << "Error: v_d must be 0 or 1." << std::endl;
      std::abort();
    }
  }
  int getvMirror() const { return v_mirror; } // Mirror force
  void setvMirror(int v) {
    if (v == 0 || v == 1) {
      v_mirror = v; // 0: no mirror, 1: mirror force
    } else {
      std::cerr << "Error: v_mirror must be 0 or 1." << std::endl;
      std::abort();
    }
  }
  int getvExB() const { return v_ExB; } // ExB drift
  void setvExB(int v) {
    if (v == 0 || v == 1) {
      v_ExB = v; // 0: no ExB, 1: ExB drift
    } else {
      std::cerr << "Error: v_ExB must be 0 or 1." << std::endl;
      std::abort();
    }
  }
  int getvEpar() const { return v_Epar; } // Epar
  void setvEpar(int v) {
    if (v == 0 || v == 1) {
      v_Epar = v; // 0: no Epar, 1: Epar drift
    } else {
      std::cerr << "Error: v_Epar must be 0 or 1." << std::endl;
      std::abort();
    }
  }
  int getIdwdt() const { return idwdt; } // 0: no dwdt, 1: dwdt
  void setIdwdt(int dwdtchoice) {
    if (dwdtchoice == 0 || dwdtchoice == 1) {
      idwdt = dwdtchoice; // 0: no dwdt, 1: dwdt
    } else {
      std::cerr << "Error: idwdt must be 0 or 1." << std::endl;
      std::abort();
    }
  }
  //
  double getTem1d(double radius) const {
    // Implement function to get temperature based on radius
    return Tem; // Placeholder
  }
  double getdlnTemdr1d(double radius) const {
    // Implement function to get the derivative of log temperature with respect
    // to radius
    return 0.0; // Placeholder
  }

  double getdens1d(double rad) const {
    double var = 0.0;
    switch (this->idensprof) {
    case 0:
      var = 0.0;
      break;
    case 1:
      var = 1.0;
      break;
    case 2:
      // ITPA EP
      var = dens_coef1d[3] *
            std::exp(-dens_coef1d[2] / dens_coef1d[1] *
                     std::tanh((rad - dens_coef1d[0]) / dens_coef1d[2]));
      break;
    default:
      // 可选：处理未知情况
      std::cerr << "Unknown idensprof = " << this->idensprof << std::endl;
      break;
    }

    return var;
  }
  double getdlndensdr1d(double rad) const {
    double var = 0.0;
    switch (this->idensprof) {
    case 0:
      var = 0.0;
      break;

    case 1:
      var = 0.0;
      break;
    case 2:
      // ITPA EP
      var = -1.0 / dens_coef1d[1] *
            (1.0 -
             std::pow(std::tanh((rad - dens_coef1d[0]) / dens_coef1d[2]), 2));
      break;
    default:
      // 可选：处理未知情况
      std::cerr << "Unknown idensprof = " << this->idensprof << std::endl;
      break;
    }

    return var;
  }

  double getdens1d_can(Equilibrium &equ, double radius, double theta,
                       double vpar, double mu) {
    // Implement function to get canonical density
    return 1.0; // Placeholder
  }
  double get_fsrcsnk(double radius) const {
    // Implement function to get source/sink term
    return 0.0; // Placeholder
  }

  void particle_onestep_euler_0(Equilibrium &equ, ParticleCoords &dxvdt,
                                double dt) {

    particle_dxvpardt_xd0vpar0_euler(
        equ, coords.partrad, coords.parttheta, coords.partvpar, partmu,
        dxvdt.partrad, dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, true);
    particle_add_coords(dxvdt, dt);
  }
  void particle_onestep_rk4_0(Equilibrium &equ, ParticleCoords &dxvdt,
                              ParticleCoords &dxv_sum, ParticleCoords &xv0,
                              double dt) {
    double dthalf = dt / 2.0;
    double dt1o6 = dt / 6.0;
    double dt2o6 = dt * 2.0 / 6.0;

    dxv_sum.setZero();
    xv0 = coords; // Save the initial coordinates
    // Step 1
    particle_dxvpardt_xd0vpar0_euler(
        equ, coords.partrad, coords.parttheta, coords.partvpar, partmu,
        dxvdt.partrad, dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, true);
    dxv_sum.axpy(dxvdt, dt1o6);

    // Step 2
    particle_add_coords(dxvdt, dthalf);
    particle_dxvpardt_xd0vpar0_euler(
        equ, coords.partrad, coords.parttheta, coords.partvpar, partmu,
        dxvdt.partrad, dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, true);
    dxv_sum.axpy(dxvdt, dt2o6);

    // Step 3
    coords = xv0; // Reset to the initial coordinates
    particle_add_coords(dxvdt, dthalf);
    particle_dxvpardt_xd0vpar0_euler(
        equ, coords.partrad, coords.parttheta, coords.partvpar, partmu,
        dxvdt.partrad, dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, true);
    dxv_sum.axpy(dxvdt, dt2o6);

    // Step 4
    coords = xv0;
    particle_add_coords(dxvdt, dt);
    particle_dxvpardt_xd0vpar0_euler(
        equ, coords.partrad, coords.parttheta, coords.partvpar, partmu,
        dxvdt.partrad, dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, true);
    dxv_sum.axpy(dxvdt, dt1o6);

    // Final step
    coords = xv0;
    particle_add_coords(dxv_sum, 1.0);
  }
  // --Part 1 of 1. xd0vpar0; 2. xEB1; 3. x01vpar1
  void particle_dxvpardt_xd0vpar0_euler(
      const Equilibrium &equ, const std::vector<double> &partrad0,
      const std::vector<double> &parttheta0,
      const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      bool isSetZero) {
    // isSetZero: helps d*dt to pass values
    // this->timer.timerClsTic(2);
    if (isSetZero) {
      std::fill(draddt.begin(), draddt.end(), 0.0);
      std::fill(dthetadt.begin(), dthetadt.end(), 0.0);
      std::fill(dphitordt.begin(), dphitordt.end(), 0.0);
      std::fill(dvpardt.begin(), dvpardt.end(), 0.0);
    }

    // Precompute constant values
    const double rhoN = equ.rhoN;
    const double Bref = equ.getBref();

    for (int i = 0; i < nptot; ++i) {
      double ptrad = partrad0[i];
      double pttheta = parttheta0[i];
      double ptvpar = partvpar0[i];
      double ptmu = partmu0[i];
      double ptRR, ptZZ, ptB, ptC1;
      double ptdBdrad, ptdBdthe, ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11,
          ptg12, ptg22, ptjaco3;

      equ.calcBJgij(ptrad, pttheta, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                    ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22);
      ptjaco3 = ptjaco2 * ptRR;

      // 1 dx/dr from vpar
      dthetadt[i] += ptvpar * ptBthe_ct / ptB;
      dphitordt[i] += ptvpar * ptBphi_ct / ptB;
      // 2 dx/dr from vd
      // note: replaced Bmaxis with Bref=1 (20230414)
      ptC1 = mass / zcharge * rhoN * (ptvpar * ptvpar + ptmu * ptB) * Bref /
             (ptB * ptB * ptB);

      draddt[i] += ptC1 * ptFF * ptdBdthe / ptjaco3;
      dthetadt[i] -= ptC1 * ptFF * ptdBdrad / ptjaco3;
      dphitordt[i] +=
          ptC1 * equ.getdpsidr(ptrad) * ptdBdrad * ptg11 / (ptRR * ptRR);
      // 3 dvpar/dt from mirror force
      dvpardt[i] -= ptmu * equ.getdpsidr(ptrad) * ptdBdthe / (ptjaco3 * ptB);
    }
    // this->timer.timerClsToc(2);
  }

  void particle_add_coords(const ParticleCoords &dcoords, const double dt) {
    if (ideltaf == 1 || ideltaf == 2) {
      coords.axpy(dcoords, dt);
    } else if (ideltaf == 3) {
      for (int fic = 0; fic < nptot; ++fic) {
        coords.partrad[fic] += dcoords.partrad[fic] * dt;
        coords.parttheta[fic] += dcoords.parttheta[fic] * dt;
        coords.partphitor[fic] += dcoords.partphitor[fic] * dt;
        coords.partvpar[fic] += dcoords.partvpar[fic] * dt;
        partfog[fic] = partfog[fic] + dcoords.partw[fic] * dt;
      }
    }
  }

  void particle_cls_track_init(int iset_track, double Bax) {
    if (itrack.empty()) {
      itrack.resize(ntrack);
    }
    for (int fic = 0; fic < ntrack; ++fic) {
      itrack[fic] = fic;
    }
    if (iset_track >= 1) {
      for (int i = 0; i < ntrack; ++i) {
        coords.parttheta[i] = (atan(1.0) * 0.0);
        coords.partphitor[i] = 0.0;
      }

      double vparmin_track = -vts * 2;
      double vparmax_track = vts * 2;

      double mumin_track = 0.15 * vts * vts / Bax;
      double mumax_track = vts * vts / Bax;

      if (iset_track == 1) { // passing particles
        for (int i = 0; i < ntrack; ++i) {
          coords.partrad[i] = rmin + rwid * 0.5;
        }

        for (int fic = 0; fic < ntrack; ++fic) {
          coords.partvpar[fic] =
              (vparmin_track +
               fic * (vparmax_track - vparmin_track) / std::max(ntrack - 1, 1));
        }

        for (int i = 0; i < ntrack; ++i) {
          partmu[i] = vts * vts * 0.04 / Bax;
        }
      } else if (iset_track == 2) { // trapped particles
        for (int i = 0; i < ntrack; ++i) {
          coords.partrad[i] = rmin + rwid * 0.8;
        }
        for (int fic = 0; fic < ntrack; ++fic) {
          partmu[fic] = mumin_track + (fic) * (mumax_track - mumin_track) /
                                          std::max(ntrack - 1, 1);
        }
        for (int i = 0; i < ntrack; ++i) {
          coords.partvpar[i] = vts * (-0.2);
        }
      }
    }
  }

  void particle_cls_test(Equilibrium &equ, int irk, int nrun, double dt_o_Ttr,
                         int iset_track) {
    // Working variables
    ParticleCoords dxvdt(nptot), dxv_sum(nptot), xv0(nptot);
    int fic;
    double dt;

    // Calculate dt
    std::cout << "Calculating dt..." << std::endl;
    std::cout << "vts: " << vts << std::endl;
    std::cout << "equ.getqloc_rt(rc, 0.0): " << equ.getqloc_rt(rc, 0.0)
              << std::endl;
    std::cout << "equ.rmaxis: " << equ.rmaxis << std::endl;
    std::cout << "rc: " << rc << std::endl;
    paux_T_transit =
        std::abs(2.0 * M_PI * equ.getqloc_rt(rc, 0.0) * equ.rmaxis / vts);
    dt = dt_o_Ttr * paux_T_transit;
    if (rank == 0) {
      std::cout << " dt set to " << dt << std::endl;
    }

    // Set passing particles
    double Bax = abs(equ.Bmaxis);
    particle_cls_track_init(iset_track, Bax);

    // Initialize coordinates
    dxvdt.setZero();
    dxv_sum.setZero();
    xv0.setZero();

    // Start loop
    if (rank == 0) {
      std::cout << "--Particle push starts--" << std::endl;
    }

    particle_cls_track_record(equ, 0);

    for (fic = 1; fic <= nrun; ++fic) {
      if (fic % static_cast<int>(std::ceil(1.0 * nrun / 6)) == 0 && rank == 0) {
        std::cout << "  step " << fic << "/" << nrun << std::endl;
      }

      if (irk == 0) {
        particle_onestep_euler_0(equ, dxvdt, dt);
      } else if (irk == 1) {
        particle_onestep_rk4_0(equ, dxvdt, dxv_sum, xv0, dt);
      } else {
        if (rank == 0) {
          std::cout << "  ERROR: set irk to 0 or 1  " << std::endl;
        }
        return;
      }

      particle_cls_track_record(equ, fic);
    }
  }

  void particle_cls_record(const Equilibrium &equ, int rank) {

    if (rank == 0) {
    }
    std::ostringstream oss_sp, oss_rank;

    oss_sp << std::setw(3) << std::setfill('0') << sp_id;
    std::string cfile = oss_sp.str();
    oss_rank << std::setw(3) << std::setfill('0') << rank;
    std::string crank = oss_rank.str();

    // 构造文件名
    std::string filename =
        "data_particle_SP" + cfile + "_RANK" + crank + ".txt";

    std::cout << "----Record to " << filename << "----" << std::endl;

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
    }

    for (int fidx = 0; fidx < nptot; ++fidx) {
      double rad = coords.partrad[fidx];
      double theta = coords.parttheta[fidx];
      double r_val = equ.getR(rad, theta);
      double z_val = equ.getZ(rad, theta);

      outfile << std::scientific << std::setprecision(16)
              << coords.partrad[fidx] << " " << coords.parttheta[fidx] << " "
              << coords.partphitor[fidx] << " " << coords.partvpar[fidx] << " "
              << partmu[fidx] << " " << coords.partw[fidx] << " "
              << partfog[fidx] << " " << r_val << " " << z_val << std::endl;
    }

    outfile.close();
  }

  void particle_cls_track_record(Equilibrium &equ, int irun) {
    // C++文件操作
    std::ofstream outfile;

    // 启动计时
    // timer.clsTic(4);

    if (rank == 0) { // 假设rank是全局变量
      int nfile = 102;
      // 使用安全的 std::ostringstream 代替 sprintf
      std::ostringstream oss;
      // oss << std::setw(10) << sp_id; // 格式化 sp_id，宽度为 10
      oss << sp_id;                  // 格式化 sp_id
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
      for (int fic = 0; fic < ntrack; ++fic) {
        int fidx = itrack[fic];
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

    // 结束计时
    // timer.clsToc(4);
  }

  void particle_cls_link2eq(const Equilibrium &equ) {}

  // Print function
  void print() const {
    std::cout << "Particle Species Info:\n";
    std::cout << "Mass: " << mass << ", Charge: " << zcharge << "\n";
    coords.print();
  }
};

class ParticleGroup {
private:
  int nsp = 0;
  std::vector<std::shared_ptr<ParticleSpecies>> speciesList; // 多种粒子种类

public:
  // Add an existing species
  void addSpecies(std::shared_ptr<ParticleSpecies> species) {
    speciesList.push_back(species);
    nsp = speciesList.size();
  }
  // Add a new species by parameters
  void addSpecies(size_t nparticles, double mass, double charge) {
    // Use std::make_shared to construct the ParticleSpecies object
    speciesList.push_back(
        std::make_shared<ParticleSpecies>(nparticles, mass, charge));
    nsp = speciesList.size();
  }

  // Get a species by index
  ParticleSpecies &getSpecies(size_t index) {
    if (index < speciesList.size()) {
      return *speciesList[index];
    }
    throw std::out_of_range(
        "Index out of bounds"); // Throw an exception if index is out of bounds
  }

  int getSpeciesCount() const { return nsp; }
  size_t getTotalParticles() const {
    size_t total = 0;
    for (const auto &species : speciesList) {
      total += species->getNptot();
    }
    return total;
  }

  // 打印所有粒子种类信息
  void print() const {
    std::cout << "Particle Group Info:\n";
    for (const auto &species : speciesList) {
      species->print();
      std::cout << "-----------------------------\n";
    }
  }
};

class Particle {
private:
  int nsp;
  Equilibrium equ;

public:
  ParticleGroup group;

  int ischeme_motion, irandom_gy;
  int v_par0 = 1, v_d = 1, v_mirror = 1, v_ExB = 1, v_Epar = 1, idwdt = 1;
  int imixvar = 0;
  int ibc_particle = 0, iset_zerosumw;
  bool irestart = false;

  double phitormin = 0.0, phitormax, phitorwid;
  int Nphimult = 1; // 1/Nphimult torus
  // double nsonN = 1.0;
  // double Tem = 1.0;

  double vparmaxovt = 5.0;
  double vparmaxovN;
  double mumaxovN2;
  double vtovN;
  std::string profilename;
  int ion_densprof;
  std::vector<double> ion_dens_coef1d = {0.49123, 0.298228, 0.198739,
                                         0.0037567601019};
  // std::vector<double> Tem_coef1d = {0.49123, 0.298228, 0.198739, 0.521298};
  // int iTemprof;
  // Assuming bspline_1d is a class or struct defined elsewhere
  // bspline_1d dens_bsp, Tem_bsp;
  int load_can, ischeme_load;
  int load_seed;
  bool use_random_seed = false; // 是否使用新的随机数生成器
  double Bax;
  double rmin = 0.2, rmax = 0.8;
  double Stot, Vtot, dens_mean, dens_intg;
  double aminor, rmaj, rwid, rmid, rc;
  double rhots;
  double rhotN;

  // std::vector<double> partrad, parttheta, partphitor, partvpar, partmu;
  // std::vector<double> partw, partfog, partg0;

  int ifilter = 0, filter_m0 = -5, filter_m1 = -2, filter_nc = 2;
  std::vector<int> filter_m1d;
  // !--perturbation--
  double pert_rc = 0.5;
  double pert_rwid = 0.125, pert_thetac = 0.1, pert_thetawid = 0.1;
  double pertamp = 4e-12, pert_lr = 0.0;
  int pert_scheme = 3, pert_m0 = -5, pert_m1 = -2, pert_nc = 2, pert_nmin,
      pert_nmax, pert_nstride, pert_lenm;
  std::vector<int> pert_m1d;
  // !--perturbation--
  double ant_rc, ant_rwid, antamp, ant_w;
  int ant_scheme = 3, ant_m0 = -5, ant_m1 = -2, ant_nc = 2, ant_lenm = 4;
  std::vector<int> ant_m1d;
  // !--neoclassical term--
  int neocls;
  // !--collision--
  int icol = 0;
  double nu_colN = 0.0;
  int pullbackN;

  int irec_track = 1, irec_Eparticle = 1;
  int irec_converge;
  bool iuseprv;
  int isrcsnk;
  std::vector<double> radsrcsnkL, radsrcsnkR;
  std::vector<double> coefsrcsnkL, coefsrcsnkR;
  double partwcut;
  //
  // timer_cls timer;
public:
  double mass_bk = 1.0; // background mass
  double getIondens1d(double rad) const {
    // Background Ion density profile function
    double var = 0.0;
    switch (this->ion_densprof) {
    case 0:
      var = 0.0;
      break;
    case 1:
      var = 1.0;
      break;
    case 2:
      // ITPA ion density profile
      var = 1.0 - ion_dens_coef1d[3] *
                      std::exp(-ion_dens_coef1d[2] / ion_dens_coef1d[1] *
                               std::tanh((rad - ion_dens_coef1d[0]) /
                                         ion_dens_coef1d[2]));
      break;
    default:
      // 可选：处理未知情况
      std::cerr << "Unknown Ion densprof = " << this->ion_densprof << std::endl;
      break;
    }

    return var;
  }
  // getters
  int getNsp() const { return nsp; }
  void readInput(const std::string &filepath, const int rank_in, int mpisize) {
    INIReaderExt reader(filepath);
    std::cout << "Reading input file: " << filepath << std::endl;
    std::cout << "Number of MPI processes: " << mpisize << std::endl;
    if (reader.ParseError() < 0) {
      std::cerr << "Error: Unable to load INI file: " << filepath << std::endl;
    }
    nsp = reader.GetInteger("MENG", "nsp", 0);
    int ischeme_motion = reader.GetInteger("MENG", "ischeme_motion", 2);
    use_random_seed = reader.GetBoolean("MENG", "use_random_seed", false);
    ion_densprof = reader.GetInteger("MENG", "ion_densprof", 1);
    ion_dens_coef1d =
        reader.GetRealList("MENG", "ion_dens_coef1d", {0.0, 0.0, 0.0, 0.0});
    //
    if (use_random_seed) {
      auto now = std::chrono::system_clock::now();
      auto time_seed = static_cast<unsigned>(
          std::chrono::duration_cast<std::chrono::microseconds>(
              now.time_since_epoch())
              .count());
      load_seed = static_cast<int>(time_seed + rank_in);
    } else {
      load_seed = 42 + rank_in; // Default seed for reproducibility
    }
    double ref_rmin = -1.0;
    double ref_rmax = -1.0;
    // Iterate over predefined sections
    for (int i = 1; i <= nsp; ++i) {
      std::string section = "Particle_" + std::to_string(i);
      std::cout << "Processing section: " << section << std::endl;

      size_t nptot_all = reader.GetInteger(section, "nptot_all", 0);
      double mass = reader.GetReal(section, "mass", 0.0);
      double zcharge = reader.GetReal(section, "zcharge", 0.0);
      std::string species_name = reader.Get(section, "species_name", "unknown");
      double rmin_input = reader.GetReal(section, "rmin", 0.0);
      double rmax_input = reader.GetReal(section, "rmax", 1.0);
      double Tem = reader.GetReal(section, "Tem", 1.0);
      double nsonN = reader.GetReal(section, "nsonN", 1.0);
      int ideltaf = reader.GetInteger(section, "ideltaf", 1);
      int ngyro = reader.GetInteger(section, "ngyro", 1);
      int idensprof = reader.GetInteger(section, "idensprof", 1);
      int iTemprof = reader.GetInteger(section, "iTemprof", 1);
      std::vector<double> dens_coef1d =
          reader.GetRealList(section, "dens_coef1d", {0.0, 0.0, 0.0, 0.0});
      std::vector<double> Tem_coef1d =
          reader.GetRealList(section, "Tem_coef1d", {0.0, 0.0, 0.0, 0.0});

      // check rmin and rmax for all species are the same
      if (i == 1) {
        ref_rmin = rmin_input;
        ref_rmax = rmax_input;
        this->rmin = rmin_input;
        this->rmax = rmax_input;
      } else {
        if (std::abs(rmin_input - ref_rmin) > 1e-10 ||
            std::abs(rmax_input - ref_rmax) > 1e-10) {
          std::cerr << "Error: Species " << species_name
                    << " has inconsistent rmin/rmax! "
                    << "(rmin = " << rmin_input << ", rmax = " << rmax_input
                    << "), expected (rmin = " << ref_rmin
                    << ", rmax = " << ref_rmax << ")\n";
          std::exit(EXIT_FAILURE);
        }
      }

      if (nptot_all > 0) {
        int nptot = nptot_all / mpisize; // calculate particles per process
        nptot_all = nptot * mpisize;     // update total number of particles
        auto species = std::make_shared<ParticleSpecies>(nptot, mass, zcharge,
                                                         species_name);
        group.addSpecies(species);
        species->rmin = rmin_input;
        species->rmax = rmax_input;
        species->setTem(Tem);
        species->setNsonN(nsonN);
        species->setSpId(i - 1);
        species->setNptotAll(nptot_all);
        species->setideltaf(ideltaf);
        species->setNgyro(ngyro);
        species->setIschemeMotion(ischeme_motion);
        species->idensprof = idensprof;
        species->iTemprof = iTemprof;
        species->dens_coef1d = dens_coef1d;
        species->Tem_coef1d = Tem_coef1d;
        std::cout << "  Species '" << species_name << "' added successfully.\n";
        // Log the parsed values
        std::cout << "  nptot: " << nptot << ", mass: " << mass
                  << ", charge: " << zcharge << std::endl;
      }
    }
    // check rmin and rmax for all species are the same
  }
  void writeProfiles(const int &rank) {
    if (rank == 0) {
      for (int isp = 0; isp < nsp; ++isp) {
        auto &species = group.getSpecies(isp);
        std::ostringstream oss;
        oss << species.getSpId();
        std::string species_id_str = oss.str();
        std::string filename = "data_species_" + species_id_str + "_dens.txt";

        std::ofstream outfile(filename); // default is std::ios::out
        if (!outfile.is_open()) {
          std::cerr << "Error opening file: " << filename << std::endl;
          return;
        }
        for (int irad = 0; irad < 100; ++irad) {
          double rad = irad * 0.01;
          double dens = species.getdens1d(rad);
          double dlndensdr = species.getdlndensdr1d(rad);
          double tem_rad = species.getTem1d(rad);
          outfile << std::scientific << std::setprecision(16) << rad << " "
                  << dens << " " << dlndensdr << " " << tem_rad << std::endl;
        }
      } // isp loop
      std::cout << "Profiles written successfully.\n";
    } // rank == 0
  }

  // Push particles by modifying their coordinates
  void pushParticle(size_t speciesIndex, double dt,
                    const ParticleCoords &velocityChange) {
    ParticleSpecies &species = group.getSpecies(speciesIndex);
    ParticleCoords &coords = species.getCoords();

    coords.axpy(velocityChange, dt); // Update position based on velocity change
  }

  void print() const { group.print(); }

  // Constructor
  Particle() {}
  Particle(Equilibrium &equ_in, int rank, int mpisize_in) : equ(equ_in) {

    std::cout << "Instance Particle @process " << rank << " of " << mpisize_in
              << ".\n";
    // 1. READINPUT
    readInput("input.ini", rank, mpisize_in);
    if (rank == 0) {
      std::cout << "Finish particle readInput...\n"
                << "nsp =" << nsp << std::endl;
      // diagnostic output
      writeProfiles(rank);
    }
    // 2. Summay of particles
    std::cout << "On process: " << rank << ", Total particles of all species: "
              << group.getTotalParticles() << std::endl;
    // 3.
    this->particle_cls_link2eq(equ);
    // 4. =================Load Markers=====================
    for (int i = 0; i < nsp; i++) {
      {
        auto &species = group.getSpecies(i);
        this->particle_cls_set_parms(equ, species, rank);
        std::cout << "Finish particle_cls_set_parms...\n";
        this->particle_cls_loadmarker(equ, species, load_seed);
        // record particle data
        if (rank == 0) {
          std::cout << "Finish particle_cls_loadmarker...\n";
        }
        species.particle_cls_record(equ, rank);
      }
    }
  }

  Particle(Equilibrium &equ_in, ParticleSpecies &sp_in, int rank,
           int mpisize_in)
      : equ(equ_in) {

    std::cout << "Instance Particle @process " << rank << " of " << mpisize_in
              << ".\n";
    // 1. READINPUT
    nsp = 1; // Only one species in this case
    group.addSpecies(std::make_shared<ParticleSpecies>(sp_in));
    if (rank == 0) {
    }
    // 2. Summay of particles
    std::cout << "On process: " << rank << ", Total particles of all species: "
              << group.getTotalParticles() << std::endl;
    // 3.
    this->particle_cls_link2eq(equ);
    // 4. =================Load Markers=====================
    for (int i = 0; i < nsp; i++) {
      {
        auto &species = group.getSpecies(i);
        this->particle_cls_set_parms(equ, species, rank);
        std::cout << "Finish particle_cls_set_parms...\n";
        this->particle_cls_loadmarker(equ, species, load_seed);
      }
    }
  }

  void particle_cls_link2eq(const Equilibrium &equ) {
    Bax = std::abs(equ.Bmaxis); // note sign!
    aminor = equ.rdim / 2;
    rmaj = equ.rmaxis;
    rhotN = equ.rhoN;
  }

  double particle_cls_getf0(double partrad, double parttheta, double partvpar,
                            double partmu, Equilibrium &equ) {
    // Implement function to get f0
    return 1.0; // Placeholder
  }

  void particle_cls_mean_dens(const Equilibrium &equ,
                              const ParticleSpecies &species, double rmin,
                              double rmax, int nrad, int nthe, double &volume,
                              double &int_dens, double &mean_dens,
                              double &Stot) {
    // Working variables
    std::vector<double> xrad(nrad), wrad(nrad);
    std::vector<double> xthe(nthe), wthe(nthe);

    // Initialize outputs
    mean_dens = 0.0;
    int_dens = 0.0;
    volume = 0.0;
    Stot = 0.0;

    // Compute quadrature points and weights
    std::cout << "Computing quadrature points and weights...\n";
    std::cout << "rmin=" << rmin << ", rmax=" << rmax << std::endl;
    std::cout << "nrad=" << nrad << ", nthe=" << nthe << std::endl;
    std::cout << "xrad.size()=" << xrad.size()
              << ", xthe.size()=" << xthe.size() << std::endl;
    std::cout << "wrad.size()=" << wrad.size()
              << ", wthe.size()=" << wthe.size() << std::endl;
    std::cout << " this->phitorwid=" << this->phitorwid << std::endl;
    UtilMath::lgwt(nrad, rmin, rmax, xrad, wrad);
    double twopi = 2.0 * M_PI;
    UtilMath::lgwt(nthe, 0.0, twopi, xthe, wthe);

    // Debug output (assuming rank == 0)
    std::cout << "sum(wrad) = "
              << std::accumulate(wrad.begin(), wrad.end(), 0.0)
              << ", sum(wthe) = "
              << std::accumulate(wthe.begin(), wthe.end(), 0.0) << std::endl;
    std::cout << "mean(xrad) = "
              << std::accumulate(xrad.begin(), xrad.end(), 0.0) / nrad
              << ", mean(xthe) = "
              << std::accumulate(xthe.begin(), xthe.end(), 0.0) / nthe
              << std::endl;

    // Nested loop for computation
    for (int fic = 0; fic < nrad; ++fic) {
      for (int fjc = 0; fjc < nthe; ++fjc) {
        double w12 = wrad[fic] * wthe[fjc];
        double jaco2 = equ.getjaco2(xrad[fic], xthe[fjc]);
        double jaco3 = jaco2 * equ.getR(xrad[fic], xthe[fjc]);

        volume += w12 * jaco3;
        int_dens += w12 * jaco3 * species.getdens1d(xrad[fic]);
        Stot += w12 * jaco2;
      }
    }

    // Final adjustments
    volume = std::abs(volume) * this->phitorwid;
    int_dens = std::abs(int_dens) * this->phitorwid;
    Stot = std::abs(Stot);

    // Compute mean density
    mean_dens = int_dens / volume;
  }

  void particle_cls_set_parms(Equilibrium &equ, ParticleSpecies &species,
                              int rank) {

    // Initialize default values
    ischeme_motion = 2;
    int ngyro = 1;
    irandom_gy = 1;
    iuseprv = false;

    isrcsnk = 0;
    radsrcsnkL = {0.0, 0.1};
    radsrcsnkR = {0.9, 1.0};
    coefsrcsnkL = {0.0, 0.0, 0.0};
    coefsrcsnkR = {0.0, 0.0, 0.0};

    partwcut = 1e-1;

    pullbackN = 0;
    neocls = 0;

    // Set Run Parameters

    // Marker
    int nptot = species.getNptot();
    int nptot_all = species.getNptotAll();
    double mass = species.getMass();
    double zcharge = species.getCharge();
    double nsonN = species.getNsonN();
    double rmin = species.rmin;
    double rmax = species.rmax;
    double Tem = species.getTem();
    int sp_id = species.getSpId();
    int ideltaf = species.getideltaf();

    // Set parameters for all species
    ischeme_motion = ischeme_motion;
    ngyro = species.getNgyro();
    irandom_gy = irandom_gy;
    imixvar = imixvar;
    // Physics
    vparmaxovt = vparmaxovt;
    // Perturbation
    pert_rc = pert_rc;
    pert_rwid = pert_rwid;
    pert_thetac = pert_thetac;
    pert_thetawid = pert_thetawid;
    pertamp = pertamp;
    pert_lr = pert_lr;
    pert_m0 = pert_m0;
    pert_m1 = pert_m1;
    pert_nc = pert_nc;
    pert_nmin = pert_nmin;
    pert_nmax = pert_nmax;
    pert_nstride = pert_nstride;
    pert_scheme = pert_scheme;
    pert_lenm = std::abs(pert_m0 - pert_m1) + 1;
    pert_m1d.resize(pert_lenm);
    int min_pert_m = std::min(pert_m0, pert_m1);
    for (int fic = 0; fic < pert_lenm; ++fic) {
      pert_m1d[fic] = min_pert_m - 1 + (fic + 1);
    }
    if (rank == 0) {
      std::cout << "lenm=" << pert_lenm << ", pert_m1d=";
      for (const auto &val : pert_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    // Antenna
    ant_rc = 0.5;
    ant_rwid = 0.125;
    antamp = 1e-5;
    double ant_wowTAE = 1.0;
    ant_scheme = 3;
    ant_m0 = -11;
    ant_m1 = -10;
    ant_nc = 6;

    ant_lenm = std::abs(ant_m0 - ant_m1) + 1;
    ant_m1d.resize(ant_lenm);
    int min_ant_m = std::min(ant_m0, ant_m1);
    for (int fic = 0; fic < ant_lenm; ++fic) {
      ant_m1d[fic] = min_ant_m - 1 + (fic + 1);
    }
    if (rank == 0) {
      std::cout << "lenm=" << ant_lenm << ", ant_m1d=";
      for (const auto &val : ant_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    // Density and Temperature Profiles
    // idensprof = 1;
    // iTemprof = 1;
    load_can = 0;
    ischeme_load = 1;
    // Assuming dens_coef1d and Tem_coef1d are already populated
    // If not, they should be read from inputs or initialized here
    rwid = rmax - rmin;
    rmid = (rmin + rmax) / 2.0;
    rc = rmid;

    ifilter = ifilter;
    filter_nc = filter_nc;
    filter_m0 = filter_m0;
    filter_m1 = filter_m1;
    int nlen_filter = std::abs(filter_m1 - filter_m0) + 1;
    if (rank == 0) {
      std::cout << "nlen=" << nlen_filter << std::endl;
    }

    std::vector<int> filter_m1d(nlen_filter);
    int min_filter_m = std::min(filter_m0, filter_m1);
    for (int fic = 0; fic < nlen_filter; ++fic) {
      filter_m1d[fic] = min_filter_m + fic;
    }
    if (rank == 0) {
      std::cout << " filter_m1d=";
      for (const auto &val : filter_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    Nphimult = Nphimult;

    irec_track = irec_track;
    irec_Eparticle = irec_Eparticle;
    irec_converge = irec_converge;
    iuseprv = iuseprv;

    ibc_particle = ibc_particle;
    iset_zerosumw = iset_zerosumw;

    isrcsnk = isrcsnk;
    radsrcsnkL = radsrcsnkL;
    radsrcsnkR = radsrcsnkR;
    coefsrcsnkL = coefsrcsnkL;
    coefsrcsnkR = coefsrcsnkR;
    pullbackN = pullbackN;
    neocls = neocls;
    partwcut = partwcut;
    profilename = profilename;

    // Derived variables
    double nu_colstar = nu_colN / std::pow(std::sqrt(rc / rmaj), 3) *
                        std::sqrt(2.0) * equ.getqloc_rt(rc, 0.0) * rmaj *
                        std::sqrt(mass / Tem);
    // Set experimental n,T profiles
    // particle_cls_read_profile1d();
    // Derived parameters
    phitormin = 0.0;
    phitormax = phitormin + 2.0 * M_PI / static_cast<double>(Nphimult);
    phitorwid = phitormax - phitormin;

    vtovN = std::sqrt(Tem / mass); // calculate thermal velocity
    vparmaxovN = vparmaxovt * vtovN;
    std::cout << "vparmaxovN=" << vparmaxovN << std::endl;
    mumaxovN2 = std::pow(vparmaxovN, 2) / (2.0 * Bax);
    std::cout << "mumaxovN2=" << mumaxovN2 << std::endl;
    species.vts = std::sqrt(Tem / mass);
    // --
    // // test lgwt
    // std::vector<double> xrad_test(10), wrad_test(10);
    // UtilMath::lgwt(10, 0.2, 0.8, xrad_test, wrad_test);
    // for (int i = 0; i < 10; i++) {
    //   std::cout << "xrad_test[" << i << "]=" << xrad_test[i] << ",
    //   wrad_test["
    //             << i << "]=" << wrad_test[i] << std::endl;
    // }

    this->particle_cls_mean_dens(equ, species, this->rmin, this->rmax, 100, 200,
                                 this->Vtot, this->dens_intg, this->dens_mean,
                                 this->Stot);
    std::cout << "Vtot=" << Vtot << ", dens_intg=" << dens_intg
              << ", dens_mean=" << dens_mean << ", Stot=" << Stot << std::endl;
    // this->Stot = M_PI * (std::pow(rmax, 2) - std::pow(rmin, 2));
    // this->Vtot = this->Stot * this->rmaj * this->phitorwid;
    // 3d spline 2021/12/08; use NPTOT_ALL 2022/03/10!
    double Cp2g = Vtot / static_cast<double>(nptot_all);
    Cp2g *= nsonN;
    Cp2g *= dens_mean / equ.rhoN / equ.rhoN;

    double Cp2g2d1f = Cp2g; // / phitorwid;
    if (rank == 0) {
      std::cout << "----set Cp2g2d1f=" << Cp2g2d1f << ", Cp2g=" << Cp2g
                << ", phitorwid=" << phitorwid << std::endl;
    }
    double Cp2g1d =
        Vtot / static_cast<double>(nptot_all); // 1d spline 2022/08/18
    Cp2g1d *= nsonN;
    Cp2g1d *= dens_mean / equ.rhoN / equ.rhoN;
    Cp2g1d /= (4.0 * std::pow(M_PI, 2) * rmaj / static_cast<double>(Nphimult));

    species.setCp2g(Cp2g);
    species.setCp2g2d1f(Cp2g2d1f);
    species.setCp2g1d(Cp2g1d);

    if (rank == 0) {
      std::cout << "  particle mean density=" << dens_mean << std::endl;
    }

    // vA^2 = Cpoisson/Campere = 1/betaN
    double wTAE = std::abs(0.5 / rmaj / equ.getqloc_rt(rc, 0.0) *
                           std::sqrt(1.0 / equ.betaN)); // from field_class
    ant_w = ant_wowTAE * wTAE;
    if (rank == 0) {
      std::cout << "ant_w=" << ant_w << ", wTAE=" << wTAE << std::endl;
    }
    for (int fic = 0; fic < ant_lenm; ++fic) {
      ant_m1d[fic] = std::min(ant_m0, ant_m1) - 1 + (fic + 1);
    }
    if (rank == 0) {
      std::cout << "lenm=" << ant_lenm << ", ant_m1d=";
      for (const auto &val : ant_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    // (Rhots only as diagnosis)
    rhots = rhotN * std::sqrt(Tem * mass) / zcharge;
    species.paux_T_transit = std::abs(2.0 * M_PI * equ.getqloc_rt(rc, 0.0) *
                                      equ.rmaxis / species.vts);

    if (rank == 0) {
      // std::cout with appropriate formatting
      std::cout << std::fixed << std::setprecision(5);
      std::cout << "imixvar=" << imixvar << ", ibc_particle=" << ibc_particle
                << ", sp_id=" << sp_id << ", ideltaf=" << ideltaf
                << ", Nphimult=" << Nphimult << ", nptot=" << nptot
                << ", nptot_all=" << nptot_all
                << ", idensprof=" << species.idensprof
                << ", iTemprof=" << species.iTemprof
                << ", ntrack=" << irec_track << ", ifilter=" << ifilter
                << ", pert_scheme=" << pert_scheme << ", v_par0=" << v_par0
                << ", v_d=" << v_d << ", v_mirror=" << v_mirror
                << ", irec_track=" << irec_track
                << ", irec_Eparticle=" << irec_Eparticle
                << ", irec_converge=" << irec_converge
                << ", filter_nc=" << filter_nc
                << ", ischeme_motion=" << ischeme_motion << ", ngyro=" << ngyro
                << ", irandom_gy=" << irandom_gy
                << ", ibc_particle=" << ibc_particle
                << ", load_can=" << load_can
                << ", ischeme_load=" << ischeme_load
                << ", pert_nmin=" << pert_nmin << ", pert_nmax=" << pert_nmax
                << ", pert_nstride=" << pert_nstride
                << ", pullbackN=" << pullbackN << ", neocls=" << neocls
                << ", icol=" << icol << ", iset_zerosumw=" << iset_zerosumw
                << std::endl;

      std::cout << "phitormin=" << phitormin << ", phitormax=" << phitormax
                << ", phitorwid=" << phitorwid << ", nsonN=" << nsonN
                << ", Tem=" << Tem << ", mass=" << mass
                << ", zcharge=" << zcharge << ", vparmaxovt=" << vparmaxovt
                << ", vparmaxovN=" << vparmaxovN << ", vtovN=" << vtovN
                << ", Bax=" << Bax << ", rmin=" << rmin << ", rmax=" << rmax
                << ", Stot=" << Stot << ", Vtot=" << Vtot
                << ", aminor=" << 0.0 // Assuming aminor is defined elsewhere
                << ", rmaj=" << rmaj << ", rwid=" << rwid << ", rmid=" << rmid
                << ", rc=" << rc << ", rhots=" << rhots << ", rhotN=" << rhotN
                << ", Cp2g=" << Cp2g << ", vts=" << species.vts
                << ", paux_T_tr=" << species.paux_T_transit
                << ", pert_rc=" << pert_rc << ", pert_rwid=" << pert_rwid
                << ", pertamp=" << pertamp << ", dens_mean=" << dens_mean
                << ", Cp2g1d=" << Cp2g1d << ", nu_colN=" << nu_colN
                << ", nu_colstar=" << nu_colstar
                << ", pert_thetac=" << pert_thetac
                << ", pert_thetawid=" << pert_thetawid
                << ", pert_lr=" << pert_lr << ", partwcut=" << partwcut
                << std::endl;

      std::cout << "isrcsnk=" << isrcsnk << std::endl;
      std::cout << "radsrcsnkL=";
      for (const auto &val : radsrcsnkL) {
        std::cout << val << " ";
      }
      std::cout << ", coefsrcsnkL=";
      for (const auto &val : coefsrcsnkL) {
        std::cout << val << " ";
      }
      std::cout << std::endl;

      std::cout << "radsrcsnkR=";
      for (const auto &val : radsrcsnkR) {
        std::cout << val << " ";
      }
      std::cout << ", coefsrcsnkR=";
      for (const auto &val : coefsrcsnkR) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
      std::cout << "iuseprv=" << iuseprv << std::endl;
    }
  }

  void
  particle_setvalue2sp_from_coords(Particle &pt,
                                   const std::vector<ParticleCoords> &xv0) {
    int nsp_tmp = xv0.size();
    if (nsp_tmp != pt.getNsp()) {
      std::cerr << "Error: wrong nsp in particle_setvalue2sp_from_coords"
                << std::endl;
    }
    for (int fsc = 0; fsc < nsp_tmp; ++fsc) {
      ParticleSpecies &species = pt.group.getSpecies(fsc);
      species.setCoords(xv0[fsc]);
    }
  }

  void particle_coords_cls_axpy2sp(std::vector<ParticleCoords> &dxv_sum,
                                   const std::vector<ParticleCoords> &dxvdt,
                                   const double dt) {
    int nsp_tmp = dxv_sum.size();
    for (int fsc = 0; fsc < nsp_tmp; ++fsc) {
      dxv_sum[fsc].axpy(dxvdt[fsc], dt);
    }
  }

  void particle_add_coords2sp(Particle &pt,
                              const std::vector<ParticleCoords> &dxvdt,
                              const double dt) {
    int nsp_tmp = dxvdt.size();
    if (nsp_tmp != pt.getNsp()) {
      std::cerr << "Error: wrong nsp in particle_add_coords2sp" << std::endl;
    }
    for (int fsc = 0; fsc < nsp_tmp; ++fsc) {
      ParticleSpecies &species = pt.group.getSpecies(fsc);
      species.particle_add_coords(dxvdt[fsc], dt);
    }
  }

  std::mt19937 gen;

  void init_random_generator(int seed = 42) { gen.seed(seed); }

  double rand_uniform_01() {
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
  }

  double rand_normal(double mu, double sigma) {
    std::normal_distribution<double> dist(mu, sigma);
    return dist(gen);
  }

  void random_vector_uniform(std::vector<double> &vec) {
    for (auto &v : vec) {
      v = rand_uniform_01();
    }
  }

  void particle_cls_loadmarker(Equilibrium &equ, ParticleSpecies &species,
                               int seed = 42) {
    std::cout << "--------LOAD MARKER (Seed = " << seed << ")--------"
              << std::endl;

    init_random_generator(seed);

    auto &coords = species.getCoords();
    auto &partrad = coords.partrad;
    auto &parttheta = coords.parttheta;
    auto &partphitor = coords.partphitor;
    auto &partvpar = coords.partvpar;
    auto &partw = coords.partw;

    auto &partmu = species.partmu;
    auto &partfog = species.partfog;
    auto &partg0 = species.partg0;
    int nptot = species.getNptot();
    int ideltaf = species.getideltaf();

    const double mass = species.getMass();

    // Init Random
    random_vector_uniform(partrad);
    random_vector_uniform(parttheta);
    random_vector_uniform(partphitor);

    double c0rnd, c1rnd;
    double Tpar, Tperp, zwidprod, BB;
    c0rnd = rmin * rmin;
    c1rnd = rmax * rmax - rmin * rmin;

    for (int fic = 0; fic < nptot; ++fic) {

      // 1. Rad
      if (ischeme_load == 1) {
        partrad[fic] = sqrt(partrad[fic] * c1rnd + c0rnd);
      } else if (ischeme_load == 2) {
        partrad[fic] = partrad[fic] * rwid + rmin;
      }

      // 2.3 Theta, Phitor
      parttheta[fic] = parttheta[fic] * 2.0 * M_PI;
      partphitor[fic] = partphitor[fic] * phitorwid + phitormin;

      // 4. Vpar : truncated normal
      Tpar = species.getTem1d(partrad[fic]);
      double stddev_vpar = sqrt(Tpar / mass * 0.5);
      double vpar;
      do {
        vpar = rand_normal(0.0, stddev_vpar);
      } while (std::abs(vpar) >= vparmaxovN);
      partvpar[fic] = vpar;

      // 5. Mu: truncated exponential
      zwidprod = phitorwid * rwid * 2.0 * M_PI;
      Tperp = species.getTem1d(partrad[fic]);
      BB = equ.getB(partrad[fic], parttheta[fic]);
      double mu;
      do {
        double rnd = rand_uniform_01();
        mu = -Tperp / mass * std::log(rnd) / BB / 2.0;
      } while (mu >= mumaxovN2);
      partmu[fic] = mu;

      // 1b. full weight
      if (load_can == 0) {
        partfog[fic] = species.getdens1d(partrad[fic]) / dens_mean;
      } else if (load_can == 1) {
        partfog[fic] = species.getdens1d_can(equ, partrad[fic], parttheta[fic],
                                             partvpar[fic], partmu[fic]) /
                       dens_mean;
      }

      double jaco = equ.getjaco3(partrad[fic], parttheta[fic]);
      if (ischeme_load == 1) {
        partfog[fic] *= jaco * zwidprod / Vtot * rmid / partrad[fic];
      } else if (ischeme_load == 2) {
        partfog[fic] *= jaco * zwidprod / Vtot;
      } else {
        std::cerr << "Error: wrong ISCHEME_LOAD value; set to 1" << std::endl;
      }

      // 1a. weight
      if (ideltaf > 0) {
        partw[fic] = 0.0;
      } else if (ideltaf == 0) {
        partw[fic] = partfog[fic];
      } else {
        std::cerr << "Error: wrong IDELTAF value; set to 1" << std::endl;
      }

      if (ideltaf == 3) {
        partg0[fic] = particle_cls_getf0(partrad[fic], parttheta[fic],
                                         partvpar[fic], partmu[fic], equ) /
                      partfog[fic];
      }
      // print particle info
      // std::cout << "partrad[" << fic << "]: " << partrad[fic] << std::endl;
      // std::cout << "parttheta[" << fic << "]: " << parttheta[fic] <<
      // std::endl; std::cout << "partphitor[" << fic << "]: " <<
      // partphitor[fic]
      //           << std::endl;
      // std::cout << "partvpar[" << fic << "]: " << partvpar[fic] <<
      // std::endl; std::cout << "partmu[" << fic << "]: " << partmu[fic] <<
      // std::endl; std::cout << "partfog[" << fic << "]: " << partfog[fic] <<
      // std::endl; std::cout << "partg0[" << fic << "]: " << partg0[fic] <<
      // std::endl; std::cout << "partw[" << fic << "]: " << partw[fic] <<
      // std::endl; std::cout << "---------------------------" << std::endl;
    }
  }
};
#endif // PARTICLE_H
