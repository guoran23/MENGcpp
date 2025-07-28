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

#ifndef FIELD_H
#define FIELD_H

#include "../inih/INIReaderExt.h"
#include "Equilibrium.h"
#include "MPIManager.h"
#include "Particle.h"
#include "SplineCubic2d1f.h"
#include "SplineCubicNd.h"
#include "util_math.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <vector>
// #include <petsc.h>

void setbuffcoef(const std::array<double, 2> &rminbuff,
                 std::array<double, 4> &coefminbuff);

class FieldCls {
private:
  int rank = 0;
  int size = 0;

public:
  // Control
  bool nl_debug = false;
  bool irestart = false;
  int imixvar = 0;
  int isolver;
  int imodel;
  int ischeme_Ohm;

  // Diagnostic
  int irandom_gy;
  int ndiag_field_h5 = 1;

  // Geometry
  double rmaj;

  // Grid
  int Nphimult = 1;
  int bc_arr[3] = {1, 0, 0};
  int nrad = 12, nthe = 32, nphi = 8;
  double radmin = 0.1, themin = 0.0, phimin = 0.0;
  double radmax = 0.9, themax = twopi, phimax = pi;
  double radwid, thewid, phiwid;

  int bcrad;
  int nradfem, nthefem, nphifem;
  int ntotfem, ntotdof;
  int ntotfem2d1f;
  int ntot12fem;

  // Filter
  int ifilter = 1, filter_m0 = -5, filter_m1 = -2, filter_nc = 2;
  std::vector<int> filter_m1d;

  // Buffer region
  int ibuff; // 0: no buffer; 1: rmax only; 2: rmin only; 3: both rmin, rmax
  double rminbuff[2], rmaxbuff[2];
  double coefminbuff[4], coefmaxbuff[4];

  // Physics parameters
  SplineCubic2d1f spc;
  double TA, wTAE, Cpoisson, Campere;

  // Implicit solver
  int icorr;

  // 2D1F data: ntor_min, ntor_max, ntor_stride from input file
  int ntor_min, ntor_max, ntor_stride, lenntor;
  std::vector<int> ntor1d;

  // Iterative field solver
  int niter_sv;

  // Matrix construction schemes
  int imatscheme;

  // Timer

  // Methods
  // Getters
  int getNtotfem2d1f() const { return ntotfem2d1f; }
  void readInput(const std::string &inputFile);
  void field_cls_init(const Equilibrium &equ);
  void field_cls_set_parms();
  double field_cls_get_fbuff(double rad, const std::array<double, 2> &rbuff,
                             const std::array<double, 4> &coef,
                             bool left) const;
  void field_cls_link2eq(const Equilibrium &equ);
  void field_cls_g2p2d1f_general(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amplitude_arr,
      const std::vector<double> &ptrad1d, const std::vector<double> &ptthe1d,
      const std::vector<double> &ptphi1d, std::vector<double> &ptf1d,
      const std::array<int, 3> &idiff, int ngyro,
      const std::vector<double> &rho1);
  void field_cls_g2p2d1f_general(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amplitude_arr,
      const double ptrad1d, const double ptthe1d, const double ptphi1d,
      double &ptf1d, const std::array<int, 3> &idiff, int ngyro, double rho1);

  void field_cls_g2p2d1f_grad(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amplitude_arr,
      const double ptrad1d, const double ptthe1d, const double ptphi1d,
      double &ptf1d100, double &ptf1d010, double &ptf1d001, int ngyro,
      double rho1) const;
  void field_cls_g2p2d1f_grad(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amplitude_arr,
      const std::vector<double> &ptrad1d, const std::vector<double> &ptthe1d,
      const std::vector<double> &ptphi1d, std::vector<double> &ptf1d100,
      std::vector<double> &ptf1d010, std::vector<double> &ptf1d001, int ngyro,
      const std::vector<double> &rho1) const;
  void field_cls_g2p2d1f_grad_complex(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amplitude_arr,
      const std::vector<double> &ptrad1d, const std::vector<double> &ptthe1d,
      const std::vector<double> &ptphi1d,
      std::vector<std::complex<double>> &ptf1d100_c,
      std::vector<std::complex<double>> &ptf1d010_c,
      std::vector<std::complex<double>> &ptf1d001_c, int ngyro,
      const std::vector<double> &rho1) const;
  void field_cls_g2p2d1f_grad_complex(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d,
      const std::vector<std::complex<double>> &amplitude_arr,
      const double ptrad1d, const double ptthe1d, const double ptphi1d,
      std::vector<std::complex<double>> &ptf1d100_c,
      std::vector<std::complex<double>> &ptf1d010_c,
      std::vector<std::complex<double>> &ptf1d001_c, int ngyro,
      double rho1) const;

  void field_cls_final();

  // default Constructor
  // Constructor definition
  FieldCls() {
    // Initialize MPI rank and shared memory parameters
    MPIManager &mpiManager = MPIManager::getInstance();
    rank = mpiManager.getRank();
    size = mpiManager.getSize();
    std::cout << "FieldCls default constructor called"
              << std::endl; // For debugging
  }
};

#endif // FIELD_H
