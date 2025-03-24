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

#include "../inih/INIReader.h"
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


void setbuffcoef(const std::array<double, 2>& rminbuff, std::array<double, 4>& coefminbuff);

class FieldCls {
private:
  Equilibrium equ;
  int rank = 0;
  int size = 0;

public:
  // Control
  bool nl_debug = false;
  bool irestart = false;
  int imixvar = 1;
  int isolver;
  int imodel;
  int ischeme_Ohm;

  // Diagnostic
  int irandom_gy;
  int ndiag_field_h5 = 1;

  // Geometry
  double rmaj;

  // Grid
  int Nphimult = 2;
  int bc_arr[3] = {1, 0, 0};
  int nrad = 12, nthe = 32, nphi = 8;
  double radmin = 0.1, themin = 0.0, phimin = 0.0;
  double radmax = 0.9, themax = twopi, phimax = pi;
  double radwid, thewid, phiwid;

  int bcrad;
  int nradfem, nthefem, nphifem;
  int ntotfem, ntotdof;

  // Filter
  int ifilter = 1, filter_m0 = -5, filter_m1 = -2, filter_nc = 2;
  std::vector<int> filter_m1d;

  // Buffer region
  int ibuff; // 0: no buffer; 1: rmax only; 2: rmin only; 3: both rmin, rmax
  double rminbuff[2], rmaxbuff[2];
  double coefminbuff[4], coefmaxbuff[4];

  // Physics parameters
  // SplineCubic3DExtCls spc;
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
  void readInput(const std::string &inputFile);
  void field_cls_init(const Equilibrium &equ);
  void field_cls_set_parms();
  double field_cls_get_fbuff(double rad, const std::array<double, 2> &rbuff,
                             const std::array<double, 4> &coef,
                             bool left) const;
  void field_cls_link2eq(const Equilibrium &equ);
  void field_cls_g2p2d1f_general(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d, const std::vector<double> &ptrad1d,
      const std::vector<double> &ptthe1d, const std::vector<double> &ptphi1d,
      std::vector<double> &ptf1d, const std::vector<int> &idiff, int ngyro,
      const std::vector<double> &rho1);

  void field_cls_g2p2d1f_grad(
      const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
      const std::vector<int> &ntor1d, const std::vector<double> &ptrad1d,
      const std::vector<double> &ptthe1d, const std::vector<double> &ptphi1d,
      std::vector<double> &ptf1d100, std::vector<double> &ptf1d010,
      std::vector<double> &ptf1d001, int ngyro,
      const std::vector<double> &rho1);
  void field_cls_final();

  // default Constructor
  // Constructor definition
  FieldCls() {
    // Initialize MPI rank and shared memory parameters
    MPIManager &mpiManager = MPIManager::getInstance();
    rank = mpiManager.getRank();
    size = mpiManager.getSize();

    int bc_arr[3] = {1, 0, 0};
    int nnode_arr[3] = {nrad, nthe, nphi};
    double zmin_arr[3] = {radmin, themin, phimin};
    double zmax_arr[3] = {radmax, themax, phimax};

    std::cout << "Initializing SplineCubic2d1f object..." << std::endl;

    spc = SplineCubic2d1f(0, bc_arr, nnode_arr, zmin_arr, zmax_arr, nl_debug);

    // Initialize field_cls with the equ object
    field_cls_init(equ);
  }
};

class Field {
private:
  double fieldStrength;

public:
  Field(double strength = 1.0);
  double getStrength() const;
  void setStrength(double strength);
  double calculateForce(double charge) const;
};

#endif // FIELD_H
