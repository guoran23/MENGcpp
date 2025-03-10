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

#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "../inih/INIReader.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <array>

class Equilibrium {
private:
  double pressure;
  double temperature;
  double Btor;
  double aminor;
  double qbar_axis;
  double qbar_edge;
  double Bref = 1.0;

  void validateInputs() {
    if (iequmodel <= 0) {
      throw std::runtime_error("Error: Invalid 'iequmodel' value.");
    }
    if (c1adhoc < 0.0) {
      throw std::runtime_error("Error: c1adhoc must be non-negative.");
    }
    if (c2adhoc < 0.0) {
      throw std::runtime_error("Error: c2adhoc must be non-negative.");
    }
  }

public:
  Equilibrium(double p, double t);
  double getPressure() const;
  double getTemperature() const;
  double getBref() const { return Bref; }
  void setBref(double B) { Bref = B; }

  void adjustEquilibrium(double dPressure, double dTemperature);
  void displayEquilibrium() const;
  // ---------------------------------------------------------
  // Constructor
  Equilibrium()
      : iequmodel(2), c1adhoc(1.71), c2adhoc(0.16), rmaxis_adhoc(10.0),
        Bmaxis_adhoc(2.0), set_zerof(false), fname("g031213.00003"),
        profilename("profile.h5") {}

  // Member function
  void readInput(const std::string &inputFile) {
    INIReader reader(inputFile);

    // Check if reading was successful
    if (reader.ParseError() < 0) {
      throw std::runtime_error("Error: Unable to open or parse input file '" +
                               inputFile + "'");
    }
    std::cout << "eq readInput: "
              << "Equilibrium" << std::endl;
    // Read key-value pairs
    iequmodel = reader.GetInteger("Equilibrium", "iequmodel", 2);
    c1adhoc = reader.GetReal("Equilibrium", "c1adhoc", 0.0);
    c2adhoc = reader.GetReal("Equilibrium", "c2adhoc", 0.0);
    rmaxis_adhoc = reader.GetReal("Equilibrium", "rmaxis_adhoc", 0.0);
    Bmaxis_adhoc = reader.GetReal("Equilibrium", "Bmaxis_adhoc", 0.0);
    set_zerof = reader.GetBoolean("Equilibrium", "set_zerof", false);
    fname = reader.Get("Equilibrium", "fname", "");
    profilename = reader.Get("Equilibrium", "profilename", "");
    rhoN = reader.GetReal("Equilibrium", "rhoN", 0.0);
    Bref = reader.GetReal("Equilibrium", "Bref", 1.0);

    rmaxis = rmaxis_adhoc;
    Bmaxis = Bmaxis_adhoc;

    // Optional: Add validation and other logic
    validateInputs();
  }

  void printState() const {
    std::cout << "iequmodel: " << iequmodel << "\n"
              << "c1adhoc: " << c1adhoc << "\n"
              << "c2adhoc: " << c2adhoc << "\n"
              << "rmaxis_adhoc: " << rmaxis_adhoc << "\n"
              << "aminor: " << aminor << "\n"
              << "Bmaxis_adhoc: " << Bmaxis_adhoc << "\n"
              << "set_zerof: " << (set_zerof ? "true" : "false") << "\n"
              << "fname: " << fname << "\n"
              << "profilename: " << profilename << "\n"
              << "Bref: " << Bref << "\n"
              << "rhoN: " << rhoN << "\n";
  }
  // ---------------------------------------------------------
  double get_qbar(double rho) const {
    double qbar = c1adhoc + c2adhoc * rho * rho;
    return qbar;
  }
  // 1.0 controls
  int iequmodel;
  bool set_zerof; // set fpol=0 manually
  double c1adhoc, c2adhoc, rmaxis_adhoc, Bmaxis_adhoc;

  // 1.1 eqdsk
  std::string fname;
  std::string profilename;
  std::string case_[6];
  int idum, nw, nh, nbbbs, limitr;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr,
      current, xdum;
  std::vector<double> fpol, pres, ffprim, pprime, qpsi, rbbbs, zbbbs, rlim,
      zlim;
  std::vector<std::vector<double>> psirz;

  // 1.2 derived variables
  double Bmaxis; // added 20200825. B on axis. bcentr--rcentr; Bmaxis--rmaxis
  double signBphi;
  double dr, dz, dpsi, psiwid;
  std::vector<double> r1d, z1d, psi1d;
  int isetpsi;
  double fBphi;

  // 1.3 spline class
  bool extrap;
  int r_order_bsp, z_order_bsp, phi_order_bsp, psi_order_bsp;
  int r_order_bsp3d, z_order_bsp3d, phi_order_bsp3d;

  // // bspline coefficients in (R,Z)
  // class Bspline2D radrz_bsp;
  // class Bspline1D F1d_bsp;
  // class Bspline1D q1d_bsp;

  // // bspline in (rad,the)
  // class Bspline1D rbbbs_bsp, zbbbs_bsp;
  // class Bspline2D rrt_bsp, zrt_bsp, Brt_bsp;
  // class Bspline2D Brad_co_bsp, Bthe_co_bsp;
  // class Bspline2D curlbrad_ct_bsp, curlrbthe_ct_bsp, curlbphi_ct_bsp;

  int nrad_sp = 32, nthe_sp = 64;
  double radmax_sp = 0.7, radmin_sp = 0.3;
  double drad_sp, dthe_sp;
  std::vector<double> rad1d_sp, the1d_sp;

  // // bspline in (rad,the,phi)
  // class Bspline3D eta3d_bsp;
  // class Bspline3D the3d_bsp;

  double nphieachfa = 4; // note: can be non-integer; set to 4 for cubic
  int nphi = 7;          // toroidal grid #
  int nphiMult = 2;      // torus wedge factor
  double phimaxfull = 2 * M_PI, phimin = 0.0;
  double phimax, phiwid, dphi;
  double phimin_fa, phimax_fa, dphi_fa;
  int nphifa_sp = 8; // toroidal grid # of 3d Spline
  int nintg_fa = 10;
  std::vector<double> phifa1d_sp;

  // 3.1 physics variables
  int isource_ref = 1;
  double rhoN, betaN, nref, Tref;

  // methods
  void equil_cls_init();
  void equil_cls_readinput();
  void equil_cls_calc_derived();

  // (R,Z) functions
  double getradRZ(double R, double Z) const;
  double gettheRZ(double R, double Z) const;

  // (rad,the) functions
  double getBdirect(double rad, double the) const;
  double getBrad_co_direct(double rad, double the) const;
  double getBthe_co_direct(double rad, double the) const;
  double getB(double rad, double the) const;
  double getdBdrad(double rad, double the) const;
  double getdBdthe(double rad, double the) const;
  double getBthe_ct(double rad, double the) const;
  double getBphi_ct(double rad, double the) const;
  double getcurlbrad_ct(double rad, double the) const;
  double getcurlbthe_ct(double rad, double the) const;
  double getcurlbphi_ct(double rad, double the) const;

  // for Rho* terms of GC motion
  double getBrad_co(double rad, double the) const;
  double getBthe_co(double rad, double the) const;
  double getBphi_co(double rad, double the) const;
  double getdBrad_codrad(double rad, double the) const;
  double getdBrad_codthe(double rad, double the) const;
  double getdBthe_codrad(double rad, double the) const;
  double getdBthe_codthe(double rad, double the) const;
  double geth_adhoc(double rad, double the, double rmaxis) const;

  void calc_curlb_ct_direct(double rad, double the);

  void calcBJgij(double rad, double the0, double &RR, double &ZZ, double &BB,
                 double &dBdrad, double &dBdthe, double &Bthe_ct,
                 double &Bphi_ct, double &jaco2, double &FF, double &g11,
                 double &g12, double &g22) const;

  double getgij_fa(double rad, double the, double phi);
  double getgij_rtp(double rad, double the, double phi);

  double getR(double rad, double the) const;
  double getZ(double rad, double the) const;
  double getdRdrad(double rad, double the) const;
  double getdZdrad(double rad, double the) const;
  double getdRdthe(double rad, double the) const;
  double getdZdthe(double rad, double the) const;
  void calcdrtdRZ(double rad, double the, double &drdR, double &drdZ,
                  double &dtdR, double &dtdZ, double &jaco2) const;

  double getjaco2(double rad, double the) const;
  double getjaco3(double rad, double the) const;

  double getqloc_rt(double rad, double the) const;
  double getdqdrloc_rt(double rad, double the) const;
  double getfpol_r(double rad) const;
  double getdfpoldr_r(double rad) const;
  double getdpsidr(double rad) const;
  double getpsi(double rad) const;

  double geteta(double rad, double the, double phi);
  double getdetadrad(double rad, double the, double phi);
  double getdetadthe(double rad, double the, double phi);
  double getdetadphi(double rad, double the, double phi);

  double getthe3d(double rad, double the, double phi);
  double getdthe3ddrad(double rad, double the, double phi);
  double getdthe3ddthe(double rad, double the, double phi);
  double getdthe3ddphi(double rad, double the, double phi);

  double getfeq();

  double getcos_straight_adhoc(double rad, double the0, double rmaxis) const;
  double getsin_straight_adhoc(double rad, double the0, double rmaxis) const;

  void equil_cls_showinfo(int rank);
  void writeWallRZ() const;

  // Initialization function
  void equil_cls_init_adhoc(int rank) {
    if (rank == 0) {
      std::cout << "----initialize iequmodel==2----" << std::endl;
    }

    // Set variables
    rmaxis = rmaxis_adhoc;
    Bmaxis = Bmaxis_adhoc;
    zmaxis = 0.0;

    // Allocate arrays for diagnosis equvar2d, 3d
    rad1d_sp.resize(nrad_sp);
    the1d_sp.resize(nthe_sp);

    // Set radial and angular arrays
    drad_sp = (radmax_sp - radmin_sp) / (nrad_sp - 1);
    double dthe_sp = 2 * M_PI / (nthe_sp - 1);

    for (int fic = 0; fic < nrad_sp; ++fic) {
      rad1d_sp[fic] = radmin_sp + fic * drad_sp;
    }

    for (int fic = 0; fic < nthe_sp; ++fic) {
      the1d_sp[fic] = fic * dthe_sp;
    }

    // Construct fa3d
    phiwid = phimaxfull - phimin;
    phimax = phiwid / nphiMult;
    dphi = phiwid / nphi;
    phimin_fa = -dphi * nphieachfa / 2;
    phimax_fa = dphi * nphieachfa / 2;
    dphi_fa = dphi * nphieachfa / (nphifa_sp - 1);

    phifa1d_sp.resize(nphifa_sp);

    for (int fic = 0; fic < nphifa_sp; ++fic) {
      phifa1d_sp[fic] = phimin_fa + fic * dphi_fa;
    }
  }
};

#endif // EQUILIBRIUM_H
