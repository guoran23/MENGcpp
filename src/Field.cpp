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

#include "Field.h"

void setbuffcoef(const std::array<double, 2> &rminbuff,
                 std::array<double, 4> &coefminbuff) {
  // Compute coefficients based on rminbuff values
  coefminbuff[3] = -2.0 / std::pow(rminbuff[0] - rminbuff[1], 3);
  coefminbuff[2] = -1.5 * (rminbuff[0] + rminbuff[1]) * coefminbuff[3];
  coefminbuff[1] = 3.0 * rminbuff[0] * rminbuff[1] * coefminbuff[3];
  coefminbuff[0] = -coefminbuff[1] * rminbuff[1] -
                   coefminbuff[2] * std::pow(rminbuff[1], 2) -
                   coefminbuff[3] * std::pow(rminbuff[1], 3);
}

// Definition of field_cls_init
void FieldCls::field_cls_init(const Equilibrium &equ) {
  if (rank == 0) {
    std::cout << "======== FieldCls Initialization ========" << std::endl;
  }

  // Set class parameters
  field_cls_set_parms();

  // Link the field class to the equilibrium class
  field_cls_link2eq(equ);
  int bc_arr[3] = {1, 0, 0};
  int nnode_arr[3] = {nrad, nthe, nphi};
  double zmin_arr[3] = {radmin, themin, phimin};
  double zmax_arr[3] = {radmax, themax, phimax};

  std::cout << "Initializing SplineCubic2d1f object..." << std::endl;

  spc = SplineCubic2d1f(0, bc_arr, nnode_arr, zmin_arr, zmax_arr, nl_debug);
  // Print memory usage details
  if (rank == 0) {
    std::cout << "Field Class size: " << sizeof(*this) << " Byte " << std::endl;
  }
}

void FieldCls::field_cls_link2eq(const Equilibrium &equ) {
  // Compute physical parameters based on equilibrium class
  Cpoisson = 1.0 / (equ.rhoN * equ.rhoN);
  Campere =
      Cpoisson * equ.betaN; // Note: No division by meomi due to N normalization

  rmaj = equ.rmaxis;
  wTAE = 0.5 / (rmaj * std::abs(equ.getqloc_rt(0.5, 0.0))) *
         std::sqrt(Cpoisson / Campere);
  TA = (2 * M_PI) / wTAE;

  // Print the computed values if running on the main rank
  if (rank == 0) {
    std::cout << "Cpoisson=" << Cpoisson << ", Campere=" << Campere
              << ", wTAE=" << wTAE << ", TA=" << TA << std::endl;
  }
}

void FieldCls::field_cls_set_parms() {
  // Declare local variables
  bool nl_debug, irestart;
  int Nphimult, nrad, nthe, nphi, ndiag_field_h5, imixvar;
  double radmin, radmax;
  int ifilter, filter_m0, filter_m1, filter_nc;
  int icorr, irandom_gy, isolver, imodel, ischeme_Ohm;
  int ntor_min, ntor_max, ntor_stride, bcrad, niter_sv;
  int ibuff; // 0: no buffer, 1: rmax only, 2: rmin only, 3: both rmin, rmax
  std::array<double, 2> rminbuff = {0.0, 0.0}, rmaxbuff = {1.0, 1.0}, rbufftmp;
  int imatscheme;

  // Set default values
  icorr = 1;
  irandom_gy = 1;
  isolver = 1;
  imodel = 0;
  ischeme_Ohm = 1;

  ntor_min = 1;
  ntor_max = 1;
  ntor_stride = 1;
  bcrad = 0;

  ibuff = 0;
  niter_sv = 0;
  imixvar = 0;
  imatscheme = 1;

  // Read parameters from file
  // Assign values to the class members
  readInput("input.ini");

  //
  //
  // const int bc_arr[3] = {bcrad, 0, 0};
  // const int nnode_arr[3] = {nrad, nthe, nphi};
  // const double zmin_arr[3] = {radmin, themin, phimin};
  // const double zmax_arr[3] = {radmax, themax, phimax};

  // Derived parameters
  rbufftmp[0] = this->rminbuff[1];
  rbufftmp[1] = this->rminbuff[0];

  // setbuffcoef(rbufftmp, this->coefminbuff);
  // setbuffcoef(this->rmaxbuff, this->coefmaxbuff);

  int nlen = std::abs(this->filter_m1 - this->filter_m0) + 1;
  this->filter_m1d.resize(nlen);
  for (int i = 0; i < nlen; ++i) {
    this->filter_m1d[i] = std::min(filter_m0, filter_m1) + i;
  }

  if (rank == 0) {
    std::cout << " field filter_m1d= ";
    for (auto val : this->filter_m1d) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
  }

  this->nradfem = this->nrad + 2;
  this->nthefem = this->nthe;
  this->nphifem = this->nphi;

  this->ntotfem = this->nradfem * this->nthefem * this->nphifem;
  this->ntotdof = this->nrad * this->nthe * this->nphi;

  this->phiwid = 2.0 * M_PI / this->Nphimult;
  this->phimax = this->phimin + this->phiwid;

  this->radwid = this->radmax - this->radmin;
  this->thewid = this->themax - this->themin;

  this->bcrad = bcrad;
  this->lenntor = (this->ntor_max - this->ntor_min) / this->ntor_stride + 1;

  this->ntor1d.resize(this->lenntor);
  for (int i = 0; i < this->lenntor; ++i) {
    this->ntor1d[i] = this->ntor_min + i * this->ntor_stride;
  }

  this->ntotfem2d1f = this->nradfem * this->nthefem * this->lenntor;
  this->ntot12fem = this->nradfem * this->nthefem;

  // Print field information
  if (rank == 0) {
    std::cout << "================ Field Info ================" << std::endl;
    std::cout << "nl_debug=" << this->nl_debug
              << ", irestart=" << this->irestart << std::endl;
    std::cout << "Nphimult=" << this->Nphimult << ", nrad=" << this->nrad
              << ", nthe=" << this->nthe << ", nphi=" << this->nphi
              << std::endl;
    std::cout << "icorr=" << this->icorr << ", irandom_gy=" << this->irandom_gy
              << ", isolver=" << this->isolver << ", imodel=" << this->imodel
              << ", ischeme_Ohm=" << this->ischeme_Ohm << std::endl;
    std::cout << "radmin=" << this->radmin << ", radmax=" << this->radmax
              << ", thewid=" << this->thewid << ", phiwid=" << this->phiwid
              << std::endl;
    std::cout << "nradfem=" << this->nradfem << ", nthefem=" << this->nthefem
              << ", nphifem=" << this->nphifem << ", ntotfem=" << this->ntotfem
              << ", ntotdof=" << this->ntotdof << std::endl;
    std::cout << "ifilter=" << this->ifilter
              << ", filter_m0=" << this->filter_m0
              << ", filter_m1=" << this->filter_m1
              << ", filter_nc=" << this->filter_nc << std::endl;
    std::cout << "ntor_min=" << this->ntor_min
              << ", ntor_max=" << this->ntor_max
              << ", lenntor=" << this->lenntor
              << ", ntor_stride=" << this->ntor_stride << std::endl;
  }
}

void FieldCls::readInput(const std::string &inputFile) {
  INIReader reader(inputFile);

  // Check if reading was successful
  if (reader.ParseError() < 0) {
    throw std::runtime_error("Error: Unable to open or parse input file '" +
                             inputFile + "'");
  }

  std::cout << "field readInput: " << std::endl;

  // 读取参数
  nl_debug = reader.GetBoolean("Field", "nl_debug", false);
  irestart = reader.GetBoolean("Field", "irestart", false);
  Nphimult = reader.GetInteger("Field", "Nphimult", 1);
  nrad = reader.GetInteger("Field", "nrad", 1);
  nthe = reader.GetInteger("Field", "nthe", 1);
  nphi = reader.GetInteger("Field", "nphi", 1);
  radmin = reader.GetReal("Field", "radmin", 0.0);
  radmax = reader.GetReal("Field", "radmax", 1.0);
  ndiag_field_h5 = reader.GetInteger("Field", "ndiag_field_h5", 0);
  ifilter = reader.GetInteger("Field", "ifilter", 0);
  filter_m0 = reader.GetInteger("Field", "filter_m0", 0);
  filter_m1 = reader.GetInteger("Field", "filter_m1", 0);
  filter_nc = reader.GetInteger("Field", "filter_nc", 0);
  icorr = reader.GetInteger("Field", "icorr", 1);
  irandom_gy = reader.GetInteger("Field", "irandom_gy", 1);
  isolver = reader.GetInteger("Field", "isolver", 1);
  imodel = reader.GetInteger("Field", "imodel", 0);
  ischeme_Ohm = reader.GetInteger("Field", "ischeme_Ohm", 1);
  ntor_min = reader.GetInteger("Field", "ntor_min", 1);
  ntor_max = reader.GetInteger("Field", "ntor_max", 1);
  ntor_stride = reader.GetInteger("Field", "ntor_stride", 1);
  bcrad = reader.GetInteger("Field", "bcrad", 0);
  ibuff = reader.GetInteger("Field", "ibuff", 0);
  rminbuff[0] = reader.GetReal("Field", "rminbuff_1", 0.0);
  rminbuff[1] = reader.GetReal("Field", "rminbuff_2", 0.0);
  rmaxbuff[0] = reader.GetReal("Field", "rmaxbuff_1", 1.0);
  rmaxbuff[1] = reader.GetReal("Field", "rmaxbuff_2", 1.0);
  niter_sv = reader.GetInteger("Field", "niter_sv", 0);
  imixvar = reader.GetInteger("Field", "imixvar", 0);
  imatscheme = reader.GetInteger("Field", "imatscheme", 1);

  std::cout << "Read parameters successfully from " << inputFile << std::endl;
}

double FieldCls::field_cls_get_fbuff(double rad,
                                     const std::array<double, 2> &rbuff,
                                     const std::array<double, 4> &coef,
                                     bool left) const {
  double var = 0.0;

  if (rad >= rbuff[0] && rad <= rbuff[1]) {
    var = coef[0] + coef[1] * rad + coef[2] * rad * rad +
          coef[3] * rad * rad * rad;
  } else if (rad > rbuff[1] && left) {
    var = 1.0;
  } else if (rad < rbuff[0] && left) {
    var = 0.0;
  } else if (rad > rbuff[1] && !left) {
    var = 0.0;
  } else if (rad < rbuff[0] && !left) {
    var = 1.0;
  }

  return var;
}

// 针对double的重载
void FieldCls::field_cls_g2p2d1f_general(
    const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp, const double ptrad1d,
    const double ptthe1d, const double ptphi1d, double &ptf1d,
    const std::array<int, 3> &idiff, int ngyro, double rho1) {
  std::vector<double> ptf1d_vec;

  field_cls_g2p2d1f_general(equ, f1d, ntor1d, amp, std::vector<double>{ptrad1d},
                            std::vector<double>{ptthe1d},
                            std::vector<double>{ptphi1d}, ptf1d_vec, idiff,
                            ngyro, std::vector<double>{rho1});

  ptf1d = ptf1d_vec[0];
}

void FieldCls::field_cls_g2p2d1f_general(
    const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp,
    const std::vector<double> &ptrad1d, const std::vector<double> &ptthe1d,
    const std::vector<double> &ptphi1d, std::vector<double> &ptf1d,
    const std::array<int, 3> &idiff, int ngyro,
    const std::vector<double> &rho1) {
  // timer.tic(3);

  size_t np = ptrad1d.size();
  ptf1d.assign(np, 0.0); // Initialize ptf1d to zeros
  std::vector<std::complex<double>> f1d_amp;
  f1d_amp.resize(ntotfem2d1f, 0.0);
  for (int itor = 0; itor < ntor1d.size(); ++itor) {
    for (int ifem = 0; ifem < ntot12fem; ++ifem) {
      int idx = ifem + itor * ntot12fem;
      f1d_amp[idx] = amp[itor] * f1d[idx];
    }
  }

  if (ngyro <= 1) {
    spc.spc_cls_sp2p2d1f(f1d_amp, ntor1d, ptrad1d, ptthe1d, ptphi1d, ptf1d,
                         idiff);
  } else {
    std::vector<double> ptf1dgy(np, 0.0), ptrad1dgy(np), ptthe1dgy(np);
    std::vector<double> ptgyangle(np), ptR(np), ptZ(np);

    // Initialize gyro angles
    if (irandom_gy == 0) {
      std::fill(ptgyangle.begin(), ptgyangle.end(), 0.0);
    } else if (irandom_gy == 1) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist(0.0, 2 * M_PI);
      for (auto &angle : ptgyangle)
        angle = dist(gen);
    } else if (irandom_gy == 2) {
      ptgyangle = ptthe1d;
    } else if (irandom_gy == 3) {
      for (size_t i = 0; i < np; ++i) {
        ptgyangle[i] = ptthe1d[i] + M_PI_2;
      }
    }

    // Compute initial R and Z
    for (size_t i = 0; i < np; ++i) {
      ptR[i] = equ.getR(ptrad1d[i], ptthe1d[i]);
      ptZ[i] = equ.getZ(ptrad1d[i], ptthe1d[i]);
    }

    // Gyro averaging
    for (int fgc = 0; fgc < ngyro; ++fgc) {
      for (size_t i = 0; i < np; ++i) {
        double Rtmp = ptR[i] + cos(ptgyangle[i]) * rho1[i];
        double Ztmp = ptZ[i] + sin(ptgyangle[i]) * rho1[i];

        ptgyangle[i] += (2 * M_PI) / ngyro;

        ptrad1dgy[i] = equ.getradRZ(Rtmp, Ztmp);
        ptthe1dgy[i] = equ.gettheRZ(Rtmp, Ztmp);
      }

      std::vector<double> temp_ptf1dgy(np, 0.0);
      spc.spc_cls_sp2p2d1f(f1d_amp, ntor1d, ptrad1dgy, ptthe1dgy, ptphi1d,
                           temp_ptf1dgy, idiff);

      for (size_t i = 0; i < np; ++i) {
        ptf1d[i] += temp_ptf1dgy[i];
      }
    }

    // Normalize by ngyro
    for (auto &val : ptf1d) {
      val /= ngyro;
    }
  }

  // timer.toc(3);
}
// 针对double的重载
void FieldCls::field_cls_g2p2d1f_grad(
    const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amplitude_arr,
    const double ptrad1d, const double ptthe1d, const double ptphi1d,
    double &ptf1d100, double &ptf1d010, double &ptf1d001, int ngyro,
    double rho1) {
  std::vector<double> ptf1d100_vec, ptf1d010_vec, ptf1d001_vec;

  field_cls_g2p2d1f_grad(
      equ, f1d, ntor1d, amplitude_arr, std::vector<double>{ptrad1d},
      std::vector<double>{ptthe1d}, std::vector<double>{ptphi1d}, ptf1d100_vec,
      ptf1d010_vec, ptf1d001_vec, ngyro, std::vector<double>{rho1});

  ptf1d100 = ptf1d100_vec[0];
  ptf1d010 = ptf1d010_vec[0];
  ptf1d001 = ptf1d001_vec[0];
}

void FieldCls::field_cls_g2p2d1f_grad(
    const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amplitude_arr,
    const std::vector<double> &ptrad1d, const std::vector<double> &ptthe1d,
    const std::vector<double> &ptphi1d, std::vector<double> &ptf1d100,
    std::vector<double> &ptf1d010, std::vector<double> &ptf1d001, int ngyro,
    const std::vector<double> &rho1) {
  // timer.tic(3);

  size_t np = ptrad1d.size();
  ptf1d100.assign(np, 0.0);
  ptf1d010.assign(np, 0.0);
  ptf1d001.assign(np, 0.0);

  if (ngyro <= 1) {
    spc.spc_cls_sp2p2d1f_grad(f1d, ntor1d, amplitude_arr, ptrad1d, ptthe1d,
                              ptphi1d, ptf1d100, ptf1d010, ptf1d001);
  } else {
    std::vector<double> ptrad1dgy(np), ptthe1dgy(np);
    std::vector<double> ptgyangle(np), ptR(np), ptZ(np);

    // Initialize gyro angles
    if (irandom_gy == 0) {
      std::fill(ptgyangle.begin(), ptgyangle.end(), 0.0);
    } else if (irandom_gy == 1) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist(0.0, 2 * M_PI);
      for (auto &angle : ptgyangle)
        angle = dist(gen);
    } else if (irandom_gy == 2) {
      ptgyangle = ptthe1d;
    } else if (irandom_gy == 3) {
      for (size_t i = 0; i < np; ++i) {
        ptgyangle[i] = ptthe1d[i] + M_PI_2;
      }
    }

    // Compute initial R and Z
    for (size_t i = 0; i < np; ++i) {
      ptR[i] = equ.getR(ptrad1d[i], ptthe1d[i]);
      ptZ[i] = equ.getZ(ptrad1d[i], ptthe1d[i]);
    }

    // Gyro averaging
    for (int fgc = 0; fgc < ngyro; ++fgc) {
      for (size_t i = 0; i < np; ++i) {
        double Rtmp = ptR[i] + cos(ptgyangle[i]) * rho1[i];
        double Ztmp = ptZ[i] + sin(ptgyangle[i]) * rho1[i];

        ptgyangle[i] += (2 * M_PI) / ngyro;

        ptrad1dgy[i] = equ.getradRZ(Rtmp, Ztmp);
        ptthe1dgy[i] = equ.gettheRZ(Rtmp, Ztmp);
      }

      std::vector<double> temp_ptf1dgy100(np, 0.0), temp_ptf1dgy010(np, 0.0),
          temp_ptf1dgy001(np, 0.0);
      spc.spc_cls_sp2p2d1f_grad(f1d, ntor1d, amplitude_arr, ptrad1dgy,
                                ptthe1dgy, ptphi1d, temp_ptf1dgy100,
                                temp_ptf1dgy010, temp_ptf1dgy001);

      for (size_t i = 0; i < np; ++i) {
        ptf1d100[i] += temp_ptf1dgy100[i];
        ptf1d010[i] += temp_ptf1dgy010[i];
        ptf1d001[i] += temp_ptf1dgy001[i];
      }
    }

    // Normalize by ngyro
    for (size_t i = 0; i < np; ++i) {
      ptf1d100[i] /= ngyro;
      ptf1d010[i] /= ngyro;
      ptf1d001[i] /= ngyro;
    }
  }

  // timer.toc(3);
}
// calc grad and return complex value
// 针对double的重载
void FieldCls::field_cls_g2p2d1f_grad_complex(
    const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp, const double ptrad1d,
    const double ptthe1d, const double ptphi1d,
    std::vector<std::complex<double>> &ptf1d100_c,
    std::vector<std::complex<double>> &ptf1d010_c,
    std::vector<std::complex<double>> &ptf1d001_c, int ngyro, double rho1) {
  std::vector<std::complex<double>> ptf1d100_vec, ptf1d010_vec, ptf1d001_vec;

  field_cls_g2p2d1f_grad_complex(
      equ, f1d, ntor1d, amp, std::vector<double>{ptrad1d},
      std::vector<double>{ptthe1d}, std::vector<double>{ptphi1d}, ptf1d100_vec,
      ptf1d010_vec, ptf1d001_vec, ngyro, std::vector<double>{rho1});

  ptf1d100_c = ptf1d100_vec;
  ptf1d010_c = ptf1d010_vec;
  ptf1d001_c = ptf1d001_vec;
}

void FieldCls::field_cls_g2p2d1f_grad_complex(
    const Equilibrium &equ, const std::vector<std::complex<double>> &f1d,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp,
    const std::vector<double> &ptrad1d, const std::vector<double> &ptthe1d,
    const std::vector<double> &ptphi1d,
    std::vector<std::complex<double>> &ptf1d100_c,
    std::vector<std::complex<double>> &ptf1d010_c,
    std::vector<std::complex<double>> &ptf1d001_c, int ngyro,
    const std::vector<double> &rho1) {
  // timer.tic(3);

  size_t np = ptrad1d.size();
  std::complex<double> zero_c(0.0, 0.0);
  ptf1d100_c.assign(np, zero_c);
  ptf1d010_c.assign(np, zero_c);
  ptf1d001_c.assign(np, zero_c);

  if (ngyro <= 1) {
    spc.spc_cls_sp2p2d1f_grad_complex(f1d, ntor1d, amp, ptrad1d, ptthe1d,
                                      ptphi1d, ptf1d100_c, ptf1d010_c,
                                      ptf1d001_c);
  } else {
    std::vector<double> ptrad1dgy(np), ptthe1dgy(np);
    std::vector<double> ptgyangle(np), ptR(np), ptZ(np);

    // Initialize gyro angles
    if (irandom_gy == 0) {
      std::fill(ptgyangle.begin(), ptgyangle.end(), 0.0);
    } else if (irandom_gy == 1) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist(0.0, 2 * M_PI);
      for (auto &angle : ptgyangle)
        angle = dist(gen);
    } else if (irandom_gy == 2) {
      ptgyangle = ptthe1d;
    } else if (irandom_gy == 3) {
      for (size_t i = 0; i < np; ++i) {
        ptgyangle[i] = ptthe1d[i] + M_PI_2;
      }
    }

    // Compute initial R and Z
    for (size_t i = 0; i < np; ++i) {
      ptR[i] = equ.getR(ptrad1d[i], ptthe1d[i]);
      ptZ[i] = equ.getZ(ptrad1d[i], ptthe1d[i]);
    }

    // Gyro averaging
    for (int fgc = 0; fgc < ngyro; ++fgc) {
      for (size_t i = 0; i < np; ++i) {
        double Rtmp = ptR[i] + cos(ptgyangle[i]) * rho1[i];
        double Ztmp = ptZ[i] + sin(ptgyangle[i]) * rho1[i];

        ptgyangle[i] += (2 * M_PI) / ngyro;

        ptrad1dgy[i] = equ.getradRZ(Rtmp, Ztmp);
        ptthe1dgy[i] = equ.gettheRZ(Rtmp, Ztmp);
      }

      std::vector<std::complex<double>> temp_ptf1dgy100(np, zero_c),
          temp_ptf1dgy010(np, zero_c), temp_ptf1dgy001(np, zero_c);
      spc.spc_cls_sp2p2d1f_grad_complex(f1d, ntor1d, amp, ptrad1dgy, ptthe1dgy,
                                        ptphi1d, temp_ptf1dgy100,
                                        temp_ptf1dgy010, temp_ptf1dgy001);

      for (size_t i = 0; i < np; ++i) {
        ptf1d100_c[i] += temp_ptf1dgy100[i];
        ptf1d010_c[i] += temp_ptf1dgy010[i];
        ptf1d001_c[i] += temp_ptf1dgy001[i];
      }
    }

    // Normalize by ngyro
    for (size_t i = 0; i < np; ++i) {
      ptf1d100_c[i] /= ngyro;
      ptf1d010_c[i] /= ngyro;
      ptf1d001_c[i] /= ngyro;
    }
  }

  // timer.toc(3);
}
