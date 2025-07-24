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

#include "Equilibrium.h"

void Equilibrium::writeWallRZ() const {

  constexpr std::size_t nbbbs = 500;
  std::array<double, nbbbs> theta;
  std::array<double, nbbbs> rbbbs, zbbbs;
  double rad = 1.0;
  double dtheta = 2.0 * M_PI / static_cast<double>(nbbbs - 1);

  for (std::size_t i = 0; i < nbbbs; ++i) {
    theta[i] = -M_PI + i * dtheta;
    rbbbs[i] = getR(rad, theta[i]);
    zbbbs[i] = getZ(rad, theta[i]);
  }

  std::ofstream file;

  // Write rzbbbs.txt
  file.open("rzbbbs.txt");
  if (!file) {
    std::cerr << "Error: Failed to open rzbbbs.txt\n";
  } else {
    for (size_t fic = 0; fic < rbbbs.size(); ++fic) {
      file << std::scientific << std::setw(16) << std::setprecision(4)
           << rbbbs[fic] << std::setw(16) << std::setprecision(4) << zbbbs[fic]
           << "\n";
    }
    file.close();
  }

  // // Write rzlim.txt
  // file.open("rzlim.txt");
  // if (!file) {
  //     std::cerr << "Error: Failed to open rzlim.txt\n";
  // } else {
  //     for (size_t fic = 0; fic < rlim.size(); ++fic) {
  //         file << std::scientific << std::setw(16) << std::setprecision(4) <<
  //         rlim[fic]
  //              << std::setw(16) << std::setprecision(4) << zlim[fic] << "\n";
  //     }
  //     file.close();
  // }
}

void Equilibrium::equil_cls_calc_derived(int rank) {
  double vA, vth;

  if (rank == 0) {
    std::cout << ">>>>CALC_DERIVED, ISOURCE_REF= " << isource_ref << std::endl;
  }

  if (isource_ref == 1) {
    // giving rhoN betaN to derive nref, Tref
    Tref = (rhoN * charge_unit * Bref) * (rhoN * charge_unit * Bref) /
           (2.0 * mass_unit * charge_unit);
    nref = betaN * Bref * Bref / (2.0 * Tref * charge_unit * mu0);

    // Extra Info
    vA = Bref / std::sqrt(mu0 * nref * mass_unit);
    vth = std::sqrt(2.0 * Tref * charge_unit / mass_unit);
  } else if (isource_ref == 2) {
    // giving nref, Tref to derive rhoN betaN
    vA = Bref / std::sqrt(mu0 * nref * mass_unit);
    vth = std::sqrt(2.0 * Tref * charge_unit / mass_unit);

    rhoN =
        std::sqrt(2.0 * mass_unit * Tref * charge_unit) / (charge_unit * Bref);
    betaN = (vth * vth) / (vA * vA);
  }
  this->vN = vth; // save vN for use

  if (rank == 0) {
    std::cout << ">> vA=" << vA << "[m/s], v_N=" << vth
              << "[m/s], rhoN=" << rhoN << ", betaN=" << betaN
              << ", nref=" << nref << "[/m^3], Tref=" << Tref << "[eV]"
              << std::endl;
  }
}

void Equilibrium::equil_cls_showinfo(int rank) {
  if (rank == 0) {
    // Formatting strings (similar to Fortran sform_r and sform_i)
    std::cout << std::setfill(' ');

    // Integer format (e.g., "iequmodel=" followed by the value)
    std::cout << std::setw(14) << "iequmodel=" << std::setw(10) << iequmodel
              << '\n';

    // Real numbers with specific formatting (scientific notation)
    std::cout << std::setw(14) << "c1adhoc=" << std::scientific
              << std::setprecision(5) << std::setw(12) << c1adhoc
              << ", c2adhoc=" << c2adhoc << ", rmaxis_adhoc=" << rmaxis_adhoc
              << ", Bmaxis_adhoc=" << Bmaxis_adhoc << '\n';

    // Integer format for 'isource_ref'
    std::cout << std::setw(14) << "isource_ref=" << std::setw(10) << isource_ref
              << '\n';

    // Real numbers with specific formatting for the other parameters
    std::cout << std::setw(14) << "rhoN=" << rhoN << ", betaN=" << betaN
              << ", nref=" << nref << ", Tref=" << Tref << ", Bref=" << Bref
              << '\n';

    // Integer format for 'isetpsi'
    std::cout << std::setw(14) << "isetpsi=" << std::setw(10) << isetpsi
              << '\n';

    // Real number for 'fBphi'
    std::cout << std::setw(14) << "fBphi=" << fBphi << '\n';
  }
}
double Equilibrium::getradRZ(double RR, double ZZ) const {
  double var = 0.0;
  int idx, idy, iflag = 0;

  if (iequmodel == 1) {
    idx = 0;
    idy = 0;
    // radrz_bsp->evaluate(RR, ZZ, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getradRZ: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    var = std::sqrt((RR - rmaxis) * (RR - rmaxis) +
                    (ZZ - zmaxis) * (ZZ - zmaxis));
  }

  return var;
}

double Equilibrium::gettheRZ(double RR, double ZZ) const {
  double var = 0.0;
  double rormaj = 0.0;
  const double twopi = 2 * M_PI; // Define twopi as 2 * pi

  if (iequmodel == 1) {
    // For model 1, the logic is simpler
    var = std::atan2(ZZ - zmaxis, RR - rmaxis);
    var = UtilMath::modulo(var, twopi); // Use fmod to calculate the modulus
  } else if (iequmodel == 2) {
    // For model 2, more complex calculation
    rormaj = getradRZ(RR, ZZ) / rmaxis;
    var = std::atan2(ZZ - zmaxis, RR - rmaxis);
    var = 2 * std::atan(std::sqrt((1.0 - rormaj) / (1.0 + rormaj)) *
                        std::tan(var / 2.0));
    var = UtilMath::modulo(var, twopi);

    // Check for the value of var and reset to 0 if not in the valid range
    if (var < 0.0 || var > twopi) {
      var = 0.0;
      std::cerr << "Warning: var reset to 0 in gettheRZ (iequmodel == 2) from "
                << var << std::endl;
    }
  }

  return var;
}
void Equilibrium::calcdrtdRZ(double rad, double the, double &drdR, double &drdZ,
                             double &dtdR, double &dtdZ, double &jaco2) const {

  double dRdrad, dRdthe, dZdrad, dZdthe, jacoinv;

  // Get the partial derivatives of R and Z with respect to rad and theta
  dRdrad = getdRdrad(rad, the);
  dRdthe = getdRdthe(rad, the);
  dZdrad = getdZdrad(rad, the);
  dZdthe = getdZdthe(rad, the);

  // Calculate the Jacobian determinant (jaco2)
  jaco2 = dRdrad * dZdthe - dRdthe * dZdrad;

  // Calculate the inverse of the Jacobian
  jacoinv = 1.0 / jaco2;

  // Calculate the derivatives of r and theta with respect to R and Z
  drdR = dZdthe * jacoinv;
  drdZ = -dRdthe * jacoinv;
  dtdR = -dZdrad * jacoinv;
  dtdZ = dRdrad * jacoinv;
}

void Equilibrium::calcBJgij(double rad, double the0, double &RR, double &ZZ,
                            double &BB, double &dBdrad, double &dBdthe,
                            double &Bthe_ct, double &Bphi_ct, double &jaco2,
                            double &FF, double &g11, double &g12,
                            double &g22) const {
  // input rad, the0
  double the =
      UtilMath::modulo(the0, 2 * M_PI); // Modulo operation to wrap the angle
  double drdR, drdZ, dtdR, dtdZ;

  // Call to another member function
  calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2);

  g11 = drdR * drdR + drdZ * drdZ; // Calculate g11
  g22 = dtdR * dtdR + dtdZ * dtdZ; // Calculate g22
  g12 = drdR * dtdR + drdZ * dtdZ; // Calculate g12

  FF = getfpol_r(rad);
  RR = getR(rad, the);
  ZZ = getZ(rad, the);

  BB = getB(rad, the);
  dBdrad = getdBdrad(rad, the);
  dBdthe = getdBdthe(rad, the);

  if (iequmodel == 1) {
    Bthe_ct = getdpsidr(rad) / (jaco2 * RR);
    Bphi_ct = FF / (RR * RR);
  } else if (iequmodel == 2) {
    Bthe_ct = getBthe_ct(rad, the0);
    Bphi_ct = FF / (RR * RR);
    // Uncomment if you want to use getBphi_ct instead:
    // Bphi_ct = getBphi_ct(rad, the0);
  }
}

// Function to calculate h_adhoc
double Equilibrium::geth_adhoc(double rad, double the, double rmaxis) const {
  double rormaj = rad / rmaxis;
  double var = (1.0 - rormaj * rormaj) / (1.0 - rormaj * std::cos(the));
  return var;
}

// Function to calculate Bthe_ct
double Equilibrium::getBthe_ct(double rad, double the0) const {
  double var = 0.0;
  double drdR, drdZ, dtdR, dtdZ, jaco2, RR;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation

  if (iequmodel == 1) {
    calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2);
    RR = getR(rad, the);
    var = getdpsidr(rad) / (jaco2 * RR);
  } else if (iequmodel == 2) {
    var = Bmaxis /
          (rmaxis * getqloc_rt(rad, the0) * geth_adhoc(rad, the0, rmaxis) *
           geth_adhoc(rad, the0, rmaxis));
  }

  return var;
}

// Function to calculate Bphi_ct
double Equilibrium::getBphi_ct(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation

  if (iequmodel == 1) {
    double RR = getR(rad, the);
    var = getfpol_r(rad) / (RR * RR);
  } else if (iequmodel == 2) {
    double h = geth_adhoc(rad, the0, rmaxis);
    var = Bmaxis / (rmaxis * h * h);
  }

  return var;
}
// Function to calculate B
double Equilibrium::getB(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;
  int idx = 0, idy = 0, iflag = 0;

  if (iequmodel == 1) {
    // For model 1, using Brt_bsp to evaluate
    // Brt_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getB: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    double rormaj_factor = 1.0 - rormaj * rormaj;
    double nu_over_q = rormaj / getqloc_rt(rad, the0);
    var = Bmaxis * (1 - rormaj * std::cos(the0)) / rormaj_factor *
          std::sqrt(nu_over_q * nu_over_q / rormaj_factor + 1.0);
  }

  return var;
}
// Function to calculate dBdrad
double Equilibrium::getdBdrad(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0, qloc = 0.0;
  int idx = 1, idy = 0, iflag = 0;

  if (iequmodel == 1) {
    // For model 1, using Brt_bsp to evaluate
    // Brt_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getdBdrad: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    // For model 2, performing calculations based on formula
    rormaj = rad / rmaxis;
    qloc = getqloc_rt(rad, the0);
    double rormaj_factor = 1.0 - rormaj * rormaj;
    double nu_over_q = rormaj / qloc;
    double cos_the0 = std::cos(the0);

    var = ((1.0 / rad - getdqdrloc_rt(rad, the0) / qloc +
            3.0 * rormaj / rormaj_factor / rmaxis) *
               nu_over_q * nu_over_q +
           2.0 * rormaj / rmaxis) *
              (1.0 - rormaj * cos_the0) * Bmaxis /
              (rormaj_factor * rormaj_factor) -
          (nu_over_q * nu_over_q / rormaj_factor + 1.0) * cos_the0 * Bmaxis /
              (rormaj_factor * rmaxis);

    var = var / std::sqrt(nu_over_q * nu_over_q / rormaj_factor + 1.0);
  }

  return var;
}
// Function to calculate dB/dthe
double Equilibrium::getdBdthe(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0, qloc = 0.0;

  if (iequmodel == 1) {
    // For model 1, using Brt_bsp to evaluate
    int idx = 0, idy = 1, iflag = 0;
    // Brt_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getdBdrad: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    // For model 2, performing calculations based on formula
    rormaj = rad / rmaxis;
    qloc = getqloc_rt(rad, the0);
    double rormaj_factor = 1.0 - rormaj * rormaj;
    var = Bmaxis * rormaj * std::sin(the0) / rormaj_factor *
          std::sqrt(rormaj * rormaj / (qloc * qloc * rormaj_factor) + 1.0);
  }

  return var;
}

// Function to calculate Brad_co
double Equilibrium::getBrad_co(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;

  if (iequmodel == 1) {
    // For model 1, using Brad_co_bsp to evaluate
    int idx = 0, idy = 0, iflag = 0;
    // Brad_co_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getBrad_co: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    double rormaj_factor = 1.0 - rormaj * rormaj;
    var = rormaj * rormaj * Bmaxis * std::sin(the0) /
          (getqloc_rt(rad, the0) * rormaj_factor * rormaj_factor);
  }

  return var;
}

// Function to calculate dBrad_co/drad
double Equilibrium::getdBrad_codrad(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;

  if (iequmodel == 1) {
    // For model 1, using Brad_co_bsp to evaluate
    int idx = 1, idy = 0, iflag = 0;
    // Brad_co_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getdBrad_codrad: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    double rormaj_factor = 1.0 - rormaj * rormaj;
    var = (-rad * getdqdrloc_rt(rad, the0) / getqloc_rt(rad, the0) + 2.0 +
           4.0 * rormaj * rormaj / rormaj_factor) *
          rormaj * Bmaxis * std::sin(the0) /
          (rmaxis * getqloc_rt(rad, the0) * rormaj_factor * rormaj_factor);
  }

  return var;
}
// Function to calculate dBrad_co/dthe
double Equilibrium::getdBrad_codthe(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;

  if (iequmodel == 1) {
    // For model 1, using Brad_co_bsp to evaluate
    int idx = 0, idy = 1, iflag = 0;
    // Brad_co_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getdBrad_codthe: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    double rormaj_factor = 1.0 - rormaj * rormaj;
    var = std::cos(the0) * rormaj * rormaj * Bmaxis /
          (getqloc_rt(rad, the0) * rormaj_factor * rormaj_factor);
  }

  return var;
}

double Equilibrium::getBthe_co(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;

  if (iequmodel == 1) {
    // For model 1, using Bthe_co_bsp to evaluate
    int idx = 0, idy = 0, iflag = 0;
    // Bthe_co_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getBthe_co: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    var = rad * rad * Bmaxis /
          (rmaxis * getqloc_rt(rad, the0) * std::pow(1.0 - rormaj * rormaj, 2));
  }

  return var;
}

// Function to calculate Bphi_co
double Equilibrium::getBphi_co(double rad, double the0) const {
  double var = 0.0;

  if (iequmodel == 1) {
    var = getfpol_r(rad);
  } else if (iequmodel == 2) {
    var = Bmaxis * rmaxis;
  }

  return var;
}

// Function to calculate dBthe_codrad
double Equilibrium::getdBthe_codrad(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;

  if (iequmodel == 1) {
    int idx = 1, idy = 0, iflag = 0;
    // Bthe_co_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getdBthe_co/drad: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    double rormaj_factor = 1.0 - rormaj * rormaj;
    var = (-rad * getdqdrloc_rt(rad, the0) / getqloc_rt(rad, the0) + 2.0 +
           2.0 * rormaj * rormaj / rormaj_factor) *
          rormaj * Bmaxis / (getqloc_rt(rad, the0) * rormaj_factor);
  }

  return var;
}

// Function to calculate dBthe_codthe
double Equilibrium::getdBthe_codthe(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = 0.0;

  if (iequmodel == 1) {
    int idx = 0, idy = 1, iflag = 0;
    // Bthe_co_bsp.evaluate(rad, the, idx, idy, var, iflag);

    // Error handling
    if (iflag != 0) {
      std::cerr << "****Error in getdBthe_co/dthe: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    rormaj = rad / rmaxis;
    var = 0.0;
  }

  return var;
}
// Function to calculate Bdirect
double Equilibrium::getBdirect(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double drdR = 0.0, drdZ = 0.0, dtdR = 0.0, dtdZ = 0.0, jaco2 = 0.0;
  double FF = 0.0, g11 = 0.0, RR = 0.0;

  if (iequmodel == 1) {
    // Calculate derivatives and jacobian
    calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2);

    // Compute g11
    g11 = drdR * drdR + drdZ * drdZ;

    // Get poloidal field and radius
    FF = getfpol_r(rad);
    RR = getR(rad, the);

    // Compute the result based on the formula
    double dpsidr = getdpsidr(rad);
    var = std::sqrt(dpsidr * dpsidr * g11 + FF * FF) / RR;
  } else if (iequmodel == 2) {
    var = getB(rad, the0);
  }

  return var;
}

// Function to calculate Brad_co_direct
double Equilibrium::getBrad_co_direct(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double drdR = 0.0, drdZ = 0.0, dtdR = 0.0, dtdZ = 0.0, jaco2 = 0.0;
  double g12 = 0.0, RR = 0.0;

  if (iequmodel == 1) {
    calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2);

    g12 = drdR * dtdR + drdZ * dtdZ;
    RR = getR(rad, the);

    var = -jaco2 * g12 * getdpsidr(rad) / RR;
  } else if (iequmodel == 2) {
    var = getBrad_co(rad, the0);
  }

  return var;
}

double Equilibrium::getBthe_co_direct(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double drdR = 0.0, drdZ = 0.0, dtdR = 0.0, dtdZ = 0.0, jaco2 = 0.0;
  double g11 = 0.0, RR = 0.0;

  if (iequmodel == 1) {
    calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2);
    g11 = drdR * drdR + drdZ * drdZ;
    RR = getR(rad, the);

    var = jaco2 * g11 * getdpsidr(rad) / RR;
  } else if (iequmodel == 2) {
    var = getBthe_co(rad, the0);
  }

  return var;
}

double Equilibrium::getcos_straight_adhoc(double rad, double the0,
                                          double rmaxis) const {
  double rormaj = rad / rmaxis;
  double var = -1.0 / rormaj + (1.0 - rormaj * rormaj) /
                                   (rormaj * (1.0 - rormaj * std::cos(the0)));
  return var;
}

double Equilibrium::getsin_straight_adhoc(double rad, double the0,
                                          double rmaxis) const {
  double rormaj = rad / rmaxis;
  double var = std::sqrt(1.0 - rormaj * rormaj) * std::sin(the0) /
               (1.0 - rormaj * std::cos(the0));
  return var;
}

// Function to get R
double Equilibrium::getR(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  int idx = 0, idy = 0, iflag = 0;

  if (iequmodel == 1) {
    // rrt_bsp->evaluate(rad, the, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getR: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    var = rmaxis + rad * getcos_straight_adhoc(rad, the0, rmaxis);
  }

  return var;
}

// Function to get the derivative of R with respect to rad
double Equilibrium::getdRdrad(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  double rormaj = rad / rmaxis;

  if (iequmodel == 1) {
    // Use rrt_bsp to evaluate the derivative with respect to rad
    int idx = 1, idy = 0, iflag = 0;
    // rrt_bsp->evaluate(rad, the, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getdRdrad: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {

    double cos_the0 = std::cos(the0);
    double sin_the0 = std::sin(the0);
    double denominator = 1.0 - rormaj * cos_the0;
    var = getcos_straight_adhoc(rad, the0, rmaxis) -
          rormaj * sin_the0 * sin_the0 / (denominator * denominator);
  }

  return var;
}

// Function to get the derivative of R with respect to the angle (theta)
double Equilibrium::getdRdthe(double rad, double the0) const {
  double var = 0.0;
  double rormaj = rad / rmaxis;

  if (iequmodel == 1) {
    // Use rrt_bsp to evaluate the derivative with respect to theta
    int idx = 0, idy = 1, iflag = 0;
    // rrt_bsp->evaluate(rad, the0, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getdRdthe: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    double denominator = 1.0 - rormaj * std::cos(the0);
    var = -rad * (1.0 - rormaj * rormaj) * std::sin(the0) /
          (denominator * denominator);
  }

  return var;
}

// Function to get Z
double Equilibrium::getZ(double rad, double the0) const {
  double var = 0.0;
  double the = UtilMath::modulo(the0, 2 * M_PI); // Modulo operation
  int idx = 0, idy = 0, iflag = 0;

  if (iequmodel == 1) {
    // Use zrt_bsp to evaluate Z
    // zrt_bsp->evaluate(rad, the0, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getZ: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    var = zmaxis + rad * getsin_straight_adhoc(rad, the0, rmaxis);
  }

  return var;
}

// Function to get the derivative of Z with respect to rad
double Equilibrium::getdZdrad(double rad, double the0) const {
  double var = 0.0;
  double rormaj = rad / rmaxis;

  if (iequmodel == 1) {
    // Use zrt_bsp to evaluate the derivative with respect to rad
    int idx = 1, idy = 0, iflag = 0;
    // zrt_bsp->evaluate(rad, the0, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getdZdrad: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    double cos_the0 = std::cos(the0);
    double sin_the0 = std::sin(the0);
    double denominator = 1.0 - rormaj * cos_the0;
    denominator *= denominator;
    denominator *= sqrt(1.0 - rormaj * rormaj);

    var = getsin_straight_adhoc(rad, the0, rmaxis) +
          rormaj * sin_the0 * (cos_the0 - rormaj) / denominator;
  }

  return var;
}

// Function to get the derivative of Z with respect to theta
double Equilibrium::getdZdthe(double rad, double the0) const {
  double var = 0.0;
  double rormaj = rad / rmaxis;

  if (iequmodel == 1) {
    // Use zrt_bsp to evaluate the derivative with respect to theta
    int idx = 0, idy = 1, iflag = 0;
    // zrt_bsp->evaluate(rad, the0, idx, idy, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getdZdthe: " << iflag << std::endl;
    }
  } else if (iequmodel == 2) {
    double cos_the0 = std::cos(the0);
    double denominator = 1.0 - rormaj * cos_the0;
    denominator *= denominator; // Square the denominator

    var = rad * std::sqrt(1.0 - rormaj * rormaj) / denominator *
          (cos_the0 - rormaj);
  }

  return var;
}

// Function to compute dqdrloc_rt
double Equilibrium::getdqdrloc_rt(double rad, double the) const {
  double var = 0.0;
  double rormaj = rad / rmaxis;

  if (iequmodel == 1) {
    std::cerr << "error: not implemented" << std::endl;
    exit(1); // Exit if model 1 is not implemented
  } else if (iequmodel == 2) {
    var = 2.0 * c2adhoc * rad + (c1adhoc + c2adhoc * rad * rad) * rormaj /
                                    ((1.0 - rormaj * rormaj) * rmaxis);
    var /= std::sqrt(1.0 - rormaj * rormaj);
  }

  return var;
}

// Function to compute qloc_rt
double Equilibrium::getqloc_rt(double rad, double the) const {
  double var = 0.0;
  double drdR, drdZ, dtdR, dtdZ, jaco2, RR;

  if (iequmodel == 1) {
    // Call calcdrtdRZ for model 1
    calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2);
    RR = getR(rad, the);
    var = jaco2 * getfpol_r(rad) / (RR * getdpsidr(rad));
  } else if (iequmodel == 2) {
    double rormaj = rad / rmaxis;
    var = c1adhoc + c2adhoc * rad * rad;
    var /= std::sqrt(1.0 - rormaj * rormaj);
  }

  return var;
}

// Function to compute psi
double Equilibrium::getpsi(double rad) const {
  double var = 0.0;

  if (iequmodel == 1) {
    var = simag + psiwid * rad * rad;
  } else if (iequmodel == 2) {
    var = Bmaxis / (2.0 * c2adhoc) *
          std::log(1.0 + (c2adhoc / c1adhoc) * rad * rad);
  }

  return var;
}

// Function to compute dpsidr
double Equilibrium::getdpsidr(double rad) const {
  double var = 0.0;

  if (iequmodel == 1) {
    var = 2.0 * psiwid * rad;
  } else if (iequmodel == 2) {
    var = Bmaxis * rad / (c1adhoc + c2adhoc * rad * rad);
  }

  return var;
}
// Function to compute fpol_r
double Equilibrium::getfpol_r(double rad) const {
  double var = 0.0;
  int idx = 0;
  int iflag = 0;

  if (iequmodel == 1) {
    // F1d_bsp.evaluate(getpsi(rad), idx, var, iflag);
    if (iflag != 0) {
      std::cerr << "****Error in getfpol_r: " << iflag << ", r=" << rad
                << ", psi=" << getpsi(rad) << std::endl;
    }
  } else if (iequmodel == 2) {
    var = Bmaxis * rmaxis;
  }

  return var;
}

// Function to compute dfpoldr_r
double Equilibrium::getdfpoldr_r(double rad) const {
  double var = 0.0;
  int idx = 1;
  int iflag = 0;

  if (iequmodel == 1) {
    // Call evaluate function of F1d_bsp for model 1
    // F1d_bsp.evaluate(getpsi(rad), idx, var, iflag);
    var = var * 2.0 * psiwid * rad;
    if (iflag != 0) {
      std::cerr << "****Error in getdfpoldr_r: " << iflag << ", r=" << rad
                << ", psi=" << getpsi(rad) << std::endl;
    }
  } else if (iequmodel == 2) {
    var = 0.0;
  }

  return var;
}

// Function to compute jaco2
double Equilibrium::getjaco2(double rad, double the) const {
  double var = 0.0;
  double dRdrad = getdRdrad(rad, the);
  double dRdthe = getdRdthe(rad, the);
  double dZdrad = getdZdrad(rad, the);
  double dZdthe = getdZdthe(rad, the);

  var = dRdrad * dZdthe - dRdthe * dZdrad;

  return var;
}

// Function to compute jaco3
double Equilibrium::getjaco3(double rad, double the) const {
  double var = 0.0;
  double dRdrad = getdRdrad(rad, the);
  double dRdthe = getdRdthe(rad, the);
  double dZdrad = getdZdrad(rad, the);
  double dZdthe = getdZdthe(rad, the);
  double RR = getR(rad, the);

  var = RR * (dRdrad * dZdthe - dRdthe * dZdrad);

  return var;
}

double Equilibrium::getcurlbrad_ct(double rad, double the) const {
  double var = 0.0;
  double BB = getB(rad, the);

  var = getBphi_co(rad, the) * getdBdthe(rad, the) /
        (getjaco3(rad, the) * BB * BB);
  return var;
}
double Equilibrium::getcurlbthe_ct(double rad, double the0) const {
  double var = 0.0;
  double BB = getB(rad, the0);
  var = -getBphi_co(rad, the0) * getdBdrad(rad, the0) /
        (getjaco3(rad, the0) * BB * BB);
  return var;
}
double Equilibrium::getcurlbphi_ct(double rad, double the0) const {
  double BB = getB(rad, the0);
  double jaco3 = getjaco3(rad, the0);
  double var =
      -(getBrad_co(rad, the0) * getdBdthe(rad, the0) -
        getBthe_co(rad, the0) * getdBdrad(rad, the0)) /
          (jaco3 * BB * BB) +
      (getdBrad_codthe(rad, the0) - getdBthe_codrad(rad, the0)) / (jaco3 * BB);
  return var;
}