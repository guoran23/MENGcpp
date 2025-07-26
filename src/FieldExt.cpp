#include "FieldExt.h"
void FieldExtCls::field_ext_cls_calc_W(
    std::vector<double> &WWW, const Equilibrium &equ, Particle &pt,
    const std::vector<std::complex<double>> &phik_c,
    const std::vector<int> &ntor1d) {
  int nsp = pt.getNsp();
  int nrad = this->nrad;
  int nthe = this->nthe;
  double rmax = this->radmax;
  double rmin = this->radmin;
  double Bref = equ.getBref();
  lenntor = ntor1d.size();
  WWW.assign(lenntor, 0.0);
  double mass_ion = pt.mass_bk;
  std::vector<std::complex<double>> amp_1(lenntor,
                                          std::complex<double>(1.0, 0.0));

  std::vector<double> xrad(100), wrad(100);
  std::vector<double> xthe(200), wthe(200);
  UtilMath::lgwt(100, rmin, rmax, xrad, wrad);
  double twopi = 2.0 * M_PI;
  UtilMath::lgwt(200, 0.0, twopi, xthe, wthe);
  // Nested loop for computation
  double volume = 0.0;
  double int_dens = 0.0;
  for (int fic = 0; fic < 100; ++fic) {
    for (int fjc = 0; fjc < 200; ++fjc) {
      double w12 = wrad[fic] * wthe[fjc];
      double ptrad = xrad[fic];
      double ptthe = xthe[fjc];
      double ptphi = 0.0; //
      double ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe;
      double ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22, dfdpar;
      equ.calcBJgij(ptrad, ptthe, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                    ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22);

      double jaco3 = ptjaco2 * ptRR;

      std::vector<std::complex<double>> dfdrad_c, dfdthe_c, dfdphi_c;
      field_cls_g2p2d1f_grad_complex(equ, phik_c, ntor1d, amp_1, ptrad, ptthe,
                                     ptphi, dfdrad_c, dfdthe_c, dfdphi_c, 1,
                                     0.0);

      double dens_rad = pt.getIondens1d(ptrad);
      double ptB2 = ptB * ptB;
      for (int itor = 0; itor < lenntor; ++itor) {
        // Calculate W for each toroidal mode
        double nabla_perp_phi =
            ptg11 * std::norm(dfdrad_c[itor]) +
            2.0 * ptg12 *
                (dfdrad_c[itor].real() * dfdthe_c[itor].real() +
                 dfdrad_c[itor].imag() * dfdthe_c[itor].imag()) +
            ptg22 * std::norm(dfdthe_c[itor]);

        nabla_perp_phi *= mass_ion * dens_rad / (ptB2);
        WWW[itor] += nabla_perp_phi * w12 * jaco3;
      }

    } // nthe
  }   // nrad

  for (int itor = 0; itor < lenntor; ++itor) {
    // integrate over toroidal angle
    WWW[itor] *= M_PI  * Bref * Bref;
  }
}



// Implementation of init() method
void FieldExtCls::init(const Equilibrium &equ, Particle &pt) {

  // Initialize MPI rank and shared memory parameters
  MPIManager &mpiManager = MPIManager::getInstance();
  rank = mpiManager.getRank();
  size = mpiManager.getSize();

  // Initialize parent class
  FieldCls::field_cls_init(equ);

  std::cout << "----init FieldExtCls----" << std::endl;

  this->ntotfem2d1f = this->lenntor * this->spc.get_ntot12fem();
  this->ntotdof2d1f = this->lenntor * this->spc.get_ntot12dof();
  std::cout << "lenntor=" << this->lenntor << std::endl;
  std::cout << "ntot12fem=" << this->spc.get_ntot12fem() << std::endl;
  std::cout << "ntot12dof=" << this->spc.get_ntot12dof() << std::endl;
  std::cout << "ntotfem2d1f=" << this->ntotfem2d1f << std::endl;
  std::cout << "ntotdof2d1f=" << this->ntotdof2d1f << std::endl;

  // Allocate index array
  this->idxdof2d1f.resize(this->ntotdof2d1f);
  this->spc.spc_cls_calc_idxdof2d1f(this->lenntor, this->idxdof2d1f);

  // Allocate field and moment arrays
  this->denskMkj.resize(this->ntotfem2d1f);
  this->phik.resize(this->ntotfem2d1f);
  this->jparkMkj.resize(this->ntotfem2d1f);
  this->apark.resize(this->ntotfem2d1f);
  this->amplitude_arr.assign(
      lenntor, std::complex<double>(0.0, 0.0)); // Default amplitude

  // Initialize perturbation fields==0
  for (int i = 0; i < this->ntotfem2d1f; ++i) {
    this->denskMkj[i] = 0.0;
    this->phik[i] = 0.0;
    this->jparkMkj[i] = 0.0;
    this->apark[i] = 0.0;
  }
  // Initialize Ana Gauss perturbation fields
  initializePerturbations(equ);
}

void FieldExtCls::initializePerturbations(const Equilibrium &equ) {
  double dtheta = 2 * M_PI / nthe; // periodic in theta
  double drad = (radmax - radmin) / (nrad - 1);
  if (rank == 0) {
    std::cout << "----Initializing perturbations with nrad=" << nrad
              << ", nthe=" << nthe << std::endl;
    std::cout << "radmin=" << radmin << ", radmax=" << radmax
              << ", themin=" << themin << ", themax=" << themax << std::endl;
    std::cout << "drad=" << drad << ", dtheta=" << dtheta << std::endl;
    std::cout << "lenntor=" << lenntor << std::endl;
    for (size_t i = 0; i < ntor1d.size(); ++i) {
      std::cout << "ntor1d[" << i << "] = " << ntor1d[i] << std::endl;
    }
  }

  // read perturbation parameters from input file
  INIReaderExt reader("input.ini");

  if (reader.ParseError() != 0) {
    std::cerr << "Can't load 'input.ini'\n";
  }
  std::vector<double> omega_arr =
      reader.GetRealList("Perturbation", "omega", {0.1}); //
  std::vector<double> amp_arr = reader.GetRealList(
      "Perturbation", "amplitude", {0.0}); // eigenmode amplitude
  std::vector<std::vector<int>> mpoloidal_arr;
  std::vector<std::vector<double>> rc_arr, rwidth_arr, mpoloidal_amp_arr;

  for (size_t i = 1; i <= omega_arr.size(); ++i) {
    // mpoloidal_mode
    std::string key_mode = "mpoloidal_mode" + std::to_string(i);
    std::vector<int> mode_list =
        reader.GetIntList("Perturbation", key_mode, {});
    if (mode_list.empty()) {
      throw std::runtime_error("Missing " + key_mode);
    }
    mpoloidal_arr.push_back(mode_list);

    // rc
    std::string key_rc = "rc_" + std::to_string(i);
    std::vector<double> rc_list =
        reader.GetRealList("Perturbation", key_rc, {});
    if (rc_list.empty())
      throw std::runtime_error("Missing or invalid " + key_rc);
    rc_arr.push_back(rc_list);

    // rwidth
    std::string key_rwidth = "rwidth_" + std::to_string(i);
    std::vector<double> rwidth_list =
        reader.GetRealList("Perturbation", key_rwidth, {});
    if (rwidth_list.empty())
      throw std::runtime_error("Missing or invalid " + key_rwidth);
    rwidth_arr.push_back(rwidth_list);

    // mpoloidal_amp
    std::string key_amp = "mpoloidal_amp" + std::to_string(i);
    std::vector<double> amp_list =
        reader.GetRealList("Perturbation", key_amp, {-1.0});
    if (amp_list.empty())
      throw std::runtime_error("Missing or invalid " + key_amp);
    mpoloidal_amp_arr.push_back(amp_list);
  }

  // 检查 omega_arr 和 amp_arr 的大小是否与 lenntor 匹配
  if (omega_arr.size() != this->lenntor || amp_arr.size() != this->lenntor ||
      rc_arr.size() != this->lenntor || rwidth_arr.size() != this->lenntor) {
    std::cerr
        << "Error: omega_arr and amp_arr must have the same size as lenntor.\n";
    return;
  }
  //
  this->omega0 = omega_arr; // Store omega0 for later use
  for (size_t i = 0; i < amp_arr.size(); ++i) {
    this->amplitude_arr[i] = std::complex<double>(amp_arr[i], 0.0);
  }

  if (rank == 0) {
    std::cout << "Initializing perturbations with rc_arr, amp_arr, "
                 "rwidth_arr, mpoloidal_arr\n";
    for (int i = 0; i < this->lenntor; ++i) {
      std::cout << "itor=" << i << ", ntor=" << ntor1d[i]
                << ", amp=" << amp_arr[i] << ", omega=" << omega_arr[i]
                << ", mpoloidal=[";
      for (size_t j = 0; j < mpoloidal_arr[i].size(); ++j) {
        std::cout << mpoloidal_arr[i][j];
        if (j != mpoloidal_arr[i].size() - 1)
          std::cout << ", ";
      }
      std::cout << "]"
                << ", rc=[";
      for (size_t j = 0; j < rc_arr[i].size(); ++j) {
        std::cout << rc_arr[i][j];
        if (j != rc_arr[i].size() - 1)
          std::cout << ", ";
      }
      std::cout << "], rwidth=[";
      for (size_t j = 0; j < rwidth_arr[i].size(); ++j) {
        std::cout << rwidth_arr[i][j];
        if (j != rwidth_arr[i].size() - 1)
          std::cout << ", ";
      }
      std::cout << "], mpoloidal_amp=[";
      for (size_t j = 0; j < mpoloidal_amp_arr[i].size(); ++j) {
        std::cout << mpoloidal_amp_arr[i][j];
        if (j != mpoloidal_amp_arr[i].size() - 1)
          std::cout << ", ";
      }
      std::cout << "]"

                << std::endl;
    }
  }
  // give a MHD perturbation delta A_//= -i/omega* \partial_// delta Phi
  for (int itor = 0; itor < this->lenntor; ++itor) {
    double omega = omega_arr[itor];
    for (int impol = 0; impol < mpoloidal_arr[itor].size(); ++impol) {
      int mpoloidal = mpoloidal_arr[itor][impol];
      double rc = rc_arr[itor][impol];
      double rwidth = rwidth_arr[itor][impol];
      double amp_nm = mpoloidal_amp_arr[itor][impol];

      for (int j = 0; j < nthefem; ++j) {
        for (int i = 0; i < nradfem; ++i) {

          int idx = itor * nradfem * nthefem + j * nradfem + i;

          double theta = themin + j * dtheta;
          double r, radial_part;

          // 高斯径向部分
          if (i == 0) {
            radial_part = 0.0; // 边界条件
          } else if (i == nradfem - 1) {
            radial_part = 0.0; // 边界条件
          } else {
            r = radmin + (i - 1) * drad;
            radial_part = amp_nm * std::exp(-std::pow((r - rc) / rwidth, 2));
          }

          // 角向复数部分
          std::complex<double> angular_part =
              std::exp(std::complex<double>(0, mpoloidal * theta));

          this->phik[idx] += radial_part * angular_part;

          // 计算对应的 A_//=-i/omega * \partial_// delta Phi
          double ww = 0.0;
          if (i == 0 || i == nradfem - 1) {
            ww = 0.0; // 边界条件
          } else {
            double ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe;
            double ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22;
            equ.calcBJgij(r, theta, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                          ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12,
                          ptg22);
            ww = (mpoloidal * ptBthe_ct + ntor1d[itor] * ptBphi_ct) / ptB /
                 omega;
          }
          this->apark[idx] +=
              std::complex<double>(ww, 0) * radial_part * angular_part;
        } // nradfem
      }   // nthefem
    }     // mpoloidal
  }       // lenntor

  if (rank == 0) {
    std::vector<double> phik_real(phik.size());
    std::transform(phik.begin(), phik.end(), phik_real.begin(),
                   [](const std::complex<double> &z) { return std::real(z); });

    UtilIO::write_arr1d_r8(phik_real, "data_phik_initial.txt", true);

    // 写出 apark 的实部
    std::vector<double> apark_real(apark.size());
    std::transform(apark.begin(), apark.end(), apark_real.begin(),
                   [](const std::complex<double> &z) { return std::real(z); });
    UtilIO::write_arr1d_r8(apark_real, "data_apark_initial.txt", true);
  }
}
void FieldExtCls::field_ext_cls_test(const Equilibrium &equ,
                                     ParticleSpecies &pt) {
  // Constants
  int igradient = 1;
  // initial amplitude
  std::vector<std::complex<double>> amp(this->lenntor, 1.0);

  // Read Fortran-like input.field
  if (rank == 0)
    std::cout << "--------Field_Ext_EM2D1F Cls_Test--------\n";

  if (rank == 0) {
    std::ofstream out("data_field2particle.txt");
    out << std::scientific << std::setprecision(2);
    out << std::setw(10) << "ptrad" << std::setw(10) << "ptthe" << std::setw(10)
        << "ptphi" << std::setw(10) << "ptRR" << std::setw(10) << "ptZZ"
        << std::setw(10) << "ptf1d" << std::setw(10) << "dfdrad"
        << std::setw(10) << "dfdthe" << std::setw(10) << "dfdphi"
        << std::setw(10) << "dfdpar" << '\n';

    for (int i3 = 0; i3 < spc.getNnodeArr()[2]; ++i3) {
      for (int i2 = 0; i2 < spc.getNnodeArr()[1]; ++i2) {
        for (int i1 = 0; i1 < spc.getNnodeArr()[0]; ++i1) {
          std::vector<double> ptrad(1, spc.get_z1d1()[i1]);
          std::vector<double> ptthe(1, spc.get_z1d2()[i2]);
          std::vector<double> ptphi(1, spc.get_z1d3()[i3]);
          std::vector<double> prho(1, 0.0);

          std::vector<double> ff(1), dfdrad(1), dfdthe(1), dfdphi(1), dfdpar(1);
          std::vector<double> ptRR(1), ptZZ(1), ptB(1), ptdBdrad(1),
              ptdBdthe(1);
          std::vector<double> ptBthe_ct(1), ptBphi_ct(1), ptjaco2(1), ptFF(1),
              ptg11(1), ptg12(1), ptg22(1);

          // std::cout << "ptrad=" << ptrad[0] << ", ptthe=" << ptthe[0]
          // << ", ptphi=" << ptphi[0] << '\n';
          field_cls_g2p2d1f_general(equ, phik, ntor1d, amp, ptrad, ptthe, ptphi,
                                    ff, std::array<int, 3>{0, 0, 0}, 1, prho);

          if (igradient == 0) {
            field_cls_g2p2d1f_general(equ, phik, ntor1d, amp, ptrad, ptthe,
                                      ptphi, dfdrad,
                                      std::array<int, 3>{1, 0, 0}, 1, prho);
            field_cls_g2p2d1f_general(equ, phik, ntor1d, amp, ptrad, ptthe,
                                      ptphi, dfdthe,
                                      std::array<int, 3>{0, 1, 0}, 1, prho);
            field_cls_g2p2d1f_general(equ, phik, ntor1d, amp, ptrad, ptthe,
                                      ptphi, dfdphi,
                                      std::array<int, 3>{0, 0, 1}, 1, prho);
          } else if (igradient == 1) {
            field_cls_g2p2d1f_grad(equ, phik, ntor1d, amp, ptrad, ptthe, ptphi,
                                   dfdrad, dfdthe, dfdphi, 1, prho);
          }

          equ.calcBJgij(ptrad[0], ptthe[0], ptRR[0], ptZZ[0], ptB[0],
                        ptdBdrad[0], ptdBdthe[0], ptBthe_ct[0], ptBphi_ct[0],
                        ptjaco2[0], ptFF[0], ptg11[0], ptg12[0], ptg22[0]);

          dfdpar[0] =
              (ptBthe_ct[0] * dfdthe[0] + ptBphi_ct[0] * dfdphi[0]) / ptB[0];

          out << std::setw(10) << ptrad[0] << std::setw(10) << ptthe[0]
              << std::setw(10) << ptphi[0] << std::setw(10) << ptRR[0]
              << std::setw(10) << ptZZ[0] << std::setw(10) << ff[0]
              << std::setw(10) << dfdrad[0] << std::setw(10) << dfdthe[0]
              << std::setw(10) << dfdphi[0] << std::setw(10) << dfdpar[0]
              << '\n';
        }
      }
    }
  }
  // write apark and gradient
  if (rank == 0) {
    std::ofstream out("data_field2particle_apark.txt");
    out << std::scientific << std::setprecision(2);
    out << std::setw(10) << "ptrad" << std::setw(10) << "ptthe" << std::setw(10)
        << "ptphi" << std::setw(10) << "ptRR" << std::setw(10) << "ptZZ"
        << std::setw(10) << "ptf1d" << std::setw(10) << "dfdrad"
        << std::setw(10) << "dfdthe" << std::setw(10) << "dfdphi"
        << std::setw(10) << "dfdpar" << '\n';

    for (int i3 = 0; i3 < spc.getNnodeArr()[2]; ++i3) {
      for (int i2 = 0; i2 < spc.getNnodeArr()[1]; ++i2) {
        for (int i1 = 0; i1 < spc.getNnodeArr()[0]; ++i1) {
          std::vector<double> ptrad(1, spc.get_z1d1()[i1]);
          std::vector<double> ptthe(1, spc.get_z1d2()[i2]);
          std::vector<double> ptphi(1, spc.get_z1d3()[i3]);
          std::vector<double> prho(1, 0.0);

          std::vector<double> ff(1), dfdrad(1), dfdthe(1), dfdphi(1), dfdpar(1);
          std::vector<double> ptRR(1), ptZZ(1), ptB(1), ptdBdrad(1),
              ptdBdthe(1);
          std::vector<double> ptBthe_ct(1), ptBphi_ct(1), ptjaco2(1), ptFF(1),
              ptg11(1), ptg12(1), ptg22(1);

          // std::cout << "ptrad=" << ptrad[0] << ", ptthe=" << ptthe[0]
          // << ", ptphi=" << ptphi[0] << '\n';
          field_cls_g2p2d1f_general(equ, apark, ntor1d, amp, ptrad, ptthe,
                                    ptphi, ff, std::array<int, 3>{0, 0, 0}, 1,
                                    prho);

          if (igradient == 0) {
            field_cls_g2p2d1f_general(equ, apark, ntor1d, amp, ptrad, ptthe,
                                      ptphi, dfdrad,
                                      std::array<int, 3>{1, 0, 0}, 1, prho);
            field_cls_g2p2d1f_general(equ, apark, ntor1d, amp, ptrad, ptthe,
                                      ptphi, dfdthe,
                                      std::array<int, 3>{0, 1, 0}, 1, prho);
            field_cls_g2p2d1f_general(equ, apark, ntor1d, amp, ptrad, ptthe,
                                      ptphi, dfdphi,
                                      std::array<int, 3>{0, 0, 1}, 1, prho);
          } else if (igradient == 1) {
            field_cls_g2p2d1f_grad(equ, apark, ntor1d, amp, ptrad, ptthe, ptphi,
                                   dfdrad, dfdthe, dfdphi, 1, prho);
          }

          equ.calcBJgij(ptrad[0], ptthe[0], ptRR[0], ptZZ[0], ptB[0],
                        ptdBdrad[0], ptdBdthe[0], ptBthe_ct[0], ptBphi_ct[0],
                        ptjaco2[0], ptFF[0], ptg11[0], ptg12[0], ptg22[0]);

          dfdpar[0] =
              (ptBthe_ct[0] * dfdthe[0] + ptBphi_ct[0] * dfdphi[0]) / ptB[0];

          out << std::setw(10) << ptrad[0] << std::setw(10) << ptthe[0]
              << std::setw(10) << ptphi[0] << std::setw(10) << ptRR[0]
              << std::setw(10) << ptZZ[0] << std::setw(10) << ff[0]
              << std::setw(10) << dfdrad[0] << std::setw(10) << dfdthe[0]
              << std::setw(10) << dfdphi[0] << std::setw(10) << dfdpar[0]
              << '\n';
        }
      }
    }
  }

  if (rank == 0) {

    std::vector<double> ptrad = pt.getCoords().partrad;
    std::vector<double> ptthe = pt.getCoords().parttheta;
    std::vector<double> ptphi = pt.getCoords().partphitor;
    std::vector<double> prho1d(ptrad.size(), 0.0);

    std::vector<double> ptf1d; // Output vector to be filled

    field_cls_g2p2d1f_general(equ, phik, ntor1d, amp, ptrad, ptthe, ptphi,
                              ptf1d, std::array<int, 3>{0, 0, 0}, 1, prho1d);
    std::vector<double> dfdrad1d, dfdthe1d,
        dfdphi1d; // Gradients output vectors
    field_cls_g2p2d1f_grad(equ, phik, ntor1d, amp, ptrad, ptthe, ptphi,
                           dfdrad1d, dfdthe1d, dfdphi1d, 1, prho1d);

    // (1) open the file – overwrite each run;
    std::ofstream fout("data_phik2particle.txt",
                       std::ios::out | std::ios::trunc);
    if (!fout) {
      throw std::runtime_error("Cannot open data_phik2particle.txt");
    }
    fout << "# idx ptrad ptthe ptphi ptRR ptZZ ptf1d dfdrad dfdthe dfdphi "
            "dfdpar\n";

    // write one line per particle
    fout << std::setprecision(10);
    for (std::size_t i = 0; i < ptf1d.size(); ++i) {
      double ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe;
      double ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22, dfdpar;
      equ.calcBJgij(ptrad[i], ptthe[i], ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                    ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22);
      dfdpar = (ptBthe_ct * dfdthe1d[i] + ptBphi_ct * dfdphi1d[i]) / ptB;
      // Write to file
      fout << i << ' ' << ptrad[i] << ' ' << ptthe[i] << ' ' << ptphi[i] << ' '
           << ptRR << ' ' << ptZZ << ' ' << ptf1d[i] << ' ' << dfdrad1d[i]
           << ' ' << dfdthe1d[i] << ' ' << dfdphi1d[i] << dfdpar << '\n';
    }

    // (4) automatically closes when fout goes out of scope
  }

  if (rank == 0)
    std::cout << "--------Field_Ext_Cls_Test Done--------\n";
}
