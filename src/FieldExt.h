#ifndef FIELDEXTCLS_H
#define FIELDEXTCLS_H

#include "Equilibrium.h"
#include "Field.h" // Base class
#include "MPIManager.h"
#include "Particle.h"
#include <complex>
#include <vector>

class FieldExtCls : public FieldCls {
private:
  int rank, size;

public:
  // Fields & moments arrays
  int ntotfem2d1f, ntotdof2d1f;
  std::vector<int> idxdof2d1f;

  std::vector<std::complex<double>> denskMkj, phik, jparkMkj, apark;
  // std::vector<std::complex<double>> aparsk, aparhk;
  // int imixvar = 0; //0, no mixvar, 1 mixing variable

  // Constructor & Destructor
  FieldExtCls() {
    std::cout << "FieldExtCls default constructor called" << std::endl;
  };
  ~FieldExtCls() = default;

  // Methods
  void init(const Equilibrium &equ, std::vector<ParticleSpecies> &pt);
  void initializePerturbations(const Equilibrium &equ);
  void restart(int flag){};
  void field_ext_cls_test(const Equilibrium &equ, ParticleSpecies &pt);
  void field_ext_cls_calc_W(std::vector<double> &WWW, const Equilibrium &equ,
                            Particle &pt,
                            const std::vector<std::complex<double>> &phik_c,
                            const std::vector<int> &ntor1d);
  void field_ext_cls_calc_T(
      std::vector<std::complex<double>> &TTT, const Equilibrium &equ,
      ParticleSpecies &species, const std::vector<double> &partrad0,
      const std::vector<double> &parttheta0,
      const std::vector<double> &partphitor0,
      const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
      const std::vector<double> &partw0, const std::vector<double> &partfog0,
      std::vector<double> &draddt, std::vector<double> &dthetadt,
      std::vector<double> &dphitordt, std::vector<double> &dvpardt,
      std::vector<double> &dwdt, int icase,
      const std::vector<std::complex<double>> &phik_c,
      const std::vector<int> &ntor1d);
};
#endif // FIELDEXTCLS_H

void FieldExtCls::field_ext_cls_calc_W(
    std::vector<double> &WWW, const Equilibrium &equ, Particle &pt,
    const std::vector<std::complex<double>> &phik_c,
    const std::vector<int> &ntor1d) {
  int nsp = pt.getNsp();
  int nrad = this->nrad;
  int nthe = this->nthe;
  double rmax = this->radmax;
  double rmin = this->radmin;
  lenntor = ntor1d.size();
  WWW.assign(lenntor, 0.0);

  for (int isp = 0; isp < nsp; ++isp) {
    ParticleSpecies &species = pt.group.getSpecies(isp);
    int nptot = species.getNptot();
    if (nptot == 0) {
      std::cout << "Warning: ParticleSpecies " << isp
                << " has no particles, skipping." << std::endl;
      continue;
    }
    double mass = species.getMass();
    std::vector<double> xrad(nrad), wrad(nrad);
    std::vector<double> xthe(nthe), wthe(nthe);
    UtilMath::lgwt(nrad, rmin, rmax, xrad, wrad);
    double twopi = 2.0 * M_PI;
    UtilMath::lgwt(nthe, 0.0, twopi, xthe, wthe);
    // Nested loop for computation
    double volume = 0.0;
    double int_dens = 0.0;
    for (int fic = 0; fic < nrad; ++fic) {
      for (int fjc = 0; fjc < nthe; ++fjc) {
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
        field_cls_g2p2d1f_grad_complex(equ, phik_c, ntor1d, ptrad, ptthe, ptphi,
                                       dfdrad_c, dfdthe_c, dfdphi_c, 1, 0.0);

        volume += w12 * jaco3;
        int_dens += w12 * jaco3 * species.getdens1d(ptrad);

        double dens_rad = species.getdens1d(ptrad);
        double ptB2 = ptB * ptB;
        for (int itor = 0; itor < lenntor; ++itor) {
          // Calculate W for each toroidal mode
          double nabla_perp_phi =
              ptg11 * std::norm(dfdrad_c[itor]) +
              2.0 * ptg12 *
                  (dfdrad_c[itor].real() * dfdthe_c[itor].real() +
                   dfdrad_c[itor].imag() * dfdthe_c[itor].imag()) +
              ptg22 * std::norm(dfdthe_c[itor]);

          nabla_perp_phi *= mass * dens_rad / (ptB2);
          WWW[itor] += nabla_perp_phi * w12 * jaco3;
        }

      } // nthe
    }   // nrad
  }     // nsp
  for (int itor = 0; itor < lenntor; ++itor) {
    // integrate over toroidal angle
    WWW[itor] *= twopi;
  }
}
// calculate T for the particles
void FieldExtCls::field_ext_cls_calc_T(
    std::vector<std::complex<double>> &TTT, const Equilibrium &equ,
    ParticleSpecies &species, const std::vector<double> &partrad0,
    const std::vector<double> &parttheta0,
    const std::vector<double> &partphitor0,
    const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
    const std::vector<double> &partw0, const std::vector<double> &partfog0,
    std::vector<double> &draddt, std::vector<double> &dthetadt,
    std::vector<double> &dphitordt, std::vector<double> &dvpardt,
    std::vector<double> &dwdt, int icase,
    const std::vector<std::complex<double>> &phik_c,
    const std::vector<int> &ntor1d) {

  int nrad = this->nrad;
  int nthe = this->nthe;
  double rmax = this->radmax;
  double rmin = this->radmin;

  std::complex<double> zero_c(0.0, 0.0);
  TTT.assign(this->lenntor, zero_c); // Initialize TTT with complex

  int nptot = species.getNptot();
  double mass = species.getMass();
  int nptot_all = species.getNptotAll();
  double zcharge = species.getCharge();
  double nsonN = species.getNsonN();
  double Tem = species.getTem();
  int sp_id = species.getSpId();
  int ideltaf = species.getideltaf();
  double Cp2g = species.getCp2g();
  // ParticleCoords &coords = species.getCoords();
  // Particle loop for deposition
  for (int fic = 0; fic < nptot; ++fic) {
    double ptrad = partrad0[fic];
    double ptthe = parttheta0[fic];
    double ptphi = partphitor0[fic];
    double ptvpar = partvpar0[fic];
    double ptw = partw0[fic];
    double ptmu = partmu0[fic];
    double ptfog = partfog0[fic];
    double ptdraddt = draddt[fic];
    double ptdthedt = dthetadt[fic];
    double ptdphidt = dphitordt[fic];

    double ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe;
    double ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22, dfdpar;
    equ.calcBJgij(ptrad, ptthe, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe, ptBthe_ct,
                  ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22);

    double jaco3 = ptjaco2 * ptRR;

    std::vector<std::complex<double>> dfdrad_c, dfdthe_c, dfdphi_c;
    field_cls_g2p2d1f_grad_complex(equ, phik_c, ntor1d, ptrad, ptthe, ptphi,
                                   dfdrad_c, dfdthe_c, dfdphi_c, 1, 0.0);

    int lenntor = this->lenntor;
    std::complex<double> phase_factor(0.0, 0.0);
    for (int itor = 0; itor < lenntor; ++itor) {
      // Calculate the perturbation term
      std::complex<double> i(0.0, 1.0); 
      std::complex<double> phase_factor =
          std::exp(-i * static_cast<double>(ntor1d[itor]) * ptphi);
      TTT[itor] +=
          ptw * phase_factor *
          (ptdraddt * std::conj(dfdrad_c[itor]) + ptdthedt * std::conj(dfdthe_c[itor]));
    }

  } // nptot

  for (int itor = 0; itor < lenntor; ++itor) {
    TTT[itor] *= Cp2g * zcharge;
  }
}

// Implementation of init() method
void FieldExtCls::init(const Equilibrium &equ,
                       std::vector<ParticleSpecies> &pt) {

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
    std::cout << "----Initializing perturbations with nthe=" << nthe
              << ", nrad=" << nrad << std::endl;
    std::cout << "radmin=" << radmin << ", radmax=" << radmax
              << ", themin=" << themin << ", themax=" << themax << std::endl;
    std::cout << "drad=" << drad << ", dtheta=" << dtheta << std::endl;
    std::cout << "lenntor=" << lenntor << std::endl;
    for (size_t i = 0; i < ntor1d.size(); ++i) {
      std::cout << "ntor1d[" << i << "] = " << ntor1d[i] << std::endl;
    }
  }
  // 直接在函数中定义并初始化
  std::vector<double> rc1_arr(this->lenntor, 0.5);
  std::vector<double> amp_arr(this->lenntor, 0.1);
  std::vector<double> rwidth_arr(this->lenntor, 0.1);
  std::vector<double> omega_arr(this->lenntor, 1.0); // real frequency
  std::vector<int> mpoloidal_arr(this->lenntor);

  for (int i = 0; i < this->lenntor; ++i) {
    mpoloidal_arr[i] = 2 + i; // 例如：每个 toroidal mode 不同的 poloidal 模数
  }
  if (rank == 0) {
    std::cout
        << "Initializing perturbations with rc1_arr, amp_arr, rwidth_arr, "
           "mpoloidal_arr\n";
    for (int i = 0; i < this->lenntor; ++i) {
      std::cout << "itor=" << i << ", rc1=" << rc1_arr[i]
                << ", amp=" << amp_arr[i] << ", rwidth=" << rwidth_arr[i]
                << ", mpoloidal=" << mpoloidal_arr[i] << "\n";
    }
  }
  // give a MHD perturbation delta A_//=i/omega* \partial_// delta Phi
  for (int itor = 0; itor < this->lenntor; ++itor) {
    double rc1 = rc1_arr[itor];
    double amp = amp_arr[itor];
    double rwidth = rwidth_arr[itor];
    int mpoloidal = mpoloidal_arr[itor];
    double omega = omega_arr[itor];
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
          radial_part = amp * std::exp(-std::pow((r - rc1) / rwidth, 2));
        }

        // 角向复数部分
        std::complex<double> angular_part =
            std::exp(std::complex<double>(0, mpoloidal * theta));

        this->phik[idx] = radial_part * angular_part;

        // 计算对应的 A_//=i/omega * \partial_// delta Phi
        double ww = 0.0;
        if (i == 0 || i == nradfem - 1) {
          ww = 0.0; // 边界条件
        } else {
          double ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe;
          double ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22;
          equ.calcBJgij(r, theta, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                        ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12,
                        ptg22);
          ww =
              -(mpoloidal * ptBthe_ct + ntor1d[itor] * ptBphi_ct) / ptB / omega;
        }

        this->apark[idx] =
            std::complex<double>(ww, 0) * radial_part * angular_part;
      }
    }
  }

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
          field_cls_g2p2d1f_general(equ, phik, ntor1d, ptrad, ptthe, ptphi, ff,
                                    std::array<int, 3>{0, 0, 0}, 1, prho);

          if (igradient == 0) {
            field_cls_g2p2d1f_general(equ, phik, ntor1d, ptrad, ptthe, ptphi,
                                      dfdrad, std::array<int, 3>{1, 0, 0}, 1,
                                      prho);
            field_cls_g2p2d1f_general(equ, phik, ntor1d, ptrad, ptthe, ptphi,
                                      dfdthe, std::array<int, 3>{0, 1, 0}, 1,
                                      prho);
            field_cls_g2p2d1f_general(equ, phik, ntor1d, ptrad, ptthe, ptphi,
                                      dfdphi, std::array<int, 3>{0, 0, 1}, 1,
                                      prho);
          } else if (igradient == 1) {
            field_cls_g2p2d1f_grad(equ, phik, ntor1d, ptrad, ptthe, ptphi,
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

    field_cls_g2p2d1f_general(equ, phik, ntor1d, ptrad, ptthe, ptphi, ptf1d,
                              std::array<int, 3>{0, 0, 0}, 1, prho1d);
    std::vector<double> dfdrad1d, dfdthe1d,
        dfdphi1d; // Gradients output vectors
    field_cls_g2p2d1f_grad(equ, phik, ntor1d, ptrad, ptthe, ptphi, dfdrad1d,
                           dfdthe1d, dfdphi1d, 1, prho1d);

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
