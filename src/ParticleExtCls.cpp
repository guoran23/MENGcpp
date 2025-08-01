#include "ParticleExtCls.h"

void ParticleExtCls::particle_ext_cls_dxvpardt123EM2d1f_2sp(
    const Equilibrium &equ, const FieldCls &fd,
    const std::vector<std::complex<double>> &phik,
    const std::vector<std::complex<double>> &apark,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp,
    std::vector<ParticleCoords> &dxvdt,
    std::vector<std::complex<double>> &TTT_allsp) {

  int nsp = this->getNsp();
  int lenntor = ntor1d.size();
  constexpr std::complex<double> zero_c(0.0, 0.0);
  TTT_allsp.resize(lenntor, zero_c);

  for (int fsc = 0; fsc < getNsp(); ++fsc) {

    ParticleSpecies &species = this->group.getSpecies(fsc);
    ParticleCoords &coords = species.getCoords();
    std::vector<std::complex<double>> T_onesp(lenntor, zero_c);

    particle_ext_cls_dxvpardt123EM2d1f(
        fsc, equ, fd, phik, apark, fd.ntor1d, amp, coords.partrad,
        coords.parttheta, coords.partphitor, coords.partvpar, species.partmu,
        coords.partw, species.partfog, dxvdt[fsc].partrad, dxvdt[fsc].parttheta,
        dxvdt[fsc].partphitor, dxvdt[fsc].partvpar, dxvdt[fsc].partw, T_onesp);

    add_vectors(TTT_allsp, T_onesp);
  }
  // print_complex_vector("TTT_allsp = ", TTT_allsp);
}
void ParticleExtCls::particle_ext_cls_dxvpardt123EM2d1f(
    const int speciesIndex, const Equilibrium &equ, const FieldCls &fd,
    const std::vector<std::complex<double>> &phik,
    const std::vector<std::complex<double>> &apark,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp, std::vector<double> &partrad0,
    std::vector<double> &parttheta0, std::vector<double> &partphitor0,
    std::vector<double> &partvpar0, std::vector<double> &partmu0,
    std::vector<double> &partw0, std::vector<double> &partfog0,
    std::vector<double> &draddt, std::vector<double> &dthetadt,
    std::vector<double> &dphitordt, std::vector<double> &dvpardt,
    std::vector<double> &dwdt, std::vector<std::complex<double>> &TTT_onesp) {
  // Equivalent to calling the general function with EM case (case_id = 2)
  this->particle_ext_cls_dxvpardt123EMgeneral(
      speciesIndex, equ, fd, partrad0, parttheta0, partphitor0, partvpar0,
      partmu0, partw0, partfog0, draddt, dthetadt, dphitordt, dvpardt, dwdt,
      2, // EM case
      phik, apark, ntor1d, amp, TTT_onesp);
}

void ParticleExtCls::particle_ext_cls_dxvpardt123EMgeneral(
    const int speciesIndex, const Equilibrium &equ, const FieldCls &fd,
    std::vector<double> &partrad0, std::vector<double> &parttheta0,
    std::vector<double> &partphitor0, std::vector<double> &partvpar0,
    std::vector<double> &partmu0, std::vector<double> &partw0,
    std::vector<double> &partfog0, std::vector<double> &draddt,
    std::vector<double> &dthetadt, std::vector<double> &dphitordt,
    std::vector<double> &dvpardt, std::vector<double> &dwdt, int icase,
    const std::vector<std::complex<double>> &phik_c,
    const std::vector<std::complex<double>> &apark_c,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp,
    std::vector<std::complex<double>> &TTT_onesp) {
  // Precompute constant values
  ParticleSpecies &species = group.getSpecies(speciesIndex);
  const double rhoN = equ.rhoN;
  const double rhotN = equ.rhoN;
  const double Bref = equ.getBref();
  const double mass = species.getMass();
  const double Tem = species.getTem();
  const double vts_max = 5.0 * std::sqrt(Tem / mass);
  const double B_axis = equ.Bmaxis_adhoc;
  const double zcharge = species.getCharge();
  const double Cp2g = species.getCp2g();
  const int nptot = species.getNptot();
  const int ngyro = species.getNgyro();
  const int ischeme_motion = species.getIschemeMotion();
  const int ideltaf = species.getideltaf();
  const int rank = species.getRank();
  double ffac = zcharge / mass;
  double cc1 = rhoN * mass * Bref / zcharge * std::sqrt(2.0);
  int lenntor = ntor1d.size();

  // 2d1f push EM particle dynamics
  // Initialization
  draddt.assign(nptot, 0.0);
  dthetadt.assign(nptot, 0.0);
  dphitordt.assign(nptot, 0.0);
  dvpardt.assign(nptot, 0.0);
  dwdt.assign(nptot, 0.0);
  //
  std::complex<double> zero_c(0.0, 0.0);
  TTT_onesp.assign(lenntor, zero_c);

  // Main particle loop
  for (int fic = 0; fic < nptot; ++fic) {
    double ptrad = partrad0[fic];
    double pttheta = parttheta0[fic];
    double ptphitor = partphitor0[fic];
    double ptvpar = partvpar0[fic];
    double ptmu = partmu0[fic];
    double ptw = partw0[fic];
    double ptfog = partfog0[fic];
    std::vector<std::complex<double>> T_onepar(lenntor, zero_c);

    if (ptrad < fd.radmin || ptrad > fd.radmax) {
      draddt[fic] = 0.0;
      dthetadt[fic] = 0.0;
      dphitordt[fic] = 0.0;
      dvpardt[fic] = 0.0;
      dwdt[fic] = 0.0;
      partw0[fic] = 0.0;
      // std::cout << "out boundary fic=" << fic << " " << ptrad << " " <<
      // pttheta
      //           << " " << ptphitor << " " << ptvpar << " " << ptmu << " " <<
      //           ptw
      //           << " " << ptfog << std::endl;

      continue; // Skip particles outside the radial range
    }
    if (std::abs(ptw) > 1.0) {
      std::cout << "weight too large "
                << "fic=" << fic << ", weight=" << ptw << std::endl;
      partw0[fic] = 0.0;
      continue;
    }
    if (std::abs(ptvpar) > vts_max) {
      std::cout << "v_|| too large "
                << "fic=" << fic << ", ptvpar=" << ptvpar << std::endl;
      partw0[fic] = 0.0;
      partvpar0[fic] = 0.0;
      continue;
    }
    if (std::abs(std::sqrt(2.0 * ptmu * B_axis)) > vts_max) {
      std::cout << "v_perp too large "
                << "fic=" << fic << ", ptmu=" << ptmu << std::endl;
      partw0[fic] = 0.0;
      partmu0[fic] = 0.0;
      continue;
    }

    double ptBthe_ct, ptBphi_ct;
    double ptBrad_co, ptBthe_co, ptBphi_co;
    double ptB, ptRR, ptZZ;
    double ptdBdrad, ptdBdthe;
    double ptjaco2, ptFF, ptg11, ptg12, ptg22;
    double ptdAdrad, ptdAdthe, ptdAdphi;
    double ptEB_dphidt;

    // Magnetic field calculations
    equ.calcBJgij(ptrad, pttheta, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                  ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22);
    double ptjaco3 = ptjaco2 * ptRR;

    double JB_inv = 1.0 / (ptB * ptjaco3);
    // --calc. Gyro
    double rho1 = cc1 * std::sqrt(ptmu / ptB);

    fd.field_cls_g2p2d1f_grad(equ, apark_c, ntor1d, amp, ptrad, pttheta,
                              ptphitor, ptdAdrad, ptdAdthe, ptdAdphi, ngyro,
                              rho1);

    // Placeholder curlB, Bstar, etc.
    double bstar_rad_ct = 0.0, bstar_the_ct = 0.0, bstar_phi_ct = 0.0,
           bstarAs_rad_ct = 0.0, bstarAs_the_ct = 0.0, bstarAs_phi_ct = 0.0,
           bstar1_rad_ct = 0.0, bstar1_the_ct = 0.0, bstar1_phi_ct = 0.0;
    double bstar_minus_b_rad_ct = 0.0, bstar_minus_b_the_ct = 0.0,
           bstar_minus_b_phi_ct = 0.0;
    double BBstar_abs = 0.0;

    if (ischeme_motion == 1) {
      // 方案1：简单Bstar计算
      bstar1_rad_ct = 0.0;
      bstar1_the_ct = 0.0;
      bstar1_phi_ct = 0.0;

      BBstar_abs = ptB;

      bstar_rad_ct = 0.0;
      bstar_the_ct = ptBthe_ct / ptB;
      bstar_phi_ct = ptBphi_ct / ptB;

      bstarAs_rad_ct = bstar_rad_ct;
      bstarAs_the_ct = bstar_the_ct;
      bstarAs_phi_ct = bstar_phi_ct;
    } else if (ischeme_motion == 2) {
      // 方案2：包含curlB项的完整计算

      // 获取共变基矢量B
      ptBrad_co = equ.getBrad_co(ptrad, pttheta);
      ptBthe_co = equ.getBthe_co(ptrad, pttheta);
      ptBphi_co = equ.getBphi_co(ptrad, pttheta);

      // Mishchenko 23JCP Eq 2.7的系数
      double facbstar = Bref * rhotN * mass / zcharge * ptvpar;

      // bstar1 = $(m_s/q_s)u_\|\nabla\times{\boldsymbol b}/B_\|^*$
      bstar1_rad_ct = facbstar * equ.getcurlbrad_ct(ptrad, pttheta);
      bstar1_the_ct = facbstar * equ.getcurlbthe_ct(ptrad, pttheta);
      bstar1_phi_ct = facbstar * equ.getcurlbphi_ct(ptrad, pttheta);

      // 计算|B*|
      BBstar_abs =
          ptB + (ptBrad_co * bstar1_rad_ct + ptBthe_co * bstar1_the_ct +
                 ptBphi_co * bstar1_phi_ct) /
                    ptB;

      // ${\boldsymbol b}^*_0={\boldsymbol
      // b}+(m_s/q_s)u_\|\nabla\times{\boldsymbol b}/B_\|^*$
      bstar_rad_ct = +bstar1_rad_ct;
      bstar_the_ct = ptBthe_ct + bstar1_the_ct;
      bstar_phi_ct = ptBphi_ct + bstar1_phi_ct;

      // bstarAs, 总的b* = b^*_0 + As修正项
      bstarAs_rad_ct =
          bstar_rad_ct + JB_inv * (ptdAdphi * ptBthe_co - ptdAdthe * ptBphi_co);
      bstarAs_the_ct =
          bstar_the_ct + JB_inv * (ptdAdrad * ptBphi_co - ptdAdphi * ptBrad_co);
      bstarAs_phi_ct =
          bstar_phi_ct + JB_inv * (ptdAdthe * ptBrad_co - ptdAdrad * ptBthe_co);

      // 对bstar进行归一化  ${\boldsymbol b}^*_0$
      bstar_rad_ct = bstar_rad_ct / BBstar_abs;
      bstar_the_ct = bstar_the_ct / BBstar_abs;
      bstar_phi_ct = bstar_phi_ct / BBstar_abs;

      // 对bstarAs进行归一化
      // bstarAs 包含 As项，是总的b*
      bstarAs_rad_ct = bstarAs_rad_ct / BBstar_abs;
      bstarAs_the_ct = bstarAs_the_ct / BBstar_abs;
      bstarAs_phi_ct = bstarAs_phi_ct / BBstar_abs;

      // 计算bstar1 = bstarAs - bstar, 是指
      // bstar1 = $+\nabla\langle\delta A_\|^{\rm s}\rangle\times{\boldsymbol
      // b}/B_\|^*$
      bstar1_rad_ct = bstarAs_rad_ct - bstar_rad_ct;
      bstar1_the_ct = bstarAs_the_ct - bstar_the_ct;
      bstar1_phi_ct = bstarAs_phi_ct - bstar_phi_ct;
      // b*-b, b equilibrium part
      bstar_minus_b_rad_ct = bstarAs_rad_ct;
      bstar_minus_b_the_ct = bstarAs_the_ct - ptBthe_ct / ptB;
      bstar_minus_b_phi_ct = bstarAs_phi_ct - ptBphi_ct / ptB;
    }

    // ptvparph = ptvpar-zcharge/mass*ptAh;
    // double ptvparph = ptvpar; // Ah=0, only As

    // === Part 1: Equilibrium motion ===
    // d(x0 + x_d) / dt and dvpar0/dt
    // ---- 1. dx/dt contribution from v_par ----
    double ptvpar0_draddt = ptvpar * bstar_rad_ct;
    double ptvpar0_dthedt = ptvpar * bstar_the_ct;
    double ptvpar0_dphidt = ptvpar * bstar_phi_ct;

    if (species.getvPar0() != 0) {
      draddt[fic] += ptvpar0_draddt;
      dthetadt[fic] += ptvpar0_dthedt;
      dphitordt[fic] += ptvpar0_dphidt;
    }

    // ---- 2. dx/dt contribution from v_drift ----
    double ptC1;
    if (ischeme_motion == 1) {
      ptC1 = mass / zcharge * rhotN * (ptvpar * ptvpar + ptmu * ptB) * Bref /
             std::pow(ptB, 3);
    } else if (ischeme_motion == 2) {
      ptC1 = mass / zcharge * rhotN * (ptmu * ptB) * Bref /
             (BBstar_abs * ptB * ptB);
    }

    double ptvd_draddt = +ptC1 * ptFF * ptdBdthe / ptjaco3;
    double ptvd_dthedt = -ptC1 * ptFF * ptdBdrad / ptjaco3;
    double ptvd_dphidt = +ptC1 * equ.getdpsidr(ptrad) *
                         (ptdBdrad * ptg11 + ptdBdthe * ptg12) / (ptRR * ptRR);

    if (species.getvD() != 0) {
      draddt[fic] += ptvd_draddt;
      dthetadt[fic] += ptvd_dthedt;
      dphitordt[fic] += ptvd_dphidt;
    }

    // ---- 3. dvpar/dt from mirror force ----
    double ptdvpardt0 =
        -ptmu * (bstar_rad_ct * ptdBdrad + bstar_the_ct * ptdBdthe);

    if (species.getvMirror() != 0) {
      dvpardt[fic] += ptdvpardt0;
    }

    // ==== Part 2
    double ptCE, ptdPdrad, ptdPdthe, ptdPdphi;
    if (ischeme_motion == 1) {
      ptCE = rhotN * Bref;
    } else if (ischeme_motion == 2) {
      ptCE = rhotN * Bref * ptB / BBstar_abs;
    }

    double ptCEthe = ptCE * equ.getdpsidr(ptrad) / std::pow(ptB * ptRR, 2);
    double ptCEphi = ptCE * ptFF / (ptjaco3 * ptB * ptB);

    fd.field_cls_g2p2d1f_grad(equ, phik_c, ntor1d, amp, ptrad, pttheta,
                              ptphitor, ptdPdrad, ptdPdthe, ptdPdphi, ngyro,
                              rho1);

    double ptdPdpar = (ptBthe_ct * ptdPdthe + ptBphi_ct * ptdPdphi) / ptB;
    double ptdAdpar = (ptBthe_ct * ptdAdthe + ptBphi_ct * ptdAdphi) / ptB;

    double ptdPAdrad = ptdPdrad - ptvpar * ptdAdrad;
    double ptdPAdthe = ptdPdthe - ptvpar * ptdAdthe;
    double ptdPAdphi = ptdPdphi - ptvpar * ptdAdphi;
    double ptdPAdpar = ptdPdpar - ptvpar * ptdAdpar;

    double ptEB_draddt = -ptCEthe * ptg11 * ptdPAdphi + ptCEphi * ptdPAdthe;
    double ptEB_dthedt = -ptCEthe * ptg12 * ptdPAdphi - ptCEphi * ptdPAdrad;

    if (ideltaf != 2 && species.getvExB() != 0) {
      ptEB_dphidt = ptCEthe * (ptg11 * ptdPAdrad + ptg12 * ptdPAdthe);

      draddt[fic] = draddt[fic] + ptEB_draddt;
      dthetadt[fic] = dthetadt[fic] + ptEB_dthedt;
      dphitordt[fic] = dphitordt[fic] + ptEB_dphidt;
    }
    // ====Part 3: dv||/dt
    double ptdvpardt1 = -ffac * (bstar_minus_b_rad_ct * ptdPdrad +
                                 bstar_minus_b_the_ct * ptdPdthe +
                                 bstar_minus_b_phi_ct * ptdPdphi);
    // ! As term in dvpar1/dt
    double ptCvpar = -rhotN * Bref * ptmu / (ptB * BBstar_abs * ptjaco3);
    ptdvpardt1 = ptdvpardt1 + ptCvpar * (ptBphi_ct * ptdBdthe * ptdAdrad -
                                         ptBphi_ct * ptdBdrad * ptdAdthe +
                                         ptBthe_ct * ptdBdrad * ptdAdphi);
    // ! Add to GC trajectory
    if (ideltaf != 2 && species.getvEpar() != 0) {
      dvpardt[fic] = dvpardt[fic] + ptdvpardt1;
    }
    // ====Weight====
    double ptwfac, ptdvpardt_tot;
    double pt_draddt, pt_dthedt, pt_dphidt;

    if (ideltaf == 1 || ideltaf == 3) { // 1:NL deltaf; 3:control variate
      ptwfac = ptfog - ptw;             // ptfog at t
    } else if (ideltaf == 2) {          // 2: linear df
      ptwfac = ptfog;
    }

    if ((ideltaf == 1 || ideltaf == 2 || ideltaf == 3) && this->idwdt != 0) {
      // 1:NL deltaf; 2:linear df

      if (ideltaf == 1 || ideltaf == 2) {
        // 1:NL deltaf; 2:linear df

        if (this->neocls == 0) {
          ptdvpardt_tot = ptdvpardt1;
          pt_draddt = ptEB_draddt;
          pt_dthedt = ptEB_dthedt;
        } else if (this->neocls == 1) {
          ptdvpardt_tot = ptdvpardt0;
          pt_draddt = ptvd_draddt;
          pt_dthedt = ptvd_dthedt + ptvpar0_dthedt;
        } else if (this->neocls == 2) {
          ptdvpardt_tot = ptdvpardt0 + ptdvpardt1;
          pt_draddt = ptEB_draddt + ptvd_draddt;
          pt_dthedt = ptEB_dthedt + ptvd_dthedt + ptvpar0_dthedt;
        } else if (this->neocls == 3) {
          ptdvpardt_tot = 0.0;
          pt_draddt = 0.0;
          pt_dthedt = 0.0;
        } else {
          if (rank == 0) {
            std::cerr << "**** error: wrong neocls value (0,1,2,3)"
                      << std::endl;
          }
          // finalize_mpi_petsc();
        }

      } else if (ideltaf == 3) {
        if (this->neocls == 0) { // note the sign '-'
          ptdvpardt_tot = -ptdvpardt0;
          pt_draddt = -ptvd_draddt;
          pt_dthedt = -ptvd_dthedt - ptvpar0_dthedt;
        } else if (this->neocls == 1) {
          ptdvpardt_tot = -ptdvpardt1;
          pt_draddt = -ptEB_draddt;
          pt_dthedt = -ptEB_dthedt;
        } else if (this->neocls == 2) {
          ptdvpardt_tot = 0.0;
          pt_draddt = 0.0;
          pt_dthedt = 0.0;
        } else if (this->neocls == 3) {
          ptdvpardt_tot = -ptdvpardt0 - ptdvpardt1;
          pt_draddt = -ptEB_draddt - ptvd_draddt;
          pt_dthedt = -ptEB_dthedt - ptvd_dthedt - ptvpar0_dthedt;
        } else {
          if (rank == 0) {
            std::cerr << "**** error: wrong neocls value (0,1,2,3)"
                      << std::endl;
          }
        }
      }

      double ptTem = species.getTem1d(ptrad);
      double ptTfac = mass / ptTem * (ptvpar * ptvpar + 2.0 * ptmu * ptB) - 1.5;
      // Note 2/T is from TN = mN*vN^2/2
      dwdt[fic] += (2.0 * ptvpar * ptdvpardt_tot * mass / ptTem -
                    pt_draddt * (species.getdlndensdr1d(ptrad) +
                                 ptTfac * species.getdlnTemdr1d(ptrad) -
                                 2.0 * mass * ptmu / ptTem * ptdBdrad) +
                    2.0 * pt_dthedt * mass * ptmu / ptTem * ptdBdthe) *
                   ptwfac;

      if (this->isrcsnk != 0) {
        dwdt[fic] += ptw * species.get_fsrcsnk(ptrad);
      }
    }
    // calc T
    if (std::abs(ptvd_draddt) > vts_max || std::abs(ptvd_dthedt) > vts_max) {
      ptvd_draddt = 0.0;
      ptvd_dthedt = 0.0;
      std::cout << "vd too large: ptvd_rad=" << ptvd_draddt
                << ", ptvd_theta=" << ptvd_dthedt << std::endl;
    }
    // T_onepar =
    //     calc_T_onePar(equ, fd, ptrad, pttheta, ptphitor, ptvpar, ptw,
    //                   ptvd_draddt + ptEB_draddt, ptvd_dthedt + ptEB_dthedt,
    //                   ptvd_dphidt + ptEB_dphidt, phik_c, ntor1d, amp);
    T_onepar = calc_T_onePar(equ, fd, ptrad, pttheta, ptphitor, ptvpar, ptw,
                             ptvd_draddt, ptvd_dthedt, ptvd_dphidt, phik_c,
                             ntor1d, amp);
    add_vectors(TTT_onesp, T_onepar);

    // std::ofstream out("T_onePar_output.txt", std::ios::app); // append 模式
    // if (out.is_open()) {
    //   out << std::imag(T_onepar[0]) << std::endl;
    //   out.close();
    // }

  } // End of particle loop

  for (int itor = 0; itor < lenntor; ++itor) {
    TTT_onesp[itor] *= Cp2g * zcharge;
  }
}

// calculate T for one particle
std::vector<std::complex<double>> ParticleExtCls::calc_T_onePar(
    const Equilibrium &equ, const FieldCls &fd, const double &partrad,
    const double &parttheta, const double &partphitor, const double &partvpar,
    const double &partw, const double &vd_rad, const double &vd_the,
    const double &vd_phi, const std::vector<std::complex<double>> &phik_c,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp) {

  int lenntor = ntor1d.size();
  std::vector<std::complex<double>> TTT;
  constexpr std::complex<double> zero_c(0.0, 0.0);
  TTT.assign(lenntor, zero_c); // Initialize TTT with complex zeros

  std::vector<std::complex<double>> dfdrad_c, dfdthe_c, dfdphi_c;
  fd.field_cls_g2p2d1f_grad_complex(equ, phik_c, ntor1d, amp, partrad,
                                    parttheta, partphitor, dfdrad_c, dfdthe_c,
                                    dfdphi_c, 1, 0.0);

  for (int itor = 0; itor < lenntor; ++itor) {
    // Calculate the perturbation term
    constexpr std::complex<double> i_c(0.0, 1.0);
    std::complex<double> phase_factor(0.0, 0.0);
    // phase_factor =
    //     std::exp(-i_c * static_cast<double>(ntor1d[itor]) * partphitor);
    TTT[itor] = partw * (vd_rad * std::conj(dfdrad_c[itor]) +
                         vd_the * std::conj(dfdthe_c[itor]));

    // if (std::norm(TTT[itor]) > 1e-16) {
    //   std::cout << "weight=" << partw << ", vd_rad=" << vd_rad
    //             << ", vd_the=" << vd_the << ", dfdrad_c =" << dfdrad_c[itor]
    //             << ", dfdthe_c =" << dfdthe_c[itor] << std::endl;
    // }
  }
  return TTT;
}

void ParticleExtCls::particle_ext_cls_test(Equilibrium &equ, FieldCls &fd,
                                           int irk, int nrun, double dt_o_Ttr,
                                           int iset_track) {
  // Working variables

  int fic;
  double dt;

  // Set passing particles
  double Bax = abs(equ.Bmaxis);

  // Start loop

  for (size_t i = 0; i < getNsp(); i++) {
    /* code */
    ParticleSpecies &species = group.getSpecies(i);
    if (species.getRank() == 0) {
      std::cout << "--Particle_ext_cls_test: --Only test field can be used by "
                   "particle"
                << std::endl;
    }
    // Get species coordinates
    auto &coords = species.getCoords(); // 先拿到引用，避免多次调用
    int nptot = species.getNptot();
    ParticleCoords dxvdt(nptot), dxv_sum(nptot), xv0(nptot);
    species.particle_cls_track_init(iset_track, Bax);

    // Initialize coordinates
    dxvdt.setZero();
    dxv_sum.setZero();
    xv0.setZero();

    std::vector<double> &partrad = coords.partrad;
    std::vector<double> &parttheta = coords.parttheta;
    std::vector<double> &partphitor = coords.partphitor;
    std::vector<double> &partvpar = coords.partvpar;
    std::vector<double> &partw = coords.partw;
    std::vector<double> &partmu = species.partmu;
    std::vector<double> &partfog = species.partfog;
    //   size_t nptot = nptot;
    auto zero_complex = std::complex<double>(0.0, 0.0);

    int ntotfem2d1f = fd.getNtotfem2d1f();
    std::vector<std::complex<double>> fk1(ntotfem2d1f, zero_complex),
        fk2(ntotfem2d1f, zero_complex);
    std::vector<std::complex<double>> amp_test(fd.lenntor, 1.0);
    std::vector<std::complex<double>> TTT(fd.lenntor, zero_complex);
    //   Only test field can be used by particle

    this->particle_ext_cls_dxvpardt123EM2d1f(
        i, equ, fd, fk1, fk2, fd.ntor1d, amp_test, partrad, parttheta,
        partphitor, partvpar, partmu, partw, partfog, dxvdt.partrad,
        dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, dxvdt.partw, TTT);
  }
}
