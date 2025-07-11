#include "ParticleExtCls.h"

void ParticleExtCls::particle_ext_cls_dxvpardt123EM2d1f_2sp(
    Equilibrium &equ, FieldCls &fd,
    const std::vector<std::complex<double>> &phik,
    const std::vector<std::complex<double>> &apark,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp,
    const std::vector<ParticleCoords> &xv0,
    std::vector<std::vector<double>> &partmu0_allsp,
    std::vector<std::vector<double>> &partfog0_allsp,
    std::vector<ParticleCoords> &dxvdt) {

  if (xv0.size() != getNsp()) {
    std::cerr << "Error: xv0 size does not match ntor1d size." << std::endl;
    return;
  }
  if (partmu0_allsp.size() != getNsp()) {
    std::cerr << "Error: partmu0_allsp size does not match ntor1d size."
              << std::endl;
    return;
  }

  for (int fsc = 0; fsc < getNsp(); ++fsc) {

    particle_ext_cls_dxvpardt123EM2d1f(
        fsc, equ, fd, phik, apark, fd.ntor1d, amp, xv0[fsc].partrad,
        xv0[fsc].parttheta, xv0[fsc].partphitor, xv0[fsc].partvpar,
        partmu0_allsp[fsc], xv0[fsc].partw, partfog0_allsp[fsc],
        dxvdt[fsc].partrad, dxvdt[fsc].parttheta, dxvdt[fsc].partphitor,
        dxvdt[fsc].partvpar, dxvdt[fsc].partw);
  }
}
void ParticleExtCls::particle_ext_cls_dxvpardt123EM2d1f(
    int speciesIndex, Equilibrium &equ, FieldCls &fd,
    const std::vector<std::complex<double>> &phik,
    const std::vector<std::complex<double>> &apark,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp,
    const std::vector<double> &partrad0, const std::vector<double> &parttheta0,
    const std::vector<double> &partphitor0,
    const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
    const std::vector<double> &partw0, const std::vector<double> &partfog0,
    std::vector<double> &draddt, std::vector<double> &dthetadt,
    std::vector<double> &dphitordt, std::vector<double> &dvpardt,
    std::vector<double> &dwdt) {
  // Equivalent to calling the general function with EM case (case_id = 2)
  this->particle_ext_cls_dxvpardt123EMgeneral(
      speciesIndex, equ, fd, partrad0, parttheta0, partphitor0, partvpar0,
      partmu0, partw0, partfog0, draddt, dthetadt, dphitordt, dvpardt, dwdt,
      2, // EM case
      phik, apark, ntor1d, amp);
}

void ParticleExtCls::particle_ext_cls_dxvpardt123EMgeneral(
    int speciesIndex, Equilibrium &equ, FieldCls &fd,
    const std::vector<double> &partrad0, const std::vector<double> &parttheta0,
    const std::vector<double> &partphitor0,
    const std::vector<double> &partvpar0, const std::vector<double> &partmu0,
    const std::vector<double> &partw0, const std::vector<double> &partfog0,
    std::vector<double> &draddt, std::vector<double> &dthetadt,
    std::vector<double> &dphitordt, std::vector<double> &dvpardt,
    std::vector<double> &dwdt, int icase,
    const std::vector<std::complex<double>> &phik_c,
    const std::vector<std::complex<double>> &apark_c,
    const std::vector<int> &ntor1d,
    const std::vector<std::complex<double>> &amp) {
  // Precompute constant values
  ParticleSpecies &species = group.getSpecies(speciesIndex);
  const double rhoN = equ.rhoN;
  const double rhotN = equ.rhoN;
  const double Bref = equ.getBref();
  const double mass = species.getMass();
  const double zcharge = species.getCharge();
  const int nptot = species.getNptot();
  const int ngyro = species.getNgyro();
  const int ischeme_motion = species.getIschemeMotion();
  const int ideltaf = species.getideltaf();
  const int rank = species.getRank();
  double ffac = zcharge / mass;
  double cc1 = rhoN * mass * Bref / zcharge * std::sqrt(2.0);

  std::vector<double> pt1ddPdrad, pt1ddPdthe, pt1ddPdphi;
  std::vector<double> pt1ddAdrad, pt1ddAdthe, pt1ddAdphi;
  std::vector<double> pt1ddAhdrad, pt1ddAhdthe, pt1ddAhdphi;
  std::vector<double> pt1ddAsdt, pt1dAorAh, pt1dtmp;
  std::vector<double> rho1arr, ptB1d;

  // 2d1f push EM particle dynamics

  if (this->batch_g2p) {
    pt1ddPdrad.resize(nptot);
    pt1ddPdthe.resize(nptot);
    pt1ddPdphi.resize(nptot);
    pt1ddAdrad.resize(nptot);
    pt1ddAdthe.resize(nptot);
    pt1ddAdphi.resize(nptot);
    pt1ddAhdrad.resize(nptot);
    pt1ddAhdthe.resize(nptot);
    pt1ddAhdphi.resize(nptot);
    pt1ddAsdt.resize(nptot);
    pt1dAorAh.resize(nptot);
    pt1dtmp.resize(nptot);
    rho1arr.resize(nptot);
    ptB1d.resize(nptot);

    // 计算B场
    for (int fic = 0; fic < nptot; ++fic) {
      ptB1d[fic] = equ.getB(partrad0[fic], parttheta0[fic]);
    }

    // 计算rho
    for (int i = 0; i < nptot; ++i) {
      rho1arr[i] = cc1 * std::sqrt(partmu0[i] / ptB1d[i]);
    }

    { // 2D1F

      // Phi梯度
      fd.field_cls_g2p2d1f_grad(equ, phik_c, ntor1d, amp, partrad0, parttheta0,
                                partphitor0, pt1ddPdrad, pt1ddPdthe, pt1ddPdphi,
                                ngyro, rho1arr);

      // A梯度
      fd.field_cls_g2p2d1f_grad(equ, apark_c, ntor1d, amp, partrad0, parttheta0,
                                partphitor0, pt1ddAdrad, pt1ddAdthe, pt1ddAdphi,
                                ngyro, rho1arr);

      //   if (this->imixvar == 0)
      {
        // A for v pullback
        fd.field_cls_g2p2d1f_general(
            equ, apark_c, ntor1d, amp, partrad0, parttheta0, partphitor0,
            pt1dAorAh, std::array<int, 3>{0, 0, 1}, ngyro, rho1arr);

        // Ah梯度无需求
        std::fill(pt1ddAhdrad.begin(), pt1ddAhdrad.end(), 0.0);
        std::fill(pt1ddAhdthe.begin(), pt1ddAhdthe.end(), 0.0);
        std::fill(pt1ddAhdphi.begin(), pt1ddAhdphi.end(), 0.0);

        // E_parallel无需求
        std::fill(pt1ddAsdt.begin(), pt1ddAsdt.end(), 0.0);
      }
    }
  }

  // Initialization
  draddt.assign(nptot, 0.0);
  dthetadt.assign(nptot, 0.0);
  dphitordt.assign(nptot, 0.0);
  dvpardt.assign(nptot, 0.0);
  dwdt.assign(nptot, 0.0);
  //

  // Main particle loop
  for (int fic = 0; fic < nptot; ++fic) {
    double ptrad = partrad0[fic];
    double pttheta = parttheta0[fic];
    double ptphitor = partphitor0[fic];
    double ptvpar = partvpar0[fic];
    double ptmu = partmu0[fic];
    double ptw = partw0[fic];
    double ptfog = partfog0[fic];

    double ptBthe_ct, ptBphi_ct;
    double ptBrad_co, ptBthe_co, ptBphi_co;
    double ptB, ptRR, ptZZ;
    double ptdBdrad, ptdBdthe;
    double ptjaco2, ptFF, ptg11, ptg12, ptg22;
    double ptdAdrad, ptdAdthe, ptdAdphi;
    double ptA;
    double ptEB_dphidt;

    // Magnetic field calculations
    equ.calcBJgij(ptrad, pttheta, ptRR, ptZZ, ptB, ptdBdrad, ptdBdthe,
                  ptBthe_ct, ptBphi_ct, ptjaco2, ptFF, ptg11, ptg12, ptg22);
    double ptjaco3 = ptjaco2 * ptRR;

    double JB_inv = 1.0 / (ptB * ptjaco3);
    // --calc. Gyro
    double rho1 = cc1 * std::sqrt(ptmu / ptB);

    if (this->batch_g2p) {
      ptdAdrad = pt1ddAdrad[fic];
      ptdAdthe = pt1ddAdthe[fic];
      ptdAdphi = pt1ddAdphi[fic];
    }

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

      // 对bstar进行归一化
      bstar_rad_ct = bstar_rad_ct / BBstar_abs;
      bstar_the_ct = bstar_the_ct / BBstar_abs;
      bstar_phi_ct = bstar_phi_ct / BBstar_abs;

      // 对bstarAs进行归一化
      bstarAs_rad_ct = bstarAs_rad_ct / BBstar_abs;
      bstarAs_the_ct = bstarAs_the_ct / BBstar_abs;
      bstarAs_phi_ct = bstarAs_phi_ct / BBstar_abs;

      // 计算bstar1 = bstarAs - bstar, 是指
      // bstar1 = $+\nabla\langle\delta A_\|^{\rm s}\rangle\times{\boldsymbol
      // b}/B_\|^*$
      bstar1_rad_ct = bstarAs_rad_ct - bstar_rad_ct;
      bstar1_the_ct = bstarAs_the_ct - bstar_the_ct;
      bstar1_phi_ct = bstarAs_phi_ct - bstar_phi_ct;
      // bstarAs 包含 As项，是总的b*
      bstar_minus_b_rad_ct = bstarAs_rad_ct;
      bstar_minus_b_the_ct = bstarAs_the_ct - ptBthe_ct / ptB;
      bstar_minus_b_phi_ct = bstarAs_phi_ct - ptBphi_ct / ptB;
    }

    //// === 调用场（field）插值 ===
    fd.field_cls_g2p2d1f_general(equ, apark_c, ntor1d, amp, ptrad, pttheta,
                                 ptphitor, ptA, std::array<int, 3>{0, 0, 0},
                                 ngyro, rho1);

    // ptvparph = ptvpar-zcharge/mass*ptAh;
    double ptvparph = ptvpar; // Ah=0, only As

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
    double ptCEphi = ptCE * ptFF / (ptjaco3 * std::pow(ptB, 2));

    if (this->batch_g2p) {
      ptdPdrad = pt1ddPdrad[fic];
      ptdPdthe = pt1ddPdthe[fic];
      ptdPdphi = pt1ddPdphi[fic];
      // pttmp = pt1dtmp[fic];  //
    } else {

      fd.field_cls_g2p2d1f_grad(equ, phik_c, ntor1d, amp, ptrad, pttheta,
                                ptphitor, ptdPdrad, ptdPdthe, ptdPdphi, ngyro,
                                rho1);
    }

    double ptdPdpar = (ptBthe_ct * ptdPdthe + ptBphi_ct * ptdPdphi) / ptB;
    double ptdAdpar = (ptBthe_ct * ptdAdthe + ptBphi_ct * ptdAdphi) / ptB;

    double ptdPAdrad = ptdPdrad - ptvparph * ptdAdrad;
    double ptdPAdthe = ptdPdthe - ptvparph * ptdAdthe;
    double ptdPAdphi = ptdPdphi - ptvparph * ptdAdphi;
    double ptdPAdpar = ptdPdpar - ptvparph * ptdAdpar;

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
    // ====Weight
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
          // finalize_mpi_petsc();
        }
      }

      double ptTem = species.getTem1d(ptrad);
      double ptTfac = mass / ptTem * (ptvpar * ptvpar + 2.0 * ptmu * ptB) - 1.5;

      dwdt[fic] +=
          (2.0 * ptvpar * ptdvpardt_tot * mass / ptTem -
           pt_draddt * (species.getdlndensdr1d(ptrad) +
                        ptTfac * species.getdlnTemdr1d(ptrad) -
                        2.0 * mass * ptmu / ptTem *
                            ptdBdrad) // Note 2/T is from TN = mN*vN^2/2
           + 2.0 * pt_dthedt * mass * ptmu / ptTem * ptdBdthe) *
          ptwfac;

      if (this->isrcsnk != 0) {
        dwdt[fic] += ptw * species.get_fsrcsnk(ptrad);
      }
    }
  } // End of particle loop
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
    const auto &coords = species.getCoords(); // 先拿到引用，避免多次调用
    int nptot = species.getNptot();
    ParticleCoords dxvdt(nptot), dxv_sum(nptot), xv0(nptot);
    species.particle_cls_track_init(iset_track, Bax);

    // Initialize coordinates
    dxvdt.setZero();
    dxv_sum.setZero();
    xv0.setZero();

    const std::vector<double> &partrad = coords.partrad;
    const std::vector<double> &parttheta = coords.parttheta;
    const std::vector<double> &partphitor = coords.partphitor;
    const std::vector<double> &partvpar = coords.partvpar;
    const std::vector<double> &partw = coords.partw;

    const std::vector<double> &partmu = species.partmu;
    const std::vector<double> &partfog = species.partfog;
    //   size_t nptot = nptot;
    auto zero_complex = std::complex<double>(0.0, 0.0);

    int ntotfem2d1f = fd.getNtotfem2d1f();
    std::vector<std::complex<double>> fk1(ntotfem2d1f, zero_complex),
        fk2(ntotfem2d1f, zero_complex);
    std::vector<std::complex<double>> amp_test(fd.lenntor, 1.0);
    //   Only test field can be used by particle

    this->particle_ext_cls_dxvpardt123EM2d1f(
        i, equ, fd, fk1, fk2, fd.ntor1d, amp_test, partrad, parttheta,
        partphitor, partvpar, partmu, partw, partfog, dxvdt.partrad,
        dxvdt.parttheta, dxvdt.partphitor, dxvdt.partvpar, dxvdt.partw);
  }
}
