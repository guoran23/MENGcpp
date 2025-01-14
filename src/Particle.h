#ifndef PARTICLE_H
#define PARTICLE_H

#include "Equilibrium.h"
#include "Field.h"
#include <cmath>
#include <functional> // std::plus, std::multiplies
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

class ParticleCoords {
private:
  int nptot;
  std::vector<double> partrad, parttheta, partphitor;
  std::vector<double> partvpar;
  std::vector<double> partw;

public:
  // Constructor
  ParticleCoords(int total_particles)
      : nptot(total_particles), partrad(total_particles),
        parttheta(total_particles), partphitor(total_particles),
        partvpar(total_particles), partw(total_particles) {}

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
  int np;
  std::string species_name;
  double mass = 1.0;
  double zcharge = -1.0;
  ParticleCoords coords;

public:
  // Constructor
  ParticleSpecies(size_t nparticles, double m, double q,
                  const std::string &name = "unknown")
      : np(nparticles), coords(nparticles), mass(m), zcharge(q),
        species_name(name) {}

  // Accessors
  double getMass() const { return mass; }
  double getCharge() const { return zcharge; }
  ParticleCoords &getCoords() { return coords; }
  const ParticleCoords &getCoords() const { return coords; }

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
  std::shared_ptr<ParticleSpecies> getSpecies(size_t index) {
    if (index < speciesList.size()) {
      return speciesList[index];
    }
    return nullptr; // Return nullptr if index is out of bounds
  }

  int getSpeciesCount() const { return nsp; }

  // 打印所有粒子种类信息
  void print() const {
    std::cout << "Particle Group Info:\n";
    for (const auto &species : speciesList) {
      species->print();
      std::cout << "-----------------------------\n";
    }
  }
};

// class Particle {
// private:
//   int nsp;
//   ParticleGroup group;

// public:
//   int ischeme_motion, ngyro, irandom_gy;
//   int imixvar = 1;
//   int ibc_particle = 0, iset_zerosumw;
//   bool irestart = false;
//   std::string species_name;
//   int sp_id = 1;
//   int ideltaf = 2;
//   double phitormin = 0.0, phitormax, phitorwid;
//   int Nphimult = 2;
//   double nsonN = 1.0;
//   double Tem = 1.0;
//   int nptot;
//   long long nptot_all = 1000000;
//   double vparmaxovt = 5.0;
//   double vparmaxovN;
//   double mumaxovN2;
//   double vtovN;
//   std::string profilename;
//   std::vector<double> dens_coef1d = {0.49123, 0.298228, 0.198739, 0.521298};
//   std::vector<double> Tem_coef1d = {0.49123, 0.298228, 0.198739, 0.521298};
//   int idensprof;
//   int iTemprof;
//   // Assuming bspline_1d is a class or struct defined elsewhere
//   // bspline_1d dens_bsp, Tem_bsp;
//   int load_can, ischeme_load;
//   double Bax;
//   double rmin = 0.2, rmax = 0.8;
//   double Stot, Vtot, dens_mean, dens_intg;
//   double aminor, rmaj, rwid, rmid, rc;
//   double rhots;
//   double rhotN;
//   double Cp2g;
//   double Cp2g2d1f;
//   double Cp2g1d;
//   std::vector<double> partrad, parttheta, partphitor, partvpar, partmu;
//   std::vector<double> partw, partfog, partg0;
//   int ntrack = 10;
//   std::vector<int> itrack;
//   double vts, paux_T_transit;
//   int ifilter = 0, filter_m0 = -5, filter_m1 = -2, filter_nc = 2;
//   std::vector<int> filter_m1d;
//   double pert_rc = 0.5;
//   double pert_rwid = 0.125, pert_thetac, pert_thetawid;
//   double pertamp = 4e-12, pert_lr = 0.0;
//   int pert_scheme = 3, pert_m0 = -5, pert_m1 = -2, pert_nc = 2, pert_nmin,
//       pert_nmax, pert_nstride, pert_lenm;
//   std::vector<int> pert_m1d;
//   double ant_rc, ant_rwid, antamp, ant_w;
//   int ant_scheme = 3, ant_m0 = -5, ant_m1 = -2, ant_nc = 2, ant_lenm = 4;
//   std::vector<int> ant_m1d;
//   int neocls;
//   int icol;
//   double nu_colN;
//   int v_par0 = 1, v_d = 0, v_mirror = 0, v_ExB = 0, v_Epar = 1, idwdt = 1;
//   int pullbackN;
//   int irec_track = 1, irec_Eparticle = 1;
//   int irec_converge;
//   bool iuseprv;
//   int isrcsnk;
//   std::vector<double> radsrcsnkL, radsrcsnkR;
//   std::vector<double> coefsrcsnkL, coefsrcsnkR;
//   double partwcut;
//   // Assuming timer_cls is a class or struct defined elsewhere
//   // timer_cls timer;
// public:
//   // Add a new species to the group
//   void addSpecies(size_t nparticles, double mass, double charge) {
//     group.addSpecies(nparticles, mass, charge);
//   }

//   // Push particles by modifying their coordinates
//   void pushParticle(size_t speciesIndex, double dt,
//                     const ParticleCoords &velocityChange) {
//     ParticleSpecies &species = group.getSpecies(speciesIndex);
//     ParticleCoords &coords = species.getCoords();

//     coords.axpy(velocityChange, dt); // Update position based on velocity change
//   }

//   void print() const { group.print(); }

//   Particle(Equilibrium &equ, std::optional<int> sp_id = std::nullopt) {
//     // Set species ID
//     this->sp_id = sp_id.value_or(1);

//     // Convert species ID to string - should be consistent with its size!!!
//     this->species_name = std::to_string(this->sp_id);

//     if (rank == 0) {
//       std::cout << "----initializing " << species_name << std::endl;
//     }

//     // 1. READINPUT
//     this->particle_cls_link2eq(equ);
//     this->particle_cls_set_parms(equ);

//     // PT Coordinates
//     this->partrad.resize(this->nptot);
//     this->parttheta.resize(this->nptot);
//     this->partphitor.resize(this->nptot);
//     this->partvpar.resize(this->nptot);
//     this->partmu.resize(this->nptot);
//     this->partw.resize(this->nptot);
//     this->partfog.resize(this->nptot);

//     if (this->ideltaf == 3) {
//       this->partg0.resize(this->nptot);
//     }

//     if (!this->irestart) {
//       // =================Load Markers=====================
//       this->particle_cls_loadmarker(equ);
//       // =================LOAD Perturbation=====================
//       this->particle_cls_loadpert(equ);
//     } else {
//       this->particle_cls_restart(0);
//     }

//     this->particle_cls_track_init(0);
//     this->particle_cls_showinfo();

//     if (rank == 0) {
//       std::cout << std::setw(20) << "Particle Class size:" << std::setw(10)
//                 << sizeof(*this) << " Byte " << std::endl;
//     }
//   }
//   void particle_cls_set_parms(Equilibrium &equ) {
//     // Default values
//     ischeme_motion = 2;
//     ngyro = 1;
//     irandom_gy = 1;
//     iuseprv = false;

//     isrcsnk = 0;
//     radsrcsnkL = {0.0, 0.1};
//     radsrcsnkR = {0.9, 1.0};
//     coefsrcsnkL = {0.0, 0.0, 0.0};
//     coefsrcsnkR = {0.0, 0.0, 0.0};

//     partwcut = 0.1;

//     pullbackN = 1;
//     neocls = 0;

//     // Set Run Parameters
//     ifile = 101;
//     std::ifstream inputFile("input");
//     if (!inputFile.is_open()) {
//       std::cerr << "Error opening file 'input'" << std::endl;
//       return;
//     }

//     for (fic = 1; fic <= sp_id; ++fic) {
//       // Read from file (assuming a similar structure to Fortran's namelist)
//       // This part needs to be adapted based on the actual file format
//     }
//     inputFile.close();

//     if (rank == 0) {
//       // Output particle settings
//       std::cout << "Species " << sp_id << ", nptot=" << nptot << std::endl;
//     }

//     // Assign values to class members
//     this->irestart = irestart;
//     this->nptot = nptot_all / mpisize;
//     this->nptot_all = this->nptot * mpisize;
//     if (rank == 0)
//       std::cout << "Species " << sp_id << ", nptot=" << this->nptot
//                 << std::endl;

//     this->ischeme_motion = ischeme_motion;
//     this->ngyro = ngyro;
//     this->irandom_gy = irandom_gy;
//     this->imixvar = imixvar;
//     this->ideltaf = ideltaf;
//     // Physics
//     this->nsonN = nsonN;
//     this->Tem = Tem;
//     this->mass = mass;
//     this->zcharge = zcharge;

//     this->vparmaxovt = vparmaxovt;
//     // Perturbation
//     this->pert_rc = pert_rc;
//     this->pert_rwid = pert_rwid;
//     this->pert_thetac = pert_thetac;
//     this->pert_thetawid = pert_thetawid;
//     this->pertamp = pertamp;
//     this->pert_lr = pert_lr;
//     this->pert_m0 = pert_m0;
//     this->pert_m1 = pert_m1;
//     this->pert_nc = pert_nc;
//     this->pert_nmin = pert_nmin;
//     this->pert_nmax = pert_nmax;
//     this->pert_nstride = pert_nstride;
//     this->pert_scheme = pert_scheme;
//     this->pert_lenm = std::abs(pert_m0 - pert_m1) + 1;
//     this->pert_m1d.resize(this->pert_lenm);
//     for (int fic = 1; fic <= this->pert_lenm; ++fic) {
//       this->pert_m1d[fic - 1] = std::min(pert_m0, pert_m1) - 1 + fic;
//     }
//     if (rank == 0) {
//       std::cout << "lenm=" << this->pert_lenm << ", pert_m1d=";
//       for (const auto &val : this->pert_m1d) {
//         std::cout << val << " ";
//       }
//       std::cout << std::endl;
//     }
//     // Antenna
//     this->ant_rc = ant_rc;
//     this->ant_rwid = ant_rwid;
//     this->antamp = antamp;
//     this->ant_m0 = ant_m0;
//     this->ant_m1 = ant_m1;
//     this->ant_nc = ant_nc;
//     this->ant_scheme = ant_scheme;

//     this->ant_lenm = std::abs(ant_m0 - ant_m1) + 1;
//     this->ant_m1d.resize(this->ant_lenm);

//     this->idensprof = idensprof;
//     this->iTemprof = iTemprof;
//     this->load_can = load_can;
//     this->ischeme_load = ischeme_load;
//     this->dens_coef1d = dens_coef1d;
//     this->Tem_coef1d = Tem_coef1d;

//     this->rmin = rmin;
//     this->rmax = rmax;
//     this->rwid = this->rmax - this->rmin;
//     this->rmid = (this->rmin + this->rmax) / 2.0;
//     this->rc = this->rmid;

//     this->ifilter = ifilter;
//     this->filter_nc = filter_nc;
//     this->filter_m0 = filter_m0;
//     this->filter_m1 = filter_m1;
//     nlen = std::abs(this->filter_m1 - this->filter_m0) + 1;
//     if (rank == 0)
//       std::cout << "nlen=" << nlen << std::endl;

//     this->filter_m1d.resize(nlen);
//     for (int fic = 1; fic <= nlen; ++fic) {
//       this->filter_m1d[fic - 1] = std::min(filter_m0, filter_m1) + (fic - 1);
//     }
//     if (rank == 0) {
//       std::cout << "filter_m1d=";
//       for (const auto &val : this->filter_m1d) {
//         std::cout << val << " ";
//       }
//       std::cout << std::endl;
//     }

//     this->Nphimult = Nphimult;
//     // Control variables
//     this->v_par0 = v_par0;
//     this->v_d = v_d;
//     this->v_mirror = v_mirror;
//     this->v_ExB = v_ExB;
//     this->v_Epar = v_Epar;
//     this->idwdt = idwdt;

//     this->irec_track = irec_track;
//     this->irec_Eparticle = irec_Eparticle;
//     this->irec_converge = irec_converge;
//     this->iuseprv = iuseprv;

//     this->ibc_particle = ibc_particle;
//     this->iset_zerosumw = iset_zerosumw;

//     this->icol = icol;
//     this->nu_colN = nu_colN;

//     this->isrcsnk = isrcsnk;
//     this->radsrcsnkL = radsrcsnkL;
//     this->radsrcsnkR = radsrcsnkR;
//     this->coefsrcsnkL = coefsrcsnkL;
//     this->coefsrcsnkR = coefsrcsnkR;
//     this->pullbackN = pullbackN;
//     this->neocls = neocls;
//     this->partwcut = partwcut;
//     this->profilename = profilename;

//     // Derived variables
//     nu_colstar = nu_colN / std::pow(std::sqrt(this->rc / this->rmaj), 3) *
//                  std::sqrt(2.0) * equ.getQlocRt(this->rc, 0.0) * this->rmaj *
//                  std::sqrt(this->mass / this->Tem);
//     // Set experimental n,T profiles
//     this->particleClsReadProfile1d();
//     // Derived parameters
//     this->phitormin = 0.0;
//     this->phitormax = this->phitormin + 2 * M_PI / this->Nphimult;
//     this->phitorwid = this->phitormax - this->phitormin;

//     this->vtovN =
//         std::sqrt(this->Tem / this->mass); // Calculate thermal velocity
//     this->vparmaxovN = this->vparmaxovt * this->vtovN;
//     this->mumaxovN2 = std::pow(this->vparmaxovN, 2) / (2 * this->Bax);
//     this->vts = std::sqrt(this->Tem / this->mass);

//     this->particleClsMeanDens(equ, this->rmin, this->rmax, 100, 200, this->Vtot,
//                               this->dens_intg, this->dens_mean, this->Stot);
//     this->Cp2g =
//         this->Vtot /
//         this->nptot_all; // 3d spline 2021/12/08; use NPTOT_ALL 2022/03/10!
//     this->Cp2g = this->Cp2g * this->nsonN * (-this->zcharge);
//     this->Cp2g = this->Cp2g * this->dens_mean;

//     this->Cp2g2d1f = this->Cp2g; // /this->phitorwid
//     if (rank == 0)
//       std::cout << "----set Cp2g2d1f=" << this->Cp2g2d1f
//                 << ",Cp2g=" << this->Cp2g << ",phitorwid=" << this->phitorwid
//                 << std::endl;
//     this->Cp2g1d = this->Vtot / this->nptot_all; // 1d spline 2022/08/18
//     this->Cp2g1d = this->Cp2g1d * this->nsonN * (-this->zcharge);
//     this->Cp2g1d = this->Cp2g1d * this->dens_mean;
//     this->Cp2g1d =
//         this->Cp2g1d / (4 * M_PI * M_PI * this->rmaj / this->Nphimult);

//     if (rank == 0)
//       std::cout << "  particle mean density=" << this->dens_mean << std::endl;

//     // vA^2=Cpoisson/Campere=1/betaN
//     wTAE = std::abs(0.5 / this->rmaj / equ.getQlocRt(this->rc, 0.0) *
//                     std::sqrt(1 / equ.betaN)); // from field_class
//     this->ant_w = ant_wowTAE * wTAE;
//     if (rank == 0)
//       std::cout << "ant_w=" << this->ant_w << ",wTAE=" << wTAE << std::endl;
//     for (int fic = 1; fic <= this->ant_lenm; ++fic) {
//       this->ant_m1d[fic - 1] = std::min(ant_m0, ant_m1) - 1 + fic;
//     }
//     if (rank == 0) {
//       std::cout << "lenm=" << this->ant_lenm << ", ant_m1d=";
//       for (const auto &val : this->pert_m1d) {
//         std::cout << val << " ";
//       }
//       std::cout << std::endl;
//     }

//     // (Rhots only as diagnosis)
//     this->rhots =
//         this->rhotN * std::sqrt(this->Tem * this->mass) / this->zcharge;
//     this->paux_T_transit = std::abs(2.0 * M_PI * equ.getQlocRt(this->rc, 0.0) *
//                                     equ.rmaxis / this->vts);

//     if (rank == 0) {
//       sform_r = "(A14,e12.5,A10,e12.5,A10,e12.5,A10,e12.5,A10,e12.5,A10,e12.5)";
//       sform_i = "(A14,I10,A14,I10,A14,I10,A14,I10,A14,I10,A14,I10)";

//       std::cout << "imixvar=" << this->imixvar
//                 << ",ibc_particle=" << this->ibc_particle
//                 << ",sp_id=" << this->sp_id << ","
//                 << "ideltaf=" << this->ideltaf << ",Nphimult=" << this->Nphimult
//                 << ",nptot=" << this->nptot << ","
//                 << "nptot_all=" << this->nptot_all
//                 << ",idensprof=" << this->idensprof
//                 << ", iTemprof=" << this->iTemprof << ","
//                 << "ntrack=" << this->ntrack << ",ifilter=" << this->ifilter
//                 << ",pert_scheme=" << this->pert_scheme << ","
//                 << "v_par0=" << this->v_par0 << ",v_d=" << this->v_d
//                 << ",v_mirror=" << this->v_mirror << ","
//                 << "irec_track=" << this->irec_track
//                 << ",irec_Eparticle=" << this->irec_Eparticle
//                 << ",irec_converge=" << this->irec_converge << ","
//                 << " filter_nc=" << this->filter_nc
//                 << " ischeme_motion=" << this->ischeme_motion
//                 << " ngyro=" << this->ngyro << " irandom_gy" << this->irandom_gy
//                 << ","
//                 << " ibc_particle=" << this->ibc_particle
//                 << " load_can=" << this->load_can
//                 << " ischeme_load=" << this->ischeme_load << ","
//                 << " pert_nmin=" << this->pert_nmin
//                 << " pert_nmax=" << this->pert_nmax
//                 << " pert_nstride=" << this->pert_nstride << ","
//                 << " pullbackN=" << this->pullbackN
//                 << " neocls=" << this->neocls << ", icol=" << this->icol
//                 << ", iset_zerosumw=" << this->iset_zerosumw << std::endl;

//       std::cout
//           << "phitormin=" << this->phitormin << ",phitormax=" << this->phitormax
//           << ",phitorwid=" << this->phitorwid << ","
//           << "nsonN=" << this->nsonN << ",Tem=" << this->Tem
//           << ",mass=" << this->mass << ","
//           << "zcharge=" << this->zcharge << ",vparmaxovt=" << this->vparmaxovt
//           << ",vparmaxovN=" << this->vparmaxovN << ","
//           << "vtovN=" << this->vtovN << ",c1=" << this->dens_coef1d[0]
//           << ",c2=" << this->dens_coef1d[1] << ","
//           << "c3=" << this->dens_coef1d[2] << ",c4=" << this->dens_coef1d[3]
//           << ",c1T=" << this->Tem_coef1d[0] << ","
//           << "c2T=" << this->Tem_coef1d[1] << ",c3T=" << this->Tem_coef1d[2]
//           << ",c4T=" << this->Tem_coef1d[3] << ","
//           << "Bax=" << this->Bax << ","
//           << "rmin=" << this->rmin << ",rmax=" << this->rmax
//           << ",Stot=" << this->Stot << ","
//           << "Vtot=" << this->Vtot << ",aminor=" << this->aminor
//           << ",rmaj=" << this->rmaj << ","
//           << "rwid=" << this->rwid << ",rmid=" << this->rmid
//           << ",rc=" << this->rc << ","
//           << "rhots=" << this->rhots << ",rhotN=" << this->rhotN
//           << ",Cp2g=" << this->Cp2g << ","
//           << "vts=" << this->vts << "paux_T_tr=" << this->paux_T_transit << ","
//           << "pert_rc=" << this->pert_rc << ",pert_rwid=" << this->pert_rwid
//           << ",pertamp=" << this->pertamp << ","
//           << "dens_mean=" << this->dens_mean << ",Cp2g1d=" << this->Cp2g1d
//           << ",nu_colN=" << this->nu_colN << ","
//           << "nu_colstar=" << nu_colstar << ",pert_thetac=" << this->pert_thetac
//           << ",pert_thetawid=" << this->pert_thetawid << ","
//           << "pert_lr=" << this->pert_lr << ",partwcut=" << this->partwcut
//           << std::endl;

//       std::cout << "isrcsnk=" << this->isrcsnk << std::endl;
//       std::cout << "radsrcsnkL=" << this->radsrcsnkL[0] << ","
//                 << this->radsrcsnkL[1]
//                 << ",coefsrcsnkL=" << this->coefsrcsnkL[0] << ","
//                 << this->coefsrcsnkL[1] << "," << this->coefsrcsnkL[2]
//                 << std::endl;
//       std::cout << "radsrcsnkR=" << this->radsrcsnkR[0] << ","
//                 << this->radsrcsnkR[1]
//                 << ",coefsrcsnkR=" << this->coefsrcsnkR[0] << ","
//                 << this->coefsrcsnkR[1] << "," << this->coefsrcsnkR[2]
//                 << std::endl;
//       std::cout << "iuseprv=" << this->iuseprv << std::endl;
//     }
//   }

//   void readInput(const std::string &inputFile) {
//     INIReader reader(inputFile);

//     // Check if reading was successful
//     if (reader.ParseError() < 0) {
//       throw std::runtime_error("Error: Unable to open or parse input file '" +
//                                inputFile + "'");
//     }

//     // Read key-value pairs
//     iequmodel = reader.GetInteger("Particle", "iequmodel", 2);
//     c1adhoc = reader.GetReal("Particle", "c1adhoc", 0.0);
//     c2adhoc = reader.GetReal("Particle", "c2adhoc", 0.0);
//     rmaxis_adhoc = reader.GetReal("Particle", "rmaxis_adhoc", 0.0);
//     Bmaxis_adhoc = reader.GetReal("Particle", "Bmaxis_adhoc", 0.0);
//     set_zerof = reader.GetBoolean("Particle", "set_zerof", false);
//     fname = reader.Get("Particle", "fname", "");
//     profilename = reader.Get("Particle", "profilename", "");

//     // Optional: Add validation and other logic
//     validateInputs();
//   }

//   void particle_cls_link2eq(Particle &thisParticle, const Equilibrium &equ) {
//     thisParticle.Bax = std::abs(equ.Bmaxis); // note sign!
//     thisParticle.aminor = equ.rdim / 2;
//     thisParticle.rmaj = equ.rmaxis;
//     thisParticle.rhotN = equ.rhoN;
//   }
// };

// class Particle {
// private:
//   double mass;
//   double zcharge;
//   std::vector<double> position;
//   std::vector<double> velocity;

// public:
//   int rank; // MPI rank of the process
//   int size; // Total number of MPI processes
//   Particle(double m, double q, std::vector<double> pos, std::vector<double>
//   vel)
//       : rank(-1), size(0){};

//   void move(const Field &field, double dt);
//   void displayState() const;

//   void initializeMPI() {
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
//   }
// };
#endif // PARTICLE_H
