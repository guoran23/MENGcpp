#ifndef PARTICLE_H
#define PARTICLE_H

#include "../inih/INIReader.h"
#include "Equilibrium.h"
#include "Field.h"
#include "util_math.h"
#include <cmath>
#include <functional> // std::plus, std::multiplies
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <random>
#include <string>
#include <vector>

class ParticleCoords {
private:
  int nptot;

public:
  std::vector<double> partrad, parttheta, partphitor;
  std::vector<double> partvpar;
  std::vector<double> partw;
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
  int nptot;
  long long nptot_all = 1000000; // particle # in all PEs
  std::string species_name;
  double mass = 1.0;
  double zcharge = -1.0;
  ParticleCoords coords;
  double Tem = 1.0;
  double nsonN = 1.0;
  int sp_id = 0; //start from 0

public:
  std::vector<double> partmu;
  std::vector<double> partfog, partg0;

  double rmin = 0.2, rmax = 0.8;
  // Constructor
  ParticleSpecies(size_t nparticles, double m, double q,
                  const std::string &name = "unknown")
      : nptot(nparticles), coords(nparticles), mass(m), zcharge(q),
        species_name(name) {

    // Initialize partmu, partfog, partg0
    partmu.resize(nptot);
    partfog.resize(nptot);
    partg0.resize(nptot);
  }

  // Accessors
  int getNptot() const { return nptot; }
  long long getNptotAll() const { return nptot_all; }
  double getMass() const { return mass; }
  double getCharge() const { return zcharge; }
  int getSpId() const { return sp_id; }
  ParticleCoords &getCoords() { return coords; }
  const ParticleCoords &getCoords() const { return coords; }
  double getTem() const { return Tem; }
  double getNsonN() const { return nsonN; }
  void setTem(double t) { Tem = t; }
  void setNsonN(double n) { nsonN = n; }
  void serSpId(int id) { sp_id = id; } 

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
  ParticleSpecies &getSpecies(size_t index) {
    if (index < speciesList.size()) {
      return *speciesList[index];
    }
    throw std::out_of_range(
        "Index out of bounds"); // Throw an exception if index is out of bounds
  }

  int getSpeciesCount() const { return nsp; }
  size_t getTotalParticles() const {
    size_t total = 0;
    for (const auto &species : speciesList) {
      total += species->getNptot();
    }
    return total;
  }

  // 打印所有粒子种类信息
  void print() const {
    std::cout << "Particle Group Info:\n";
    for (const auto &species : speciesList) {
      species->print();
      std::cout << "-----------------------------\n";
    }
  }
};

class Particle {
private:
  int nsp;
  ParticleGroup group;
  Equilibrium equ;

public:
  int ischeme_motion, ngyro, irandom_gy;
  int imixvar = 1;
  int ibc_particle = 0, iset_zerosumw;
  bool irestart = false;
  std::string species_name;
  int ideltaf = 2;
  double phitormin = 0.0, phitormax, phitorwid;
  int Nphimult = 2;  // 1/Nphimult torus
  // double nsonN = 1.0;
  // double Tem = 1.0;

  double vparmaxovt = 5.0;
  double vparmaxovN;
  double mumaxovN2;
  double vtovN;
  std::string profilename;
  std::vector<double> dens_coef1d = {0.49123, 0.298228, 0.198739, 0.521298};
  std::vector<double> Tem_coef1d = {0.49123, 0.298228, 0.198739, 0.521298};
  int idensprof;
  int iTemprof;
  // Assuming bspline_1d is a class or struct defined elsewhere
  // bspline_1d dens_bsp, Tem_bsp;
  int load_can, ischeme_load;
  double Bax;
  double rmin = 0.2, rmax = 0.8;
  double Stot, Vtot, dens_mean, dens_intg;
  double aminor, rmaj, rwid, rmid, rc;
  double rhots;
  double rhotN;
  double Cp2g;
  double Cp2g2d1f;
  double Cp2g1d;
  // std::vector<double> partrad, parttheta, partphitor, partvpar, partmu;
  // std::vector<double> partw, partfog, partg0;

  int ntrack = 10;
  std::vector<int> itrack;
  double vts, paux_T_transit;
  int ifilter = 0, filter_m0 = -5, filter_m1 = -2, filter_nc = 2;
  std::vector<int> filter_m1d;
  // !--perturbation--
  double pert_rc = 0.5;
  double pert_rwid = 0.125, pert_thetac, pert_thetawid;
  double pertamp = 4e-12, pert_lr = 0.0;
  int pert_scheme = 3, pert_m0 = -5, pert_m1 = -2, pert_nc = 2, pert_nmin,
      pert_nmax, pert_nstride, pert_lenm;
  std::vector<int> pert_m1d;
  // !--perturbation--
  double ant_rc, ant_rwid, antamp, ant_w;
  int ant_scheme = 3, ant_m0 = -5, ant_m1 = -2, ant_nc = 2, ant_lenm = 4;
  std::vector<int> ant_m1d;
  // !--neoclassical term--
  int neocls;
  // !--collision--
  int icol = 0;
  double nu_colN = 0.0;
  int v_par0 = 1, v_d = 0, v_mirror = 0, v_ExB = 0, v_Epar = 1, idwdt = 1;
  int pullbackN;

  int irec_track = 1, irec_Eparticle = 1;
  int irec_converge;
  bool iuseprv;
  int isrcsnk;
  std::vector<double> radsrcsnkL, radsrcsnkR;
  std::vector<double> coefsrcsnkL, coefsrcsnkR;
  double partwcut;
  //
  // timer_cls timer;
public:
  void readInput(const std::string &filepath, int mpisize) {
    INIReader reader(filepath);
    std::cout << "Reading input file: " << filepath << std::endl;
    std::cout << "Number of MPI processes: " << mpisize << std::endl;
    if (reader.ParseError() < 0) {
      std::cerr << "Error: Unable to load INI file: " << filepath << std::endl;
    }
    nsp = reader.GetInteger("MENG", "nsp", 0);
    // Iterate over predefined sections
    for (int i = 1; i <= nsp; ++i) {
      std::string section = "Particle_" + std::to_string(i);
      std::cout << "Processing section: " << section << std::endl;

      size_t nptot_all = reader.GetInteger(section, "nptot_all", 0);
      double mass = reader.GetReal(section, "mass", 0.0);
      double zcharge = reader.GetReal(section, "zcharge", 0.0);
      std::string species_name = reader.Get(section, "species_name", "unknown");
      double rmin = reader.GetReal(section, "rmin", 0.0);
      double rmax = reader.GetReal(section, "rmax", 0.0);
      double Tem = reader.GetReal(section, "Tem", 1.0);
      double nsonN = reader.GetReal(section, "nsonN", 1.0);

      if (nptot_all > 0) {
        int nptot = nptot_all / mpisize; // calculate particles per process
        nptot_all = nptot * mpisize;     // update total number of particles
        auto species = std::make_shared<ParticleSpecies>(nptot, mass, zcharge,
                                                         species_name);
        group.addSpecies(species);
        species->rmin = rmin;
        species->rmax = rmax;
        species->setTem(Tem);
        species->setNsonN(nsonN);
        species->serSpId(i-1); 
        std::cout << "  Species '" << species_name << "' added successfully.\n";
        // Log the parsed values
        std::cout << "  nptot: " << nptot << ", mass: " << mass
                  << ", charge: " << zcharge << std::endl;
      }
    }
  }

  // Push particles by modifying their coordinates
  void pushParticle(size_t speciesIndex, double dt,
                    const ParticleCoords &velocityChange) {
    ParticleSpecies &species = group.getSpecies(speciesIndex);
    ParticleCoords &coords = species.getCoords();

    coords.axpy(velocityChange, dt); // Update position based on velocity change
  }

  void print() const { group.print(); }

  Particle(Equilibrium &equ_in, int rank, int mpisize_in) : equ(equ_in) {

    std::cout << "Instance Particle @process " << rank << " of " << mpisize_in
              << ".\n";
    // 1. READINPUT
    readInput("input.ini", mpisize_in);
    if (rank == 0) {
      std::cout << "Finish particle readInput...\n"
                << "nsp =" << nsp << std::endl;
    }
    // 2. Summay of particles
    std::cout << "On process: " << rank << ", Total particles of all species: "
              << group.getTotalParticles() << std::endl;
    // 3.
    this->particle_cls_link2eq(equ);
    // 4. =================Load Markers=====================
    for (int i = 0; i < nsp; i++) {
      {
        auto &species = group.getSpecies(i);
        this->particle_cls_set_parms(equ, species, rank);
        std::cout << "Finish particle_cls_set_parms...\n";
        this->particle_cls_loadmarker(equ, species);
      }
    }
  }

  void particle_cls_link2eq(const Equilibrium &equ) {
    Bax = std::abs(equ.Bmaxis); // note sign!
    aminor = equ.rdim / 2;
    rmaj = equ.rmaxis;
    rhotN = equ.rhoN;
  }

  void random_seed() {
    // Initialize random seed
    srand(time(0));
  }

  void random_number(std::vector<double> &vec) {
    for (auto &v : vec) {
      v = static_cast<double>(rand()) / RAND_MAX;
    }
  }

  double getTem1d(double radius) {
    // Implement function to get temperature based on radius
    return 1.0; // Placeholder
  }

  double getdens1d(double radius) {
    // Implement function to get density
    return 1.0; // Placeholder
  }

  double getdens1d_can(Equilibrium &equ, double radius, double theta,
                       double vpar, double mu) {
    // Implement function to get canonical density
    return 1.0; // Placeholder
  }

  double gennor(double mu, double sigma) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(mu, sigma);
    return dis(gen);
  }

  double particle_cls_getf0(double partrad, double parttheta, double partvpar,
                            double partmu, Equilibrium &equ) {
    // Implement function to get f0
    return 1.0; // Placeholder
  }

  

  void particle_cls_mean_dens(const Equilibrium &equ, double rmin, double rmax,
                              int nrad, int nthe, double &volume,
                              double &int_dens, double &mean_dens,
                              double &Stot) {
    // Working variables
    std::vector<double> xrad(nrad), wrad(nrad);
    std::vector<double> xthe(nthe), wthe(nthe);

    // Initialize outputs
    mean_dens = 0.0;
    int_dens = 0.0;
    volume = 0.0;
    Stot = 0.0;

    // Compute quadrature points and weights
    std::cout << "Computing quadrature points and weights...\n";
    std::cout << "rmin=" << rmin << ", rmax=" << rmax << std::endl;
    std::cout << "nrad=" << nrad << ", nthe=" << nthe << std::endl;
    std::cout << "xrad.size()=" << xrad.size()
              << ", xthe.size()=" << xthe.size() << std::endl;
    std::cout << "wrad.size()=" << wrad.size()
              << ", wthe.size()=" << wthe.size() << std::endl;
    std::cout << " this->phitorwid=" << this->phitorwid << std::endl;
    UtilMath::lgwt(nrad, rmin, rmax, xrad, wrad);
    double twopi = 2.0 * M_PI;
    UtilMath::lgwt(nthe, 0.0, twopi, xthe, wthe);

    // Debug output (assuming rank == 0)
    std::cout << "sum(wrad) = "
              << std::accumulate(wrad.begin(), wrad.end(), 0.0)
              << ", sum(wthe) = "
              << std::accumulate(wthe.begin(), wthe.end(), 0.0) << std::endl;
    std::cout << "mean(xrad) = "
              << std::accumulate(xrad.begin(), xrad.end(), 0.0) / nrad
              << ", mean(xthe) = "
              << std::accumulate(xthe.begin(), xthe.end(), 0.0) / nthe
              << std::endl;

    // Nested loop for computation
    for (int fic = 0; fic < nrad; ++fic) {
      for (int fjc = 0; fjc < nthe; ++fjc) {
        double w12 = wrad[fic] * wthe[fjc];
        double jaco2 = equ.getjaco2(xrad[fic], xthe[fjc]);
        double jaco3 = jaco2 * equ.getR(xrad[fic], xthe[fjc]);

        volume += w12 * jaco3;
        int_dens += w12 * jaco3 * this->getdens1d(xrad[fic]);
        Stot += w12 * jaco2;
      }
    }

    // Final adjustments
    volume = std::abs(volume) * this->phitorwid;
    int_dens = std::abs(int_dens) * this->phitorwid;
    Stot = std::abs(Stot);

    // Compute mean density
    mean_dens = int_dens / volume;
  }

  void particle_cls_set_parms(Equilibrium &equ, ParticleSpecies &species,
                              int rank) {

    // Initialize default values
    ischeme_motion = 2;
    ngyro = 1;
    irandom_gy = 1;
    iuseprv = false;

    isrcsnk = 0;
    radsrcsnkL = {0.0, 0.1};
    radsrcsnkR = {0.9, 1.0};
    coefsrcsnkL = {0.0, 0.0, 0.0};
    coefsrcsnkR = {0.0, 0.0, 0.0};

    partwcut = 1e-1;

    pullbackN = 1;
    neocls = 0;

    // Set Run Parameters
    // int ifile = 101;
    // std::ifstream infile("input");
    // if (!infile.is_open()) {
    //   std::cerr << "Failed to open input file." << std::endl;
    //   return;
    // }
    // for (int fic = 1; fic <= sp_id; ++fic) {
    //   // Read namelist variables from file
    //   // This requires parsing logic based on the file structure
    //   // Placeholder for reading logic
    // }
    // infile.close();

    // Todo: Read from input file
    //       parameters for each species and parameters for all species

    // Marker
    int nptot = species.getNptot();
    int nptot_all = species.getNptotAll();
    double mass = species.getMass();
    double zcharge = species.getCharge();
    double nsonN = species.getNsonN();
    double rmin = species.rmin;
    double rmax = species.rmax;
    double Tem = species.getTem();
    int sp_id = species.getSpId();

    // Set parameters for all species
    ischeme_motion = ischeme_motion;
    ngyro = ngyro;
    irandom_gy = irandom_gy;
    imixvar = imixvar;
    ideltaf = ideltaf;
    // Physics
    vparmaxovt = vparmaxovt;
    // Perturbation
    pert_rc = pert_rc;
    pert_rwid = pert_rwid;
    pert_thetac = pert_thetac;
    pert_thetawid = pert_thetawid;
    pertamp = pertamp;
    pert_lr = pert_lr;
    pert_m0 = pert_m0;
    pert_m1 = pert_m1;
    pert_nc = pert_nc;
    pert_nmin = pert_nmin;
    pert_nmax = pert_nmax;
    pert_nstride = pert_nstride;
    pert_scheme = pert_scheme;
    pert_lenm = std::abs(pert_m0 - pert_m1) + 1;
    pert_m1d.resize(pert_lenm);
    int min_pert_m = std::min(pert_m0, pert_m1);
    for (int fic = 0; fic < pert_lenm; ++fic) {
      pert_m1d[fic] = min_pert_m - 1 + (fic + 1);
    }
    if (rank == 0) {
      std::cout << "lenm=" << pert_lenm << ", pert_m1d=";
      for (const auto &val : pert_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    // Antenna
    ant_rc = 0.5;
    ant_rwid = 0.125;
    antamp = 1e-5;
    double ant_wowTAE = 1.0;
    ant_scheme = 3;
    ant_m0 = -11;
    ant_m1 = -10;
    ant_nc = 6;

    ant_lenm = std::abs(ant_m0 - ant_m1) + 1;
    ant_m1d.resize(ant_lenm);
    int min_ant_m = std::min(ant_m0, ant_m1);
    for (int fic = 0; fic < ant_lenm; ++fic) {
      ant_m1d[fic] = min_ant_m - 1 + (fic + 1);
    }
    if (rank == 0) {
      std::cout << "lenm=" << ant_lenm << ", ant_m1d=";
      for (const auto &val : ant_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    // Density and Temperature Profiles
    idensprof = 1;
    iTemprof = 1;
    load_can = 0;
    ischeme_load = 1;
    // Assuming dens_coef1d and Tem_coef1d are already populated
    // If not, they should be read from inputs or initialized here
    rwid = rmax - rmin;
    rmid = (rmin + rmax) / 2.0;
    rc = rmid;

    ifilter = ifilter;
    filter_nc = filter_nc;
    filter_m0 = filter_m0;
    filter_m1 = filter_m1;
    int nlen_filter = std::abs(filter_m1 - filter_m0) + 1;
    if (rank == 0) {
      std::cout << "nlen=" << nlen_filter << std::endl;
    }

    std::vector<int> filter_m1d(nlen_filter);
    int min_filter_m = std::min(filter_m0, filter_m1);
    for (int fic = 0; fic < nlen_filter; ++fic) {
      filter_m1d[fic] = min_filter_m + fic;
    }
    if (rank == 0) {
      std::cout << " filter_m1d=";
      for (const auto &val : filter_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    Nphimult = Nphimult;
    // Control variables
    v_par0 = v_par0;
    v_d = v_d;
    v_mirror = v_mirror;
    v_ExB = v_ExB;
    v_Epar = v_Epar;
    idwdt = idwdt;

    irec_track = irec_track;
    irec_Eparticle = irec_Eparticle;
    irec_converge = irec_converge;
    iuseprv = iuseprv;

    ibc_particle = ibc_particle;
    iset_zerosumw = iset_zerosumw;

    isrcsnk = isrcsnk;
    radsrcsnkL = radsrcsnkL;
    radsrcsnkR = radsrcsnkR;
    coefsrcsnkL = coefsrcsnkL;
    coefsrcsnkR = coefsrcsnkR;
    pullbackN = pullbackN;
    neocls = neocls;
    partwcut = partwcut;
    profilename = profilename;

    // Derived variables
    double nu_colstar = nu_colN / std::pow(std::sqrt(rc / rmaj), 3) *
                        std::sqrt(2.0) * equ.getqloc_rt(rc, 0.0) * rmaj *
                        std::sqrt(mass / Tem);
    // Set experimental n,T profiles
    // particle_cls_read_profile1d();
    // Derived parameters
    phitormin = 0.0;
    phitormax = phitormin + 2.0 * M_PI / static_cast<double>(Nphimult);
    phitorwid = phitormax - phitormin;

    vtovN = std::sqrt(Tem / mass); // calculate thermal velocity
    vparmaxovN = vparmaxovt * vtovN;
    std::cout << "vparmaxovN=" << vparmaxovN << std::endl;
    mumaxovN2 = std::pow(vparmaxovN, 2) / (2.0 * Bax);
    std::cout << "mumaxovN2=" << mumaxovN2 << std::endl;
    vts = std::sqrt(Tem / mass);
    // --
    // // test lgwt
    // std::vector<double> xrad_test(10), wrad_test(10);
    // UtilMath::lgwt(10, 0.2, 0.8, xrad_test, wrad_test);
    // for (int i = 0; i < 10; i++) {
    //   std::cout << "xrad_test[" << i << "]=" << xrad_test[i] << ", wrad_test["
    //             << i << "]=" << wrad_test[i] << std::endl;
    // }

    this->particle_cls_mean_dens(equ, this->rmin, this->rmax, 100, 200,
                                 this->Vtot, this->dens_intg, this->dens_mean,
                                 this->Stot);
    std::cout << "Vtot=" << Vtot << ", dens_intg=" << dens_intg
              << ", dens_mean=" << dens_mean << ", Stot=" << Stot << std::endl;
    // this->Stot = M_PI * (std::pow(rmax, 2) - std::pow(rmin, 2));
    // this->Vtot = this->Stot * this->rmaj * this->phitorwid;
    Cp2g = Vtot /
           static_cast<double>(
               nptot_all); // 3d spline 2021/12/08; use NPTOT_ALL 2022/03/10!
    Cp2g *= nsonN * (-zcharge);
    Cp2g *= dens_mean;

    Cp2g2d1f = Cp2g; // / phitorwid;
    if (rank == 0) {
      std::cout << "----set Cp2g2d1f=" << Cp2g2d1f << ", Cp2g=" << Cp2g
                << ", phitorwid=" << phitorwid << std::endl;
    }
    Cp2g1d = Vtot / static_cast<double>(nptot_all); // 1d spline 2022/08/18
    Cp2g1d *= nsonN * (-zcharge);
    Cp2g1d *= dens_mean;
    Cp2g1d /= (4.0 * std::pow(M_PI, 2) * rmaj / static_cast<double>(Nphimult));

    if (rank == 0) {
      std::cout << "  particle mean density=" << dens_mean << std::endl;
    }

    // vA^2 = Cpoisson/Campere = 1/betaN
    double wTAE = std::abs(0.5 / rmaj / equ.getqloc_rt(rc, 0.0) *
                           std::sqrt(1.0 / equ.betaN)); // from field_class
    ant_w = ant_wowTAE * wTAE;
    if (rank == 0) {
      std::cout << "ant_w=" << ant_w << ", wTAE=" << wTAE << std::endl;
    }
    for (int fic = 0; fic < ant_lenm; ++fic) {
      ant_m1d[fic] = std::min(ant_m0, ant_m1) - 1 + (fic + 1);
    }
    if (rank == 0) {
      std::cout << "lenm=" << ant_lenm << ", ant_m1d=";
      for (const auto &val : ant_m1d) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }

    // (Rhots only as diagnosis)
    rhots = rhotN * std::sqrt(Tem * mass) / zcharge;
    paux_T_transit =
        std::abs(2.0 * M_PI * equ.getqloc_rt(rc, 0.0) * equ.rmaxis / vts);

    if (rank == 0) {
      std::string sform_r =
          "(A14,e12.5,A10,e12.5,A10,e12.5,A10,e12.5,A10,e12.5,A10,e12.5)";
      std::string sform_i = "(A14,I10,A14,I10,A14,I10,A14,I10,A14,I10,A14,I10)";

      // Since C++ doesn't support formatted write like Fortran, we use
      // std::cout with appropriate formatting
      std::cout << std::fixed << std::setprecision(5);
      std::cout << "imixvar=" << imixvar << ", ibc_particle=" << ibc_particle
                << ", sp_id=" << sp_id << ", ideltaf=" << ideltaf
                << ", Nphimult=" << Nphimult << ", nptot=" << nptot
                << ", nptot_all=" << nptot_all << ", idensprof=" << idensprof
                << ", iTemprof=" << iTemprof << ", ntrack=" << irec_track
                << ", ifilter=" << ifilter << ", pert_scheme=" << pert_scheme
                << ", v_par0=" << v_par0 << ", v_d=" << v_d
                << ", v_mirror=" << v_mirror << ", irec_track=" << irec_track
                << ", irec_Eparticle=" << irec_Eparticle
                << ", irec_converge=" << irec_converge
                << ", filter_nc=" << filter_nc
                << ", ischeme_motion=" << ischeme_motion << ", ngyro=" << ngyro
                << ", irandom_gy=" << irandom_gy
                << ", ibc_particle=" << ibc_particle
                << ", load_can=" << load_can
                << ", ischeme_load=" << ischeme_load
                << ", pert_nmin=" << pert_nmin << ", pert_nmax=" << pert_nmax
                << ", pert_nstride=" << pert_nstride
                << ", pullbackN=" << pullbackN << ", neocls=" << neocls
                << ", icol=" << icol << ", iset_zerosumw=" << iset_zerosumw
                << std::endl;

      std::cout << "phitormin=" << phitormin << ", phitormax=" << phitormax
                << ", phitorwid=" << phitorwid << ", nsonN=" << nsonN
                << ", Tem=" << Tem << ", mass=" << mass
                << ", zcharge=" << zcharge << ", vparmaxovt=" << vparmaxovt
                << ", vparmaxovN=" << vparmaxovN << ", vtovN=" << vtovN
                << ", c1=" << dens_coef1d[0] << ", c2=" << dens_coef1d[1]
                << ", c3=" << dens_coef1d[2] << ", c4=" << dens_coef1d[3]
                << ", c1T=" << Tem_coef1d[0] << ", c2T=" << Tem_coef1d[1]
                << ", c3T=" << Tem_coef1d[2] << ", c4T=" << Tem_coef1d[3]
                << ", Bax=" << Bax << ", rmin=" << rmin << ", rmax=" << rmax
                << ", Stot=" << Stot << ", Vtot=" << Vtot
                << ", aminor=" << 0.0 // Assuming aminor is defined elsewhere
                << ", rmaj=" << rmaj << ", rwid=" << rwid << ", rmid=" << rmid
                << ", rc=" << rc << ", rhots=" << rhots << ", rhotN=" << rhotN
                << ", Cp2g=" << Cp2g << ", vts=" << vts
                << ", paux_T_tr=" << paux_T_transit << ", pert_rc=" << pert_rc
                << ", pert_rwid=" << pert_rwid << ", pertamp=" << pertamp
                << ", dens_mean=" << dens_mean << ", Cp2g1d=" << Cp2g1d
                << ", nu_colN=" << nu_colN << ", nu_colstar=" << nu_colstar
                << ", pert_thetac=" << pert_thetac
                << ", pert_thetawid=" << pert_thetawid
                << ", pert_lr=" << pert_lr << ", partwcut=" << partwcut
                << std::endl;

      std::cout << "isrcsnk=" << isrcsnk << std::endl;
      std::cout << "radsrcsnkL=";
      for (const auto &val : radsrcsnkL) {
        std::cout << val << " ";
      }
      std::cout << ", coefsrcsnkL=";
      for (const auto &val : coefsrcsnkL) {
        std::cout << val << " ";
      }
      std::cout << std::endl;

      std::cout << "radsrcsnkR=";
      for (const auto &val : radsrcsnkR) {
        std::cout << val << " ";
      }
      std::cout << ", coefsrcsnkR=";
      for (const auto &val : coefsrcsnkR) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
      std::cout << "iuseprv=" << iuseprv << std::endl;
    }
  }

  void particle_cls_loadmarker(Equilibrium &equ, ParticleSpecies &species) {
    std::cout << "--------LOAD MARKER--------" << std::endl;

    auto &coords = species.getCoords();
    auto &partrad = coords.partrad;
    auto &parttheta = coords.parttheta;
    auto &partphitor = coords.partphitor;
    auto &partvpar = coords.partvpar;
    auto &partw = coords.partw;

    auto &partmu = species.partmu;
    auto &partfog = species.partfog;
    auto &partg0 = species.partg0;
    int nptot = species.getNptot();

    const double mass = species.getMass();

    // Init Random
    random_seed();
    random_number(partrad);
    std::cout << "partrad[0]: " << partrad[0] << std::endl;

    double c0rnd, c1rnd;
    double Tpar, Tperp, zwidprod, BB;
    int fic;

    if (ischeme_load == 1) {
      for (fic = 0; fic < nptot; ++fic) {
        c0rnd = rmin * rmin;
        c1rnd = rmax * rmax - rmin * rmin;
        partrad[fic] = sqrt(partrad[fic] * c1rnd + c0rnd);
      }
    } else if (ischeme_load == 2) {
      for (fic = 0; fic < nptot; ++fic) {
        partrad[fic] = partrad[fic] * rwid + rmin;
      }
    }

    // 2.3 Theta, Phitor
    random_number(parttheta);
    random_number(partphitor);
    for (fic = 0; fic < nptot; ++fic) {
      parttheta[fic] = parttheta[fic] * 2.0 * M_PI;
      partphitor[fic] = partphitor[fic] * phitorwid + phitormin;
    }
    std::cout << "parttheta[0]: " << parttheta[0] << std::endl;
    // 4. Vpar
    double vpar0 = 0.0;
    double rndmu4 = vpar0;
    double rndsd4;
    for (fic = 0; fic < nptot; ++fic) {
      Tpar = getTem1d(partrad[fic]);
      rndsd4 = sqrt(Tpar / mass) * sqrt(1.0 / 2.0);
      double rndf8 = gennor(rndmu4, rndsd4);
      while (abs(rndf8) >= vparmaxovN) {
        rndf8 = gennor(rndmu4, rndsd4);
      }
      partvpar[fic] = rndf8;
    }
    std::cout << "partvpar[0]: " << partvpar[0] << std::endl;
    // 5. Mu
    zwidprod = phitorwid * rwid * 2.0 * M_PI;

    for (fic = 0; fic < nptot; ++fic) {
      Tperp = getTem1d(partrad[fic]);
      double rndf8 = static_cast<double>(rand()) / RAND_MAX;

      BB = equ.getB(partrad[fic], parttheta[fic]);
      // std::cout << "partmu[" << fic << "]: " << partmu[fic] << std::endl;

      partmu[fic] = -Tperp / mass * std::log(rndf8) / BB / 2.0;

      while (partmu[fic] >= mumaxovN2) {
        rndf8 = static_cast<double>(rand()) / RAND_MAX;
        partmu[fic] = -Tperp / mass * std::log(rndf8) / BB / 2.0;
      }
      
      // 1b. full weight
      if (load_can == 0) {
        partfog[fic] = getdens1d(partrad[fic]) / dens_mean;
      } else if (load_can == 1) {
        partfog[fic] = getdens1d_can(equ, partrad[fic], parttheta[fic],
                                     partvpar[fic], partmu[fic]) /
                       dens_mean;
      }

      if (ischeme_load == 1) {
        partfog[fic] *= equ.getjaco3(partrad[fic], parttheta[fic]) * zwidprod /
                        Vtot * rmid / partrad[fic];
      } else if (ischeme_load == 2) {
        partfog[fic] *=
            equ.getjaco3(partrad[fic], parttheta[fic]) * zwidprod / Vtot;
      } else {
        std::cerr << "Error: wrong ISCHEME_LOAD value; set to 1" << std::endl;
      }

      // 1a. weight
      if (ideltaf > 0) {
        partw[fic] = 0.0;
      } else if (ideltaf == 0) {
        partw[fic] = partfog[fic];
      } else {
        std::cerr << "Error: wrong IDELTAF value; set to 1" << std::endl;
      }

      if (ideltaf == 3) {
        partg0[fic] = particle_cls_getf0(partrad[fic], parttheta[fic],
                                         partvpar[fic], partmu[fic], equ) /
                      partfog[fic];
      }

    std::cout << "partmu[" << fic << "]: " << partmu[fic] << std::endl;  
    }
  }
};

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
