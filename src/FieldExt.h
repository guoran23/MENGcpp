#ifndef FIELDEXTCLS_H
#define FIELDEXTCLS_H

#include "Field.h"  // Base class
#include "Equilibrium.h"
#include "Particle.h"
#include "MPIManager.h"
#include <vector>
#include <complex>

class FieldExtCls : public FieldCls {
    private:
        int rank, size;
    public:
        // Fields & moments arrays
        int ntotfem2d1f, ntotdof2d1f;
        std::vector<int> idxdof2d1f;
    
        std::vector<std::complex<double>> denskMkj, phik, jparkMkj, apark, aparsk, aparhk;
    
        // Solvers
        // SolverExtCls svPoisson, svAmpere, svAhcorr, svOhm;
        // DataC16ConvergeCls cvg;
    
        // Constructor & Destructor
        FieldExtCls() {
            std::cout << "FieldExtCls default constructor called"
            << std::endl;
        };
        ~FieldExtCls() = default;
    
        // Methods
        void init(const Equilibrium& equ, std::vector<ParticleSpecies>& pt);
        void finalize();
        void restart(int flag){};
        void test();
    
        void g2p2g();
        void g2p2g_2sp();
    
        void solve_dAdt();
        void solve_poisson();
        void solve_ampere();
        void solve_poisson_ampere();
    
        void solve_Acorr();
        void pullback();
    
        void apply_vec_bc();
        void apply_buff();
    
        void record2h5();
        void record1d();
    };
    

// Implementation of init() method
void FieldExtCls::init(const Equilibrium& equ, std::vector<ParticleSpecies>& pt) {
    
    // Initialize MPI rank and shared memory parameters
    MPIManager &mpiManager = MPIManager::getInstance();
    rank = mpiManager.getRank();
    size = mpiManager.getSize();

    // Initialize parent class
    FieldCls::field_cls_init(equ);

    std::cout << "----init FieldExtCls----" << std::endl;

    this->ntotfem2d1f = this->lenntor * this->spc.get_ntot12fem();
    this->ntotdof2d1f = this->lenntor * this->spc.get_ntot12dof();

    // Allocate index array
    this->idxdof2d1f.resize(this->ntotdof2d1f);
    this->spc.spc_cls_calc_idxdof2d1f(this->lenntor, this->idxdof2d1f);

    // Allocate field and moment arrays
    this->denskMkj.resize(this->ntotfem2d1f);
    this->phik.resize(this->ntotfem2d1f);
    this->jparkMkj.resize(this->ntotfem2d1f);
    this->apark.resize(this->ntotfem2d1f);

    if (this->imixvar == 1) {
        this->aparsk.resize(this->ntotfem2d1f, 0.0);
        this->aparhk.resize(this->ntotfem2d1f, 0.0);

        if (this->irestart) {
            restart(0);
        }
    }

    // Initialize solvers
    if (rank == 0) {
        std::cout << "----init Poisson, Ampere, Ohm solver in FieldExtCls----" << std::endl;
    }

    // this->svPoisson.init(this->ntotdof2d1f, false, 0);

    // // Poisson solver setup
    // if (this->isolver == 1) {
    //     std::vector<int> ieq1dL_poi = {-14, -14, -14, -14, -14, 17, 18, 18, 19};
    //     std::vector<int> igij1dL_poi = {1, 2, 4, 5, 9, 0, 0, 0, 0};
    //     std::vector<std::vector<int>> idiff2dL_poi = {
    //         {1, 1, 0, 0, 0, 0}, {1, 0, 0, 1, 0, 0}, {0, 1, 1, 0, 0, 0},
    //         {0, 0, 1, 1, 0, 0}, {0, 0, 1, 0, 0, 1}, {0, 0, 0, 1, 1, 0},
    //         {0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0}
    //     };

    //     this->svPoisson.assembly_general2d1f(this->spc, equ, 1, idiff2dL_poi, ieq1dL_poi, igij1dL_poi, this->bcrad, this->ntor1d, pt);
    // } else {
    //     if (rank == 0) {
    //         std::cerr << "****Error: wrong FD%ISOLVER value" << std::endl;
    //     }
    // }

    // // Ampere solver initialization
    // if (this->imixvar == 0) {
    //     this->svAmpere.init(this->ntotdof2d1f, false, 0);
    // } else {
    //     this->svAmpere.init(this->ntotdof2d1f, true, 2);
    // }

    // if (this->isolver == 1) {
    //     std::vector<int> ieq1dL_amp = {-1, -1, -1, -1, -1, 27, 28, 28, 29, 116};
    //     std::vector<int> igij1dL_amp = {1, 2, 4, 5, 9, 0, 0, 0, 0, 0};
    //     std::vector<std::vector<int>> idiff2dL_amp = {
    //         {1, 1, 0, 0, 0, 0}, {1, 0, 0, 1, 0, 0}, {0, 1, 1, 0, 0, 0},
    //         {0, 0, 1, 1, 0, 0}, {0, 0, 1, 0, 0, 1}, {0, 0, 0, 1, 1, 0},
    //         {0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0}
    //     };

    //     this->svAmpere.assembly_general2d1f(this->spc, equ, 1, idiff2dL_amp, ieq1dL_amp, igij1dL_amp, this->bcrad, this->ntor1d, pt);
    // }
}
#endif // FIELDEXTCLS_H