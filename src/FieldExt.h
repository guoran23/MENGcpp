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
        void initializePerturbations();
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

    if (this->imixvar == 1) {
        this->aparsk.resize(this->ntotfem2d1f, 0.0);
        this->aparhk.resize(this->ntotfem2d1f, 0.0);

        if (this->irestart) {
            restart(0);
        }
    }

    // Initialize perturbation fields==0
    for (int i = 0; i < this->ntotfem2d1f; ++i) {
        this->denskMkj[i] = 0.0;
        this->phik[i] = 0.0;
        this->jparkMkj[i] = 0.0;
        this->apark[i] = 0.0;
        if (this->imixvar == 1) {
            this->aparsk[i] = 0.0;
            this->aparhk[i] = 0.0;
        }
    }
    // Initialize Ana Gauss perturbation fields
    initializePerturbations();

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

void FieldExtCls::initializePerturbations() {
    double dtheta = 2 * M_PI / nthe;
    double drad = (radmax - radmin) / nrad;

     // 直接在函数中定义并初始化
    std::vector<double> rc1_arr(this->lenntor, 0.5);
    std::vector<double> amp_arr(this->lenntor, 0.1);
    std::vector<double> rwidth_arr(this->lenntor, 0.1);
    std::vector<int> mpoloidal_arr(this->lenntor);

    for (int i = 0; i < this->lenntor; ++i) {
        mpoloidal_arr[i] = 5 + i;  // 例如：每个 toroidal mode 不同的 poloidal 模数
    }

    for (int itor = 0; itor < this->lenntor; ++itor) {
        double rc1 = rc1_arr[itor];
        double amp = amp_arr[itor];
        double rwidth = rwidth_arr[itor];
        int mpoloidal = mpoloidal_arr[itor];

        for (int i = 0; i < nradfem; ++i) {
            for (int j = 0; j < nthefem; ++j) {
                int idx = itor * nradfem * nthefem + i * nthefem + j;

                double r = radmin + i * drad;
                double theta = themin + j * dtheta;

                // 高斯径向部分
                double radial_part = amp * std::exp( - std::pow((r - rc1) / rwidth, 2) );

                // 角向复数部分
                std::complex<double> angular_part = std::exp(std::complex<double>(0, mpoloidal * theta));

                this->phik[idx] = radial_part * angular_part;
            }
        }
    }
}
#endif // FIELDEXTCLS_H