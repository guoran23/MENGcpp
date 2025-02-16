#ifndef SPLINECUBICND_H
#define SPLINECUBICND_H

#include <mpi.h>
// #include <petscksp.h>
#include "MPIManager.h"
#include "util_math.h"
#include "util_io.h"
#include <optional>
#include <vector>

// Assuming util_math, mpi_class, solver_class are translated into C++ headers
// #include "util_math.h"
// #include "mpi_class.h"
// #include "solver_class.h"

const double TWOPI = 2.0 * M_PI;

class SplineCubicNdCls {
private:
  int ndim;      // Number of dimensions
  int norder;    // Spline order (e.g., 4 for cubic)
  int nintg;     // Integration points
  int nintg_tot; // Total integration points
  double xshift; // Coordinate shift

  std::vector<double> xknot;  // Knot vector
  std::vector<double> xx_tot; // Total coordinates
  std::vector<double> ww_tot; // Total weights

  // Arrays for nodal properties (1D directions)
  std::vector<int> nnode_arr; // Number of nodes per direction
  std::vector<int> nfem_arr;  // FE count per direction
  std::vector<int> ndof_arr;  // DOF count per direction
  std::vector<int> bc_arr;    // Boundary condition type
  std::vector<int> ext_arr;   // Extrapolation method

  // Geometric properties per direction
  std::vector<double> zmin_arr; // Minimum coordinate
  std::vector<double> zmax_arr; // Maximum coordinate
  std::vector<double> zwid_arr; // Width
  std::vector<double> zmid_arr; // Midpoint
  std::vector<double> dz_arr;   // Grid spacing

  std::vector<std::vector<double>> z1d_arr; // Ragged array of 1D coordinates

  // Boundary and DOF indices
  int ntotbc;              // Total boundary nodes
  int ntotdof;             // Total free DOFs
  int ntotfem;             // Total finite elements
  int ntotnode;            // Total nodes
  std::vector<int> idxbc;  // Boundary node indices
  std::vector<int> idxdof; // Free DOF indices

  // Control parameters
  bool nl_debug; // Debug flag

  // Spline storage
  std::vector<double> fval; // Function values
  std::vector<double> fspl; // Spline coefficients

  int rank, size;
  template <typename T>
  void printArray(const std::string &label, const std::vector<T> &arr) const {
    std::cout << std::setw(8) << label << "=";
    for (size_t i = 0; i < std::min<size_t>(20, arr.size()); ++i) {
      std::cout << std::setw(10) << std::fixed << std::setprecision(4)
                << arr[i];
    }
    std::cout << std::endl;
  }

public:
  // Constructor/Destructor
  SplineCubicNdCls();
  ~SplineCubicNdCls() = default;

  // getter
  int getDimension() const { return ndim; }
  double getZMin(int idim) const { return zmin_arr[idim]; }
  double getDz(int idim) const { return dz_arr[idim]; }
  double getZWidth(int idim) const { return zwid_arr[idim]; }
  std::vector<double> getZ1d(int idim) const { return z1d_arr[idim]; }
  int getBC(int idim) const { return bc_arr[idim]; }
  int getExt(int idim) const { return ext_arr[idim]; }
  std::vector<double> getFval() const { return fval; }

  // Initialization methods
  void spc_cls_init_bc();
  void initialize1d(const std::vector<double> &knots);
  void initialize2d(const std::vector<double> &knots_x,
                    const std::vector<double> &knots_y);
  void initialize3d(const std::vector<double> &knots_x,
                    const std::vector<double> &knots_y,
                    const std::vector<double> &knots_z);

  // Boundary condition handling
  void spc_cls_apply_bc_vec(std::vector<double> &f1d);

  // Spline fitting
  void spc_cls_calc_splinefit_ijv1d(int dim, std::vector<int> &i,
                                    std::vector<int> &j,
                                    std::vector<double> &v);

  // Evaluation methods
  double evaluate2d(double x, double y) const;

  // Utility methods

  void spc_cls_init(int input_source, std::optional<int> dim = std::nullopt,
                    std::vector<int> bc = {}, std::vector<int> ext = {},
                    std::vector<int> nodes = {}, std::vector<double> zmin = {},
                    std::vector<double> zmax = {},
                    std::optional<bool> debug = std::nullopt) {
    // control parameters
    nl_debug = true;

    // constants
    norder = 4;
    nintg = 5;
    nintg_tot = nintg * norder;
    xshift = 2.0;

    // Allocate Gauss Quadrature Arrays
    xx_tot.resize(nintg_tot, 0.0);
    ww_tot.resize(nintg_tot, 0.0);
    UtilMath::lgwt(nintg, 0.0, 1.0, xx_tot, ww_tot);
    // Expand Gaussian Quadrature Points
    for (int fic = 1; fic < norder; fic++) {
      for (int j = 0; j < nintg; j++) {
        xx_tot[fic * nintg + j] = xx_tot[j] + static_cast<double>(fic);
        ww_tot[fic * nintg + j] = ww_tot[j];
      }
    }

    // Apply coordinate shift
    for (double &val : xx_tot) {
      val -= xshift;
    }

    // Debug Output
    if (nl_debug && rank == 0) {
      std::cout << "xx = ";
      for (const auto &val : xx_tot) {
        std::cout << std::fixed << std::setprecision(16) << val << " ";
      }
      std::cout << std::endl;

      std::cout << "ww = ";
      for (const auto &val : ww_tot) {
        std::cout << std::fixed << std::setprecision(16) << val << " ";
      }
      std::cout << std::endl;
    }
    // Set up basis functions
    xknot.resize(norder + 1);
    for (int i = 0; i <= norder; i++) {
      xknot[i] = static_cast<double>(i + 1) - 3.0;
    }

    // !--input_source: 0 from arguments; 1 from file "input.spline", namelist
    // "spline_cubicNd_general/_array"
    //     !========Adjustable parameters========
    // Initialize from arguments or file
    if (input_source == 0) {
      // !--bc:0 for periodic; 1 for non-periodic
      if (rank == 0)
        std::cout << "--set spline_cubicNd from arguments--" << std::endl;
      if (dim.has_value()) {
        ndim = dim.value();
      }
      if (debug.has_value()) {
        nl_debug = debug.value();
      }

      // Resize integer arrays
      bc_arr.resize(ndim);
      ext_arr.resize(ndim);
      nnode_arr.resize(ndim);
      nfem_arr.resize(ndim);
      ndof_arr.resize(ndim);

      // Resize double arrays
      zmin_arr.resize(ndim);
      zmax_arr.resize(ndim);
      zwid_arr.resize(ndim);
      dz_arr.resize(ndim);

      // Assign arrays only if provided
      if (!bc.empty())
        bc_arr = bc;
      if (!ext.empty())
        ext_arr = ext;
      if (!nodes.empty())
        nnode_arr = nodes;
      if (!zmin.empty())
        zmin_arr = zmin;
      if (!zmax.empty())
        zmax_arr = zmax;

    } else if (input_source == 1) {
      if (rank == 0)
        std::cout << "--set spline_cubicNd from input file--" << std::endl;
      readInput("input.spline");
    } else {
      if (rank == 0)
        std::cerr << "========Error in SPC_CLS_INIT========" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Derived parameters
    zwid_arr.resize(ndim);
    zmid_arr.resize(ndim);
    dz_arr.resize(ndim);
    nfem_arr.resize(ndim);
    ndof_arr.resize(ndim);

    ntotnode = 1;
    ntotdof = 1;
    ntotfem = 1;

    for (int i = 0; i < ndim; i++) {
      zwid_arr[i] = zmax_arr[i] - zmin_arr[i];
      zmid_arr[i] = (zmax_arr[i] + zmin_arr[i]) / 2.0;
      ndof_arr[i] = nnode_arr[i]; // True for cubic splines

      if (bc_arr[i] == 0) { // Periodic
        nfem_arr[i] = nnode_arr[i];
        dz_arr[i] = zwid_arr[i] / nnode_arr[i];
      } else { // Non-periodic
        nfem_arr[i] = nnode_arr[i] + norder - 2;
        dz_arr[i] = zwid_arr[i] / (nnode_arr[i] - 1);
      }

      ntotnode *= nnode_arr[i];
      ntotdof *= ndof_arr[i];
      ntotfem *= nfem_arr[i];
    }

    ntotbc = ntotfem - ntotdof;

    // Generate 1D arrays
    z1d_arr.resize(ndim);
    for (int i = 0; i < ndim; i++) {
      z1d_arr[i].resize(nnode_arr[i]);
      for (int j = 0; j < nnode_arr[i]; j++) {
        z1d_arr[i][j] = zmin_arr[i] + j * dz_arr[i];
      }

      // !make sure the last value is identical to zmax
      if (bc_arr[i] >= 1) {
        z1d_arr[i][nnode_arr[i] - 1] = zmax_arr[i];
      }
    }

    // Boundary conditions
    initBC();

    // Display info
    showInfo();
  }

  void spc_cls_test() {
    // Allocate memory for fval and fspl
    fval.resize(ntotnode, 0.0);
    fspl.resize(ntotfem, 0.0); // Initialize to zero

    // Set values for fval
    spc_cls_set_ms(fval);

    // Perform spline fitting
    spc_cls_splinefitNd(fval, fspl);

    // Write results to files if rank == 0
    if (rank == 0) {
      UtilIO::write_arr1d_r8(fval, "data_fval.txt", true);
      UtilIO::write_arr1d_r8(fspl, "data_fspl.txt", true);
    }
  }

  

  void initialize1d(const std::vector<double> &x1d,
                    const std::vector<double> &f1d, int bcx, int extx) {
    constexpr int ndim = 1;
    std::array<int, ndim> bc_arr, ext_arr, nnode_arr;
    std::array<double, ndim> zmin_arr, zmax_arr;
    int idim = 0; // Zero-based index in C++
    int input_source = 0;

    // Ensure x1d and f1d have the same size
    if (x1d.size() != f1d.size()) {
      std::cerr << "========error: wrong size in initialize1d========"
                << std::endl;
      return;
    }

    bc_arr[idim] = bcx;
    ext_arr[idim] = extx;
    nnode_arr[idim] = x1d.size();
    zmin_arr[idim] = x1d.front(); // First element of x1d

    if (bc_arr[idim] == 0) {
      zmax_arr[idim] = x1d.back() + (x1d[1] - x1d[0]);
    } else if (bc_arr[idim] >= 1) {
      zmax_arr[idim] = x1d.back();
    }

    // Call the class method spc_cls_init using this-> implicitly
    this->spc_cls_init(input_source, ndim,
                       std::vector<int>(bc_arr.begin(), bc_arr.end()),
                       std::vector<int>(ext_arr.begin(), ext_arr.end()),
                       std::vector<int>(nnode_arr.begin(), nnode_arr.end()),
                       std::vector<double>(zmin_arr.begin(), zmin_arr.end()),
                       std::vector<double>(zmax_arr.begin(), zmax_arr.end()));

    // Resize fspl instead of allocating memory manually
    fspl.resize(ntotfem);

    // Perform spline fitting
    this->spc_cls_splinefitNd(f1d, fspl);

    // Clear idxbc and idxdof instead of deallocation
    idxbc.clear();
    idxdof.clear();
  }

  void evaluate1d(double XX, int idiffx, double &output_fval) const {
    int idim = 0; // Fixed for 1D case
    output_fval = 0.0;

    int ibas1ref, ibas1;
    double x1ref, x1;

    // Convert global coordinate to local
    spc_cls_glb2loc(XX, idim, 0, ibas1ref, x1ref);

    // Summation loop over basis functions
    for (int i1shift = 0; i1shift <= 3; ++i1shift) {
      ibas1 = ibas1ref + i1shift;
      x1 = x1ref - i1shift;

      // Debug output
      std::cout << ibas1 << " " << x1 << " "
                << spc_cls_get_fbas_ix(ibas1, x1, idiffx, idim) << std::endl;

      // Accumulate weighted spline values
      output_fval += spc_cls_get_fbas_ix(ibas1, x1, idiffx, idim) * fspl[ibas1];
    }
  }

  void spc_cls_set_ms(std::vector<double> &fval) {
    fval.assign(fval.size(), 0.0); // Initialize fval to zeros

    if (ndim == 1) {
      for (int fic = 0; fic < nnode_arr[0]; ++fic) {
        fval[fic] = std::cos(TWOPI * dz_arr[0] * fic / zwid_arr[0]);
      }
    } else if (ndim == 2) {
      for (int fic = 0; fic < nnode_arr[0]; ++fic) {
        for (int fjc = 0; fjc < nnode_arr[1]; ++fjc) {
          int idx = fic + fjc * nnode_arr[0];
          fval[idx] = std::cos(TWOPI * dz_arr[0] * fic / zwid_arr[0]) *
                      std::cos(TWOPI * dz_arr[1] * fjc / zwid_arr[1]);
        }
      }
    } else if (ndim == 3) {
      for (int fic = 0; fic < nnode_arr[0]; ++fic) {
        for (int fjc = 0; fjc < nnode_arr[1]; ++fjc) {
          for (int fkc = 0; fkc < nnode_arr[2]; ++fkc) {
            int idx =
                fic + fjc * nnode_arr[0] + fkc * nnode_arr[0] * nnode_arr[1];
            fval[idx] = std::cos(TWOPI * dz_arr[0] * fic / zwid_arr[0]) *
                        std::cos(TWOPI * dz_arr[1] * fjc / zwid_arr[1]) *
                        std::cos(TWOPI * dz_arr[2] * fkc / zwid_arr[2]);
          }
        }
      }
    } else {
      std::cerr << "======== Error: ndim > 3 not implemented in spc_cls_set_ms "
                   "========"
                << std::endl;
    }
  }

  void spc_cls_splinefitNd(const std::vector<double> &fval,
                           std::vector<double> &fspl) {
    if (ndim == 1) {
      spc_cls_splinefit1d(fval, fspl); // Call 1D spline fitting
    } else if (ndim == 2 || ndim == 3) {
      spc_cls_splinefit2d(fval, fspl); // Call 2D or 3D spline fitting
    }
  }
  // To do
  void spc_cls_splinefit1d(const std::vector<double> &fval,
                           std::vector<double> &fspl) {
    // Resize fspl to match the size of fval if needed
    if (fspl.size() != fval.size()) {
      fspl.resize(fval.size());
    }

    // Copy values from fval to fspl
    fspl = fval;
  }
  void spc_cls_splinefit2d(const std::vector<double> &fval,
                           std::vector<double> &fspl) {
    // Resize fspl to match the size of fval if needed
    if (fspl.size() != fval.size()) {
      fspl.resize(fval.size());
    }

    // Copy values from fval to fspl
    fspl = fval;
  }

  void initBC() {
    std::vector<bool> isbc_arr(ndim, false);
    int ith_bc = 0, ith_dof = 0, ith_any = 0;

    // 分配 idxdof 和 idxbc
    idxdof.resize(ntotdof);
    idxbc.resize(ntotbc);

    if (ndim == 1) {
      for (int fic = 1; fic <= nfem_arr[0]; ++fic) {
        isbc_arr[0] = (bc_arr[0] >= 1) && (fic == 1 || fic == nfem_arr[0]);

        // 设置 idx 值
        ith_any++;
        if (isbc_arr[0]) {
          idxbc[ith_bc++] = ith_any;
        } else {
          idxdof[ith_dof++] = ith_any;
        }
      }
    } else if (ndim == 2) {
      for (int fjc = 1; fjc <= nfem_arr[1]; ++fjc) {
        isbc_arr[1] = (bc_arr[1] >= 1) && (fjc == 1 || fjc == nfem_arr[1]);

        for (int fic = 1; fic <= nfem_arr[0]; ++fic) {
          isbc_arr[0] = (bc_arr[0] >= 1) && (fic == 1 || fic == nfem_arr[0]);

          // 设置 idx 值
          ith_any++;
          if (isbc_arr[0] || isbc_arr[1]) {
            idxbc[ith_bc++] = ith_any;
          } else {
            idxdof[ith_dof++] = ith_any;
          }
        }
      }
    } else if (ndim == 3) {
      for (int fkc = 1; fkc <= nfem_arr[2]; ++fkc) {
        isbc_arr[2] = (bc_arr[2] >= 1) && (fkc == 1 || fkc == nfem_arr[2]);

        for (int fjc = 1; fjc <= nfem_arr[1]; ++fjc) {
          isbc_arr[1] = (bc_arr[1] >= 1) && (fjc == 1 || fjc == nfem_arr[1]);

          for (int fic = 1; fic <= nfem_arr[0]; ++fic) {
            isbc_arr[0] = (bc_arr[0] >= 1) && (fic == 1 || fic == nfem_arr[0]);

            // 设置 idx 值
            ith_any++;
            if (isbc_arr[0] || isbc_arr[1] || isbc_arr[2]) {
              idxbc[ith_bc++] = ith_any;
            } else {
              idxdof[ith_dof++] = ith_any;
            }
          }
        }
      }
    } else {
      std::cerr << "======== Error in initBC: ndim > 3 not implemented ========"
                << std::endl;
      return;
    }

    // Debug 信息
    if (nl_debug) {
      std::cout << "idxdof = ";
      for (size_t i = 0; i < std::min(100, ntotdof); ++i) {
        std::cout << idxdof[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "idxbc = ";
      for (size_t i = 0; i < std::min(50, ntotbc); ++i) {
        std::cout << idxbc[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "ith_dof = " << ith_dof << ", ith_bc = " << ith_bc
                << std::endl;
    }
  }

  void showInfo() const {
    if (nl_debug) {
      std::cout << "==================Init SplineNd======================"
                << std::endl;
      std::cout << std::setw(8) << "NDIM=" << std::setw(6) << ndim
                << ", NORDER=" << std::setw(6) << norder
                << ", NINTG=" << std::setw(6) << nintg << std::endl;
      std::cout << std::setw(8) << "NFEM=" << std::setw(6) << ntotfem
                << ", NDOF=" << std::setw(6) << ntotdof
                << ", NNODE=" << std::setw(6) << ntotnode
                << ", NBC=" << std::setw(6) << ntotbc << std::endl;

      for (int fic = 0; fic < ndim; ++fic) {
        std::cout << std::setw(8) << "z1d_arr=";
        for (const auto &val : z1d_arr[fic]) {
          std::cout << std::setw(10) << std::fixed << std::setprecision(4)
                    << val;
        }
        std::cout << std::endl;
      }

      printArray("BC", bc_arr);
      printArray("EXTRAP", ext_arr);
      printArray("Nnode", nnode_arr);
      printArray("Ndof", ndof_arr);
      printArray("Nfem", nfem_arr);

      printArray("Zmin", zmin_arr);
      printArray("Zmax", zmax_arr);
      printArray("Zwid", zwid_arr);
      printArray("Zmid", zmid_arr);
      printArray("dZ", dz_arr);

      std::cout
          << "==================End of Init SplineNd======================"
          << std::endl;
    }

    std::cout << std::setw(20) << "Spline Class size:" << sizeof(*this)
              << " Byte" << std::endl;
  }

  void readInput(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
      std::cerr << "Error: Unable to open " << filename << std::endl;
      return;
    }

    std::string line;
    std::cout << "--------Readinput starts--------" << std::endl;

    // Read general parameters
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      std::string key;
      iss >> key;

      if (key == "ndim") {
        iss >> ndim;
      } else if (key == "nl_debug") {
        int debug_flag;
        iss >> debug_flag;
        nl_debug = (debug_flag != 0);
      } else if (key == "bc_arr") {
        bc_arr.resize(ndim);
        for (int &val : bc_arr)
          iss >> val;
      } else if (key == "ext_arr") {
        ext_arr.resize(ndim);
        for (int &val : ext_arr)
          iss >> val;
      } else if (key == "nnode_arr") {
        nnode_arr.resize(ndim);
        for (int &val : nnode_arr)
          iss >> val;
      } else if (key == "zmin_arr") {
        zmin_arr.resize(ndim);
        for (double &val : zmin_arr)
          iss >> val;
      } else if (key == "zmax_arr") {
        zmax_arr.resize(ndim);
        for (double &val : zmax_arr)
          iss >> val;
      }
    }
    infile.close();

    // Allocate additional arrays
    nfem_arr.resize(ndim);
    ndof_arr.resize(ndim);
    zwid_arr.resize(ndim);
    dz_arr.resize(ndim);

    // Copy values
    std::cout << "ndim: " << ndim << std::endl;
    std::cout << "nl_debug: " << nl_debug << std::endl;

    std::cout << "--------Readinput ends--------" << std::endl;
  }

  void spc_cls_glb2loc(double XX0, int idirect, int ishift, int &ibas,
                       double &xloc) const {
    double XX = XX0 - zmin_arr[idirect]; // Adjust XX to start from 0

    if (bc_arr[idirect] == 0) {
      // Periodic boundary condition
      XX = std::fmod(XX, zwid_arr[idirect]);
      if (XX < 0)
        XX += zwid_arr[idirect]; // Ensure positive modulus

      ibas = static_cast<int>(std::floor(XX / dz_arr[idirect])) + ishift;
      xloc = XX / dz_arr[idirect] - (ibas - 1);
      ibas = (ibas - 1) % nfem_arr[idirect] + 1;

    } else if (bc_arr[idirect] >= 1) {
      // Handle special case at zmax
      if (XX0 == zmax_arr[idirect]) {
        ibas = nnode_arr[idirect];
        xloc = 1.0;
        return;
      }
      
      ibas = static_cast<int>(std::floor(XX / dz_arr[idirect])) + ishift + 1;

      if (ibas <= 0) {
        ibas = 1;
        xloc = (ext_arr[idirect] == 1) ? 0.0 : 1.0; // Out-of-bounds behavior

      } else if (ibas >= nnode_arr[idirect]) {
        ibas = nnode_arr[idirect];
        xloc = (ext_arr[idirect] == 1) ? 2.0 : 1.0; // Out-of-bounds behavior

      } else {
        xloc = XX / dz_arr[idirect] - (ibas - 1) + 1;
      }
    }
  }

  //
  void evaluate3d(double XX, double YY, double ZZ, int idiffx, int idiffy,
                  int idiffz) {
    double x1, x2, x3, fbas1, fbas2, fbas3, x1ref, x2ref, x3ref;
    int ibas1, i1shift, ibas1ref, ibas2ref, ibas3ref;
    int ibas2, i2shift, ibas3, i3shift;
    constexpr int idim1 = 0, idim2 = 1, idim3 = 2;
    int idx;

    // Reset fval array to zero
    std::fill(fval.begin(), fval.end(), 0.0);

    // Convert global to local indices
    spc_cls_glb2loc(ZZ, idim3, 0, ibas3ref, x3ref);
    spc_cls_glb2loc(YY, idim2, 0, ibas2ref, x2ref);
    spc_cls_glb2loc(XX, idim1, 0, ibas1ref, x1ref);

    for (i3shift = 0; i3shift <= 3; i3shift++) {
      ibas3 = ibas3ref + i3shift;
      x3 = x3ref - i3shift;
      fbas3 = spc_cls_get_fbas_ix(ibas3, x3, idiffz, idim3);

      for (i2shift = 0; i2shift <= 3; i2shift++) {
        ibas2 = ibas2ref + i2shift;
        x2 = x2ref - i2shift;
        fbas2 = spc_cls_get_fbas_ix(ibas2, x2, idiffy, idim2);

        for (i1shift = 0; i1shift <= 3; i1shift++) {
          ibas1 = ibas1ref + i1shift;
          x1 = x1ref - i1shift;
          fbas1 = spc_cls_get_fbas_ix(ibas1, x1, idiffx, idim1);

          // Compute index for accessing fspl array
          idx = ibas1 + (ibas2 - 1) * nfem_arr[idim1] +
                (ibas3 - 1) * nfem_arr[idim1] * nfem_arr[idim2];

          // Ensure fval has enough space
          if (idx >= fval.size()) {
            fval.resize(idx + 1, 0.0);
          }

          fval[idx] += fbas1 * fbas2 * fbas3 * fspl[idx];
        }
      }
    }
  }

  // some basic functions
  double spc_cls_get_fbas_ix(int ibas, double x, int idiff, int idirec) const {
    if (idirec >= 3 || idirec < 0) {
      std::cerr << "====Error: wrong idirec in Spc_Cls_Get_Fbas_ix===="
                << std::endl;
      return 0.0;
    }

    int nfem = nfem_arr[idirec];
    if (bc_arr[idirec] == 0) {
      return spc_cls_fbas_in(x, idiff);
    } else if (bc_arr[idirec] >= 1) {
      if (ibas >= 4 && ibas <= nfem - 3) {
        return spc_cls_fbas_in(x, idiff);
      } else if (ibas == 1) {
        return spc_cls_fbas_ed(x, idiff, 1, true);
      } else if (ibas == 2) {
        return spc_cls_fbas_ed(x, idiff, 2, true);
      } else if (ibas == 3) {
        return spc_cls_fbas_ed(x, idiff, 3, true);
      } else if (ibas == nfem) {
        return spc_cls_fbas_ed(x, idiff, 1, false);
      } else if (ibas == nfem - 1) {
        return spc_cls_fbas_ed(x, idiff, 2, false);
      } else if (ibas == nfem - 2) {
        return spc_cls_fbas_ed(x, idiff, 3, false);
      }
    }
    return 0.0;
  }

  double spc_cls_fbas_in(double x, int idiff) const {
    if (x >= -2.0 && x < -1.0) {
      if (idiff == 0) {
        return 4.0 / 3.0 + 2.0 * x + x * x + x * x * x / 6.0;
      } else if (idiff == 1) {
        return 2.0 + 2.0 * x + x * x / 2.0;
      }
    } else if (x >= -1.0 && x < 0.0) {
      if (idiff == 0) {
        return 2.0 / 3.0 - x * x - x * x * x / 2.0;
      } else if (idiff == 1) {
        return -2.0 * x - 1.5 * x * x;
      }
    } else if (x >= 0.0 && x < 1.0) {
      if (idiff == 0) {
        return 2.0 / 3.0 - x * x + x * x * x / 2.0;
      } else if (idiff == 1) {
        return -2.0 * x + 1.5 * x * x;
      }
    } else if (x >= 1.0 && x <= 2.0) {
      if (idiff == 0) {
        return 4.0 / 3.0 - 2.0 * x + x * x - x * x * x / 6.0;
      } else if (idiff == 1) {
        return -2.0 + 2.0 * x - x * x / 2.0;
      }
    }
    return 0.0;
  }

  double spc_cls_fbas_ed(double x, int idiff, int iedge, bool nl_left) const {
    double var = nl_left ? spc_cls_fbas_ed_left(x, idiff, iedge)
                         : spc_cls_fbas_ed_left(-x, idiff, iedge) *
                               (idiff % 2 == 0 ? 1 : -1);
    return var;
  }

  double spc_cls_fbas_ed_left(double x, int idiff, int iedge) const {
    if (iedge == 1) {
      if (x >= 1.0 && x <= 2.0) {
        if (idiff == 0) {
          return 8.0 - 12.0 * x + 6.0 * x * x - x * x * x;
        } else if (idiff == 1) {
          return -12.0 + 12.0 * x - 3.0 * x * x;
        }
      }
    } else if (iedge == 2) {
      if (x >= 0.0 && x < 1.0) {
        if (idiff == 0) {
          return 2.0 * x - 3.0 * x * x + 7.0 * x * x * x / 6.0;
        } else if (idiff == 1) {
          return 2.0 - 6.0 * x + 7.0 * x * x / 2.0;
        }
      } else if (x >= 1.0 && x <= 2.0) {
        if (idiff == 0) {
          return 4.0 / 3.0 - 2.0 * x + x * x - x * x * x / 6.0;
        } else if (idiff == 1) {
          return -2.0 + 2.0 * x - x * x / 2.0;
        }
      }
    } else if (iedge == 3) {
      if (x >= -1.0 && x < 0.0) {
        if (idiff == 0) {
          return 2.0 / 3.0 - x * x - x * x * x / 3.0;
        } else if (idiff == 1) {
          return -2.0 * x - x * x;
        }
      } else if (x >= 0.0 && x < 1.0) {
        if (idiff == 0) {
          return 2.0 / 3.0 - x * x + x * x * x / 2.0;
        } else if (idiff == 1) {
          return -2.0 * x + 1.5 * x * x;
        }
      } else if (x >= 1.0 && x <= 2.0) {
        if (idiff == 0) {
          return 4.0 / 3.0 - 2.0 * x + x * x - x * x * x / 6.0;
        } else if (idiff == 1) {
          return -2.0 + 2.0 * x - x * x / 2.0;
        }
      }
    }
    return 0.0;
  }

private:
  // Internal helper functions
  int spc_cls_loc2glb(const std::vector<int> &loc_idx) const;
  double fbas_ix(int i, double x) const;
  double fbas_in(double x) const;
  double fbas_ed(double x) const;
  double fbas_ed_left(double x) const;
};

// Example implementation of selected methods
SplineCubicNdCls::SplineCubicNdCls() : ndim(0), norder(4), nl_debug(false) {
  // Default constructor initialization
  MPIManager &mpiManager = MPIManager::getInstance(); // 获取MPIManager的实例
  rank = mpiManager.getRank(); // 获取当前进程的rank
  size = mpiManager.getSize(); // 获取总的进程数};
}

void SplineCubicNdCls::initialize1d(const std::vector<double> &knots) {
  xknot = knots;
  // Additional 1D initialization
}

#endif // SPLINECUBICND_H