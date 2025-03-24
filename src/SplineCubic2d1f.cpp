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
#include "SplineCubic2d1f.h"
// Default constructor (just an empty implementation or whatever initialization
// is needed)
SplineCubic2d1f::SplineCubic2d1f() {
  // Empty constructor or you can initialize default values
  // Example: spline2d = SomeDefaultInitialization;
  std::cout << "SplineCubic2d1f default constructor called"
            << std::endl; // For debugging
}
// Constructor, initializing the spline2d object
SplineCubic2d1f::SplineCubic2d1f(int input_source, const int bc_arr[3],
                                 const int nnode_arr[3],
                                 const double zmin_arr[3],
                                 const double zmax_arr[3], bool nl_debug) {
  // Initialize vectors using the first two elements of each array
  std::vector<int> bc_xy_arr(bc_arr, bc_arr + 2);
  std::vector<int> nodes_xy_arr(nnode_arr, nnode_arr + 2);
  std::vector<double> zmin_xy_arr(zmin_arr, zmin_arr + 2);
  std::vector<double> zmax_xy_arr(zmax_arr, zmax_arr + 2);

  int spline2d_ndim = 2;
  std::vector<int> ext_xy_arr = {0, 0}; // r,the Extrapolation method
  // Initialize the spline2d object
  spline2d.spc_cls_init(input_source, spline2d_ndim, bc_xy_arr, ext_xy_arr,
                        nodes_xy_arr, zmin_xy_arr, zmax_xy_arr, nl_debug);
  spline2d.initialize2d(spline2d.getZ1d(0), spline2d.getZ1d(1),
                        spline2d.getFval(), spline2d.getBC(0),
                        spline2d.getBC(1), spline2d.getExt(0),
                        spline2d.getExt(1));

  // Set the derived parameters
  // Assuming these member variables are in the class and are of type
  // std::vector<double>
  this->zmin_arr = std::vector<double>(zmin_arr, zmin_arr + 3);
  this->zmax_arr = std::vector<double>(zmax_arr, zmax_arr + 3);
  this->zwid_arr = std::vector<double>(3, 0.0);
  this->zmid_arr = std::vector<double>(3, 0.0);
  this->bc_arr = std::vector<int>(bc_arr, bc_arr + 3);
  this->nnode_arr = std::vector<int>(nnode_arr, nnode_arr + 3);
  this->nfem_arr = std::vector<int>(3, 0);
  this->ndof_arr = std::vector<int>(3, 0);
  this->dz_arr = std::vector<double>(3, 0.0);
  this->ntotdof = 0;
  this->ntotfem = 0;
  this->ntotbc = 0;
  this->ntot12fem = 0;
  this->ntot12dof = 0;
  this->nl_debug = nl_debug;
  for (size_t i = 0; i < this->zmax_arr.size(); ++i) {
    this->zwid_arr[i] = this->zmax_arr[i] - this->zmin_arr[i];
    this->zmid_arr[i] = (this->zmax_arr[i] + this->zmin_arr[i]) / 2.0;
  }

  // Resizing ndof_arr to the size of nnode_arr and assigning values
  this->ndof_arr = std::vector<int>(nnode_arr, nnode_arr + 3);

  // Initialize ntotdof and ntotfem
  this->ntotdof = 1;
  this->ntotfem = 1;

  //   Loop over ndim
  std::cout << "ndim=" << this->ndim << std::endl;
  std::cout << "nnode_arr[0]=" << this->nnode_arr[0] << std::endl;
  std::cout << "bc_arr[0]=" << this->bc_arr[0] << std::endl;
  std::cout << "bc_arr[1]=" << this->bc_arr[1] << std::endl;
  std::cout << "norder=" << this->norder << std::endl;
  for (int fic = 0; fic < this->ndim; ++fic) {
    if (this->bc_arr[fic] == 0) {
      this->nfem_arr[fic] = this->nnode_arr[fic];
      this->dz_arr[fic] = this->zwid_arr[fic] / this->nnode_arr[fic];
    } else {
      this->nfem_arr[fic] = this->nnode_arr[fic] + this->norder - 2;
      this->dz_arr[fic] = this->zwid_arr[fic] / (this->nnode_arr[fic] - 1);
    }

    this->ntotdof *= this->ndof_arr[fic];
    this->ntotfem *= this->nfem_arr[fic];
  }

  this->ntotbc = this->ntotfem - this->ntotdof;
  this->ntot12fem = this->nfem_arr[0] * this->nfem_arr[1];
  this->ntot12dof = this->ndof_arr[0] * this->ndof_arr[1];
  std::cout << "bc_arr = " << this->bc_arr[0] << " " << this->bc_arr[1]
            << " " << this->bc_arr[2] << std::endl;
  std::cout << "nnode_arr = " << this->nnode_arr[0] << " " << this->nnode_arr[1]
            << " " << this->nnode_arr[2] << std::endl;
  std::cout << "nfem_arr = " << this->nfem_arr[0] << " " << this->nfem_arr[1]
            << " " << this->nfem_arr[2] << std::endl;
  std::cout << "ntot12fem = " << this->ntot12fem << std::endl;
  std::cout << "ntot12dof = " << this->ntot12dof << std::endl;
}

void SplineCubic2d1f::spc_cls_calc_idxdof2d1f(int lenntor,
                                              std::vector<int> &idxdof2d1f) {
  int ith_bc = 0, ith_dof = 0, ith_any = 0;

  if (idxdof2d1f.size() != ntot12dof * lenntor) {
    std::cerr << "----error of idxdof2d1f size in spc_cls_calc_idxdof2d1f----"
              << std::endl;
    return;
  }

  for (int fkc = 1; fkc <= lenntor; ++fkc) {
    std::vector<bool> isbc_arr(3, false);

    for (int fjc = 1; fjc <= nfem_arr[1]; ++fjc) {
      if (bc_arr[1] == 1 && (fjc == 1 || fjc == nfem_arr[1])) {
        isbc_arr[1] = true;
      } else {
        isbc_arr[1] = false;
      }

      for (int fic = 1; fic <= nfem_arr[0]; ++fic) {
        if (bc_arr[0] == 1 && (fic == 1 || fic == nfem_arr[0])) {
          isbc_arr[0] = true;
        } else {
          isbc_arr[0] = false;
        }

        // Set idx values
        ith_any++;
        if (isbc_arr[0] || isbc_arr[1] || isbc_arr[2]) {
          ith_bc++;
          // idxbc2d1f[ith_bc] = ith_any;  // Uncomment if using idxbc2d1f
        } else {
          ith_dof++;
          idxdof2d1f[ith_dof - 1] = ith_any; // Adjust for 0-based indexing
        }
      }
    }
  }

  if (nl_debug) {
    std::cout << "idxdof2d1f=" << std::endl;
    for (int i = 0; i < std::min(20, ntot12dof * lenntor); ++i) {
      std::cout << idxdof2d1f[i] << " ";
    }
    std::cout << std::endl;
  }
}