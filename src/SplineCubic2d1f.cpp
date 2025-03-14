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
// Default constructor (just an empty implementation or whatever initialization is needed)
SplineCubic2d1f::SplineCubic2d1f() {
    // Empty constructor or you can initialize default values
    // Example: spline2d = SomeDefaultInitialization;
    std::cout << "SplineCubic2d1f default constructor called" << std::endl;  // For debugging
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
  spline2d.spc_cls_init(input_source, spline2d_ndim, bc_xy_arr, ext_xy_arr,
                        nodes_xy_arr, zmin_xy_arr, zmax_xy_arr, nl_debug);
  spline2d.initialize2d(spline2d.getZ1d(0), spline2d.getZ1d(1),
                        spline2d.getFval(), spline2d.getBC(0),
                        spline2d.getBC(1), spline2d.getExt(0),
                        spline2d.getExt(1));
}