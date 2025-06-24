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
#ifndef SPLINE_CUBIC_2D1F_H
#define SPLINE_CUBIC_2D1F_H

#include "SplineCubicNd.h"
#include <complex>

// 2D1F: 2d in (r,theta), Fourier (1f) in phi
class SplineCubic2d1f {
private:
  SplineCubicNdCls spline2d; // 组合一个 Nd 版本的样条类
  std::vector<int> ntor1d;   // toroidal modes
  int ndim = 3;              // Number of dimensions
  int norder = 4;            // Spline order (e.g., 4 for cubic)
  int nintg;                 // Integration points
  int nintg_tot;             // Total integration points
  double xshift;             // Coordinate shift

  std::vector<double> xknot;  // Knot vector
  std::vector<double> xx_tot; // Total coordinates
  std::vector<double> ww_tot; // Total weights

  // Arrays for nodal properties (1D directions)
  std::vector<int> nnode_arr; // Number of nodes per direction
  std::vector<int> nfem_arr;  // FE count per direction
  std::vector<int> ndof_arr;  // DOF count per direction
  std::vector<int> bc_arr;    // Boundary condition type

  // Geometric properties per direction
  std::vector<double> zmin_arr; // Minimum coordinate
  std::vector<double> zmax_arr; // Maximum coordinate
  std::vector<double> zwid_arr; // Width
  std::vector<double> zmid_arr; // Midpoint
  std::vector<double> dz_arr;   // Grid spacing

  std::vector<std::vector<double>> z1d_arr; // Ragged array of 1D coordinates
  std::vector<double> z1d1, z1d2, z1d3;

  // Boundary and DOF indices
  int ntotbc;              // Total boundary nodes
  int ntotdof;             // Total free DOFs
  int ntotfem;             // Total finite elements
  int ntotnode;            // Total nodes
  std::vector<int> idxbc;  // Boundary node indices
  std::vector<int> idxdof; // Free DOF indices

  // index of boundary and free nodes
  int ntot12fem, ntot12dof;

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
  // geters
  const std::vector<int> &getNtor1d() const { return ntor1d; }
  const std::vector<int> &getNnodeArr() const { return nnode_arr; }
  const std::vector<double> &get_z1d1() const { return z1d1; }
  const std::vector<double> &get_z1d2() const { return z1d2; }
  const std::vector<double> &get_z1d3() const { return z1d3; }
  const int get_ntot12fem() const { return ntot12fem; }
  const int get_ntot12dof() const { return ntot12dof; }
  const int get_ntotdof() const { return ntotdof; }
  const int get_ntotfem() const { return ntotfem; }
  // Constructor, initializing the spline2d object
  SplineCubic2d1f(); // Declaration of the default constructor
  SplineCubic2d1f(int input_source, const int bc_arr[3], const int nnode_arr[3],
                  const double zmin_arr[3], const double zmax_arr[3],
                  bool nl_debug);

  void spc_cls_calc_idxdof2d1f(int lenntor, std::vector<int> &idxdof2d1f);

  // SHM: f1d(1:ntotfem); assume: 1--rad, 2--the, 3--phi
  // prad1d,pthe1d,pphi1d, 1d array; idiff:
  // (d^i1/drad^i1,d^i2/dthe^i2,d^i3/dphi^i3) note output pf1d is real*8,
  // pf=2*real(f_{k/ne0}*exp{inphi})+real(f_{k=0})
  void spc_cls_sp2p2d1f(const std::vector<std::complex<double>> &f1d,
                        const std::vector<int> &ntor1d,
                        const std::vector<double> &prad1d,
                        const std::vector<double> &pthe1d,
                        const std::vector<double> &pphi1d,
                        std::vector<double> &pf1d,
                        const std::array<int, 3> &idiff) {

    int dimrad{0}, dimthe{1};
    pf1d.assign(pf1d.size(), 0.0);
    int np = prad1d.size();
    int lenntor = ntor1d.size();

    for (int fpc = 0; fpc < np; ++fpc) {
      double rad1 = prad1d[fpc];
      if (rad1 > zmax_arr[0] || rad1 < zmin_arr[0])
        continue;

      double the1 = UtilMath::modulo(pthe1d[fpc] - zmin_arr[1], zwid_arr[1]) +
                    zmin_arr[1];
      double phi1 = UtilMath::modulo(pphi1d[fpc] - zmin_arr[2], zwid_arr[2]) +
                    zmin_arr[2];

      int ithefloorhere, iradfloorhere; // from 1 to ngrids
      double xthe;
      std::complex<double> pfradthe_c16 = 0.0;
      std::complex<double> pfrad_c16 = 0.0;
      double xrad;
      for (int itor = 0; itor < lenntor; ++itor) {
        pfradthe_c16 = 0.0;
        for (int itheshift = 0; itheshift <= 3; ++itheshift) {
          spc_cls_glb2loc(the1, dimthe, itheshift, ithefloorhere, xthe);
          // std::cout << "ithefloorhere=" << ithefloorhere << std::endl;
          pfrad_c16 = 0.0;
          if (ithefloorhere >= 1) {
            for (int iradshift = 0; iradshift <= 3; ++iradshift) {
              spc_cls_glb2loc(rad1, dimrad, iradshift, iradfloorhere, xrad);

              // std::cout << "iradfloorhere=" << iradfloorhere << std::endl;
              if (iradfloorhere >= 1) {
                int idx123 = iradfloorhere + nfem_arr[0] * (ithefloorhere - 1) +
                             ntot12fem * (itor);

                pfrad_c16 +=
                    spc_cls_get_fbas_ix(iradfloorhere, xrad, idiff[0], dimrad) *
                    f1d[idx123 - 1];
              } // inside rad
            }   // iradshift
            pfradthe_c16 += pfrad_c16 * spc_cls_get_fbas_ix(ithefloorhere, xthe,
                                                            idiff[1], dimthe);
          } // inside the
        }   // itheshift

        double f1or2 = (ntor1d[itor] != 0) ? 2.0 : 1.0;
        std::complex comp =
            std::pow(std::complex<double>(0.0, ntor1d[itor]), idiff[2]);
        pf1d[fpc] += f1or2 * std::real(pfradthe_c16 * comp *
                                       std::exp(std::complex<double>(
                                           0.0, 1.0 * ntor1d[itor] * phi1)));
      } // itor
    }   // fpc

    double yfac =
        std::pow(dz_arr[0], -idiff[0]) * std::pow(dz_arr[1], -idiff[1]);
    for (double &val : pf1d) {
      val *= yfac;
    }
  }

  void spc_cls_sp2p2d1f_grad(const std::vector<std::complex<double>> &f1d,
                             const std::vector<int> &ntor1d,
                             const std::vector<double> &prad1d,
                             const std::vector<double> &pthe1d,
                             const std::vector<double> &pphi1d,
                             std::vector<double> &pdfdrad1d,
                             std::vector<double> &pdfdthe1d,
                             std::vector<double> &pdfdphi1d) {
    int dimrad{0}, dimthe{1};
    std::vector<double> fNrad0_arr(4, 0.0), fNrad1_arr(4, 0.0);
    double fNthe0, fNthe1;

    pdfdrad1d.assign(pdfdrad1d.size(), 0.0);
    pdfdthe1d.assign(pdfdthe1d.size(), 0.0);
    pdfdphi1d.assign(pdfdphi1d.size(), 0.0);

    int np = prad1d.size();
    int lenntor = ntor1d.size();

    for (int fpc = 0; fpc < np; ++fpc) {
      double rad1 = prad1d[fpc];
      if (rad1 < zmin_arr[0] || rad1 > zmax_arr[0])
        continue;

      double the1 = UtilMath::modulo(pthe1d[fpc] - zmin_arr[1], zwid_arr[1]) +
                    zmin_arr[1];
      double phi1 = UtilMath::modulo(pphi1d[fpc] - zmin_arr[2], zwid_arr[2]) +
                    zmin_arr[2];

      for (int itor = 0; itor < lenntor; ++itor) {
        std::complex<double> pfradthe_c16[3] = {0.0, 0.0, 0.0};
        for (int itheshift = 0; itheshift <= 3; ++itheshift) {
          int ithefloorhere;
          double xthe;
          spc_cls_glb2loc(the1, dimthe, itheshift, ithefloorhere, xthe);
          std::complex<double> pfrad_c16[3] = {0.0, 0.0, 0.0};
          if (ithefloorhere >= 1) { // edge marker resides in <=3 elements
            for (int iradshift = 0; iradshift <= 3; ++iradshift) {
              int iradfloorhere;
              double xrad;
              spc_cls_glb2loc(rad1, dimrad, iradshift, iradfloorhere, xrad);

              if (iradfloorhere >= 1) {
                int idx123 = iradfloorhere + nfem_arr[0] * (ithefloorhere - 1) +
                             ntot12fem * (itor);
                //  --sp2p versus p2sp: need switch pf1d & f1d
                if (itheshift == 0) {
                  fNrad0_arr[iradshift] =
                      spc_cls_get_fbas_ix(iradfloorhere, xrad, 0, dimrad);
                  fNrad1_arr[iradshift] =
                      spc_cls_get_fbas_ix(iradfloorhere, xrad, 1, dimrad);
                }
                pfrad_c16[0] += fNrad1_arr[iradshift] * f1d[idx123 - 1];
                pfrad_c16[1] += fNrad0_arr[iradshift] * f1d[idx123 - 1];
                // pfrad_c16(2) is same as pfrad_c16(1)
              } // inside iradfloorhere >= 1
            }   // iradshift
            pfrad_c16[2] = pfrad_c16[1];
            fNthe0 = spc_cls_get_fbas_ix(ithefloorhere, xthe, 0, dimthe);
            fNthe1 = spc_cls_get_fbas_ix(ithefloorhere, xthe, 1, dimthe);
            pfradthe_c16[0] += pfrad_c16[0] * fNthe0;
            pfradthe_c16[1] += pfrad_c16[1] * fNthe1;
            pfradthe_c16[2] += pfrad_c16[2] * fNthe0;
          } // inside ithefloorhere >= 1
        }   // itheshift

        double f1or2 = (ntor1d[itor] != 0) ? 2.0 : 1.0;
        std::complex<double> fNphi0_c16 =
            std::exp(std::complex<double>(0.0, 1.0 * ntor1d[itor] * phi1));
        std::complex<double> fNphi1_c16 =
            fNphi0_c16 * std::complex<double>(0.0, ntor1d[itor]);

        pdfdrad1d[fpc] += f1or2 * std::real(pfradthe_c16[0] * fNphi0_c16);
        pdfdthe1d[fpc] += f1or2 * std::real(pfradthe_c16[1] * fNphi0_c16);
        pdfdphi1d[fpc] += f1or2 * std::real(pfradthe_c16[2] * fNphi1_c16);
      } // itor
    }   // fpc

    // Convert division to multiplication by the reciprocal
    double reciprocal = 1.0 / dz_arr[0];
    std::transform(pdfdrad1d.begin(), pdfdrad1d.end(), pdfdrad1d.begin(),
                   [reciprocal](double val) { return val * reciprocal; });
    reciprocal = 1.0 / dz_arr[1];
    // Use std::transform to multiply each element by reciprocal
    std::transform(pdfdthe1d.begin(), pdfdthe1d.end(), pdfdthe1d.begin(),
                   [reciprocal](double val) { return val * reciprocal; });
    // no factor for pdfdphi1d
  }
  // --------Utility subroutines--------
  void spc_cls_glb2loc(double XX0, int idirect, int ishift, int &ibas,
                       double &xloc) const {
    // Reset XX to [0, zmax)
    // ibas from 1 to nfem_arr[idirect]
    double XX = XX0 - zmin_arr[idirect];

    if (bc_arr[idirect] == 0) {
      XX = UtilMath::modulo(
          XX, zwid_arr[idirect]); // Equivalent to Fortran's modulo()

      ibas = static_cast<int>(std::floor(XX / dz_arr[idirect])) + ishift;
      xloc = XX / dz_arr[idirect] - (ibas - 1);
      ibas = UtilMath::modulo((ibas - 1), nfem_arr[idirect]) + 1;

    } else if (bc_arr[idirect] == 1) {
      if (XX >= 0 && XX < zwid_arr[idirect]) {
        ibas = static_cast<int>(std::floor(XX / dz_arr[idirect])) + ishift + 1;
        xloc = XX / dz_arr[idirect] - (ibas - 1) + 1;

        if (ibas > nfem_arr[idirect]) {
          ibas = nfem_arr[idirect];
          xloc = 100;
        }
      } else { // Make N to be 0
        ibas = 1;
        xloc = -100;
      }
    }
  }
  double spc_cls_get_fbas_ix(int ibas, double x, int idiff, int idirec) {
    // idirec: 0--rad, 1--the
    // ibas: 1,2,3,...,nfem_arr[idirec]
    // idiff: 0--fbas, 1--dfbas/d coord
    if (idirec >= 3 || idirec < 0) {
      std::cerr << "====Error: SplineCubic2d1f: wrong idirec in "
                   "Spc_Cls_Get_Fbas_ix===="
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

  double spc_cls_fbas_in(double x, int idiff) {
    if (x >= -2.0 && x < -1.0) {
      return (idiff == 0) ? (4.0 / 3.0 + 2 * x + x * x + x * x * x / 6.0)
                          : (2.0 + 2.0 * x + x * x / 2.0);
    } else if (x >= -1.0 && x < 0.0) {
      return (idiff == 0) ? (2.0 / 3.0 - x * x - x * x * x / 2.0)
                          : (-2.0 * x - 1.5 * x * x);
    } else if (x >= 0.0 && x < 1.0) {
      return (idiff == 0) ? (2.0 / 3.0 - x * x + x * x * x / 2.0)
                          : (-2.0 * x + 1.5 * x * x);
    } else if (x >= 1.0 && x <= 2.0) {
      return (idiff == 0) ? (4.0 / 3.0 - 2 * x + x * x - x * x * x / 6.0)
                          : (-2.0 + 2.0 * x - x * x / 2.0);
    }
    return 0.0;
  }

  double spc_cls_fbas_ed(double x, int idiff, int iedge, bool nl_left) {
    double var = nl_left
                     ? spc_cls_fbas_ed_left(x, idiff, iedge)
                     : spc_cls_fbas_ed_left(-x, idiff, iedge) * pow(-1, idiff);
    return var;
  }

  double spc_cls_fbas_ed_left(double x, int idiff, int iedge) {
    if (iedge == 1) {
      if (x >= 1.0 && x <= 2.0) {
        return (idiff == 0) ? (8.0 - 12.0 * x + 6.0 * x * x - x * x * x)
                            : (-12.0 + 12.0 * x - 3.0 * x * x);
      }
    } else if (iedge == 2) {
      if (x >= 0.0 && x < 1.0) {
        return (idiff == 0) ? (2.0 * x - 3.0 * x * x + 7.0 * x * x * x / 6.0)
                            : (2.0 - 6.0 * x + 7.0 * x * x / 2.0);
      } else if (x >= 1.0 && x <= 2.0) {
        return (idiff == 0) ? (4.0 / 3.0 - 2.0 * x + x * x - x * x * x / 6.0)
                            : (-2.0 + 2.0 * x - x * x / 2.0);
      }
    } else if (iedge == 3) {
      if (x >= -1.0 && x < 0.0) {
        return (idiff == 0) ? (2.0 / 3.0 - x * x - x * x * x / 3.0)
                            : (-2.0 * x - x * x);
      } else if (x >= 0.0 && x < 1.0) {
        return (idiff == 0) ? (2.0 / 3.0 - x * x + x * x * x / 2.0)
                            : (-2.0 * x + 3.0 * x * x / 2.0);
      } else if (x >= 1.0 && x <= 2.0) {
        return (idiff == 0) ? (4.0 / 3.0 - 2.0 * x + x * x - x * x * x / 6.0)
                            : (-2.0 + 2.0 * x - x * x / 2.0);
      }
    }
    return 0.0;
  }
  // --------end of Utility subroutines--------
};

#endif // SPLINE_CUBIC_2D1F_H
