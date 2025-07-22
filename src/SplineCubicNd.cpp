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
#include "SplineCubicNd.h"
// Example implementation of selected methods
SplineCubicNdCls::SplineCubicNdCls() : ndim(0), norder(4), nl_debug(false) {
    // Default constructor initialization
    MPIManager &mpiManager = MPIManager::getInstance(); // 获取MPIManager的实例
    rank = mpiManager.getRank(); // 获取当前进程的rank
    size = mpiManager.getSize(); // 获取总的进程数
  }
  
  void SplineCubicNdCls::initialize1d(const std::vector<double> &knots) {
    xknot = knots;
    // Additional 1D initialization
  }
  