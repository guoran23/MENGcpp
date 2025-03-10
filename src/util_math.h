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

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace UtilMath {
    int modulo(int a, int p);
    void lgwt(int n, double a, double b, std::vector<double> &x1d, std::vector<double> &w1d);
}


#endif // UTIL_MATH_H
