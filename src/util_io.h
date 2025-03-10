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

#ifndef UTIL_IO_H
#define UTIL_IO_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

class UtilIO {
public:
    static void write_arr1d_r8(const std::vector<double>& arr1d, 
                               const std::string& sfile, 
                               bool nl_showinfo = false);
};

#endif // UTIL_IO_H
