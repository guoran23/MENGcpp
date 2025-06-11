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

#include <array>

template <typename T, std::size_t N>
void print_vector(const std::string& name, const std::array<T, N>& arr) {
    std::cout << name << " = [";
    for (size_t i = 0; i < N; ++i) {
        std::cout << arr[i];
        if (i < N - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

template<typename T>
void print_vector(const std::string& name, const std::vector<T>& vec) {
    std::cout << name << " = [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

class UtilIO {
public:
    static void write_arr1d_r8(const std::vector<double>& arr1d, 
                               const std::string& sfile, 
                               bool nl_showinfo = false);
};

#endif // UTIL_IO_H
