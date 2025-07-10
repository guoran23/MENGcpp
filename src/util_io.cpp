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

#include "util_io.h"

std::vector<double> parseArray(const std::string& input) {
    std::vector<double> result;
    std::stringstream ss(input);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
        try {
            result.push_back(std::stod(token));
        } catch (...) {
            std::cerr << "Invalid number in array: " << token << "\n";
        }
    }
    return result;
}

void UtilIO::write_arr1d_r8(const std::vector<double>& arr1d, 
                            const std::string& sfile, 
                            bool nl_showinfo) {
    std::size_t narr = arr1d.size(); // Get array size

    // Print info if requested
    if (nl_showinfo) {
        std::cout << "----write data to " << sfile << ", len= " << narr << std::endl;
    }

    // Open file for writing
    std::ofstream outfile(sfile);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << sfile << " for writing." << std::endl;
        return;
    }

    // Write data to file
    for (const auto& value : arr1d) {
        outfile << std::scientific << std::setprecision(3) << std::setw(16) << value << std::endl;
    }

    outfile.close(); // Close file
}
