#include "util_io.h"

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
