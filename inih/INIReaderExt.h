#pragma once

#include "INIReader.h"
#include <vector>
#include <string>
#include <sstream>

class INIReaderExt : public INIReader {
public:
    using INIReader::INIReader;

    std::vector<double> GetRealList(const std::string& section,
                                    const std::string& key,
                                    const std::vector<double>& default_val) const {
        std::string raw = Get(section, key, "");
        if (raw.empty()) return default_val;

        std::vector<double> result;
        std::stringstream ss(raw);
        std::string token;

        while (std::getline(ss, token, ',')) {
            try {
                result.push_back(std::stod(token));
            } catch (...) {
                // 跳过无效项
            }
        }

        return result;
    }
    
    std::vector<int> GetIntList(const std::string& section,
                            const std::string& key,
                            const std::vector<int>& default_val) const {
    std::string raw = Get(section, key, "");
    if (raw.empty()) return default_val;

    std::vector<int> result;
    std::stringstream ss(raw);
    std::string token;

    while (std::getline(ss, token, ',')) {
        try {
            result.push_back(std::stoi(token));
        } catch (...) {
            // 跳过无效项
        }
    }

    return result;
}

};
