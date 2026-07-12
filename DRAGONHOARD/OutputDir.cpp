//
//  HDF5Output.cpp
//  DRAGONHOARD
//
//  Created by Bobbie Markwick on 03/07/2026.
//

#include "DragonHoard.hpp"
#include "HDF5_Attrs.hpp"
#include <filesystem>
#include <vector>
#include <cctype>



//MARK: Verify the output directory exists
void DRAGONHOARD::verifyOutputDirectory(){
    const std::filesystem::path& dir = CONFIG::output_dir;
    if (std::filesystem::exists(dir)) {
        if (!std::filesystem::is_directory(dir)) {
            throw std::runtime_error(dir.string() + " exists but is not a directory");
        }
    } else {
        std::filesystem::create_directories(dir);
    }
}
//MARK: Identify the restart file
static bool verifyFileType(const std::string& s, const std::string& ext =  file_ext) {
    std::size_t end = s.rfind('.');
    if (end == std::string::npos) end = s.size(); //File extension missing
    return s.substr(end) == ext;
}
static int extractNumber(const std::string& s) {
    const std::size_t end = s.find('.');
    const std::size_t pos = s.rfind('_', end);
    if (pos == std::string::npos || pos + 1 >= end) return -1; //String isn't in the format [NAME]_#####, skip it
   
    for (std::size_t i = pos + 1; i < end; ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) {
            return -1; //String isn't in the format [NAME]_#####, skip it
        }
    }
    
    return std::stoi(s.substr(pos+1, end-pos-1));
}


std::string DRAGONHOARD::restartFileName(){
    std::string filename = "";
    int max_frame_num = -1;
    
    const std::filesystem::path& dir = CONFIG::output_dir;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        auto str = entry.path().filename().string();
        if(!verifyFileType(str)) continue;
        
        int n = extractNumber(str);
        #if RESTART_FRAME < 0
        if(n > max_frame_num){
            filename = str;
            max_frame_num = n;
        }
        #else
        if (n==RESTART_FRAME) return str;
        #endif
    }
    return filename;
}
