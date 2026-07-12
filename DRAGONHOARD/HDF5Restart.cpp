//
//  HDF5Output.cpp
//  DRAGON/IO
//
//  Created by Bobbie Markwick on 03/07/2026.
//

#include "DragonHoard.hpp"
#include "HDF5_Attrs.hpp"
#include "Config.h"
#include <H5Cpp.h>
#include <vector>
#include <format>
#include <stdexcept>

//MARK: Helpers
namespace{

std::vector<double> readArray(H5::H5File& file, const std::string& key) {
    H5::DataSet dataset = file.openDataSet(key);
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    std::vector<hsize_t> dims(rank);
    dataspace.getSimpleExtentDims(dims.data());

    size_t total_size = 1;
    for (hsize_t d : dims) total_size *= static_cast<size_t>(d);

    std::vector<double> data(total_size);
    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    return data;
}

int readIntAttribute(H5::H5File& file, const std::string& name) {
    int value;
    H5::Attribute attr = file.openAttribute(name);
    attr.read(H5::PredType::NATIVE_INT, &value);
    return value;
}

double readDoubleAttribute(H5::H5File& file, const std::string& name) {
    double value;
    H5::Attribute attr = file.openAttribute(name);
    attr.read(H5::PredType::NATIVE_DOUBLE, &value);
    return value;
}

}
std::string checkExtension(const std::string& filename) {
    if (filename.size() <= file_ext.size()) return filename + file_ext; //Shorter = definitley missing
    if (filename.substr(filename.size() - file_ext.size()) == file_ext) return filename;
    return filename + file_ext;
}

//MARK: Dispatch
void DRAGONHOARD::loadFromFile(Grid& grid, double& t, int& cycle, const std::string& filename){
    Grid3D* grid3D = dynamic_cast<Grid3D*>(&grid);
    if(grid3D){
        loadFromFile(*grid3D, t, cycle, filename);
        return;
    }
    Grid2D* grid2D = dynamic_cast<Grid2D*>(&grid);
    if(grid2D){
        loadFromFile(*grid2D, t, cycle, filename);
        return;
    }
    Grid1D* grid1D = dynamic_cast<Grid1D*>(&grid);
    if(grid1D){
        loadFromFile(*grid1D, t, cycle, filename);
        return;
    }
    
    throw std::runtime_error("IO attempted with unsupported Grid type");
}

//MARK: 1D
void DRAGONHOARD::loadFromFile(Grid1D& grid, double& t, int& cycle, const std::string& filename){
    std::string path = CONFIG::output_dir + "/" + checkExtension(filename);
    H5::H5File file(path, H5F_ACC_RDONLY);
    
    //Verify compatibility
    const int fmt = readIntAttribute(file, key_fmt);
    if(fmt > 1) throw std::runtime_error("Please update to the latest version of DRAGON to read this file");
    const int dim = readIntAttribute(file, key_dim);
    if(dim != 1) throw std::runtime_error(std::format("Expected 1 dimension but got {}",dim));
    const bool mhd = readIntAttribute(file, key_mhd);
    #ifndef MHD
    if(mhd) throw std::runtime_error("File was saved using MHD. Please enable MHD in Config.h");
    #endif

    
    //Extract size
    const int nx = readIntAttribute(file, key_nx);
    if(nx != grid.getSize())
        throw std::runtime_error("File size does not match the grid size");
    int ng = readIntAttribute(file, key_ng);
    if(ng > grid.getGhosts()) ng = grid.getGhosts();
    const int i0 = -ng, in = nx+ng;

    //Other metadata
    grid.dx = readDoubleAttribute(file, key_dx);
    t = readDoubleAttribute(file, key_time);
    cycle = readIntAttribute(file, key_cyc);
    
    const int write_opt = readIntAttribute(file, key_wrt_opt);
    
    //Write-option-indepependent values
    std::vector<double> rho = readArray(file, key_rho);
    #ifdef MHD
    std::vector<double> Bx = mhd ? readArray(file, key_Bx) : std::vector<double>();
    std::vector<double> By = mhd ? readArray(file, key_By) : std::vector<double>();
    std::vector<double> Bz = mhd ? readArray(file, key_Bz) : std::vector<double>();
    #endif
    for(int i = i0; i<in; i++){
        size_t n = (i-i0);
        grid[i].rho = rho[n];
        #ifdef MHD
        if(mhd) {
            grid[i].B = {Bx[n],By[n],Bz[n]};
        } else {
            grid[i].B = {0,0,0};
        }
        #endif
    }
    
    if(write_opt & HDF5_WRITE_PRIMITIVE) { //Load from Primitive Variables
        std::vector<double> p = readArray(file, key_p);
        std::vector<double> vx = readArray(file, key_vx);
        std::vector<double> vy = readArray(file, key_vy);
        std::vector<double> vz = readArray(file, key_vz);
        for(int i = i0; i<in; i++){
            size_t n = (i-i0);
            grid[i].p = p[n];
            grid[i].v = {vx[n],vy[n],vz[n]};
        }
    } else { //Load from Conservative Variables
        std::vector<double> E = readArray(file, key_E);
        std::vector<double> px = readArray(file, key_px);
        std::vector<double> py = readArray(file, key_py);
        std::vector<double> pz = readArray(file, key_pz);
        for(int i = i0; i<in; i++){
            size_t n = (i-i0);
            ConservativeState U;
            U.rho = grid[i].rho;
            U.p = {px[n],py[n],pz[n]};
            U.E = E[n];
            #ifdef MHD
            U.B = grid[i].B;
            #endif
            grid[i] = U;
        }
    }
}

//MARK: 2D
void DRAGONHOARD::loadFromFile(Grid2D& grid, double& t, int& cycle, const std::string& filename){
    std::string path = CONFIG::output_dir + "/" + checkExtension(filename);
    H5::H5File file(path, H5F_ACC_RDONLY);
    
    //Verify compatibility
    const int fmt = readIntAttribute(file, key_fmt);
    if(fmt > 1) throw std::runtime_error("Please update to the latest version of DRAGON to read this file");
    const int dim = readIntAttribute(file, key_dim);
    if(dim != 2) throw std::runtime_error(std::format("Expected 2 dimensions but got {}",dim));
    const bool mhd = readIntAttribute(file, key_mhd);
    #ifndef MHD
    if(mhd) throw std::runtime_error("File was saved using MHD. Please enable MHD in Config.h");
    #endif
    
    //Extract size
    const int nx = readIntAttribute(file, key_nx);
    const int ny = readIntAttribute(file, key_ny);
    if(nx != grid.getSizeX() || ny != grid.getSizeY())
        throw std::runtime_error("File size does not match the grid size");
    int ng = readIntAttribute(file, key_ng);
    if(ng > grid.getGhosts()) ng = grid.getGhosts();
    const int i0 = -ng, in = nx+ng, j0 = -ng, jn = ny+ng;

    //Other metadata
    grid.dx = readDoubleAttribute(file, key_dx);
    grid.dy = readDoubleAttribute(file, key_dy);
    t = readDoubleAttribute(file, key_time);
    cycle = readIntAttribute(file, key_cyc);
    
    const int write_opt = readIntAttribute(file, key_wrt_opt);
    
    //Write-option-indepependent values
    std::vector<double> rho = readArray(file, key_rho);
    #ifdef MHD
    std::vector<double> Bx = mhd ? readArray(file, key_Bx) : std::vector<double>();
    std::vector<double> By = mhd ? readArray(file, key_By) : std::vector<double>();
    std::vector<double> Bz = mhd ? readArray(file, key_Bz) : std::vector<double>();
    std::vector<double> Az = mhd ? readArray(file, key_Az) : std::vector<double>();
    #endif
    for(int i = i0; i<in; i++){
        for(int j=j0; j<jn; j++){
            size_t n = (j-j0)*(in-i0) + (i-i0);
            grid[i,j].rho = rho[n];
            #ifdef MHD
            if(mhd) {
                grid[i,j].B = {Bx[n],By[n],Bz[n]};
            } else {
                grid[i,j].B = {0,0,0};
            }
            #endif
        }
    }
    #ifdef MHD
    for(int i = i0; i<=in; i++){
        for(int j=j0; j<=jn; j++){
            if(mhd) {
                size_t n = (j-j0)*(in+1-i0) + (i-i0);
                grid._A()[i,j] = {0,0,Az[n]};
            } else {
                grid._A()[i,j] = {0,0,0};
            }
        }
    }
    #endif
    
    if(write_opt & HDF5_WRITE_PRIMITIVE) { //Load from Primitive Variables
        std::vector<double> p = readArray(file, key_p);
        std::vector<double> vx = readArray(file, key_vx);
        std::vector<double> vy = readArray(file, key_vy);
        std::vector<double> vz = readArray(file, key_vz);
        for(int i = i0; i<in; i++){
            for(int j=j0; j<jn; j++){
                size_t n = (j-j0)*(in-i0) + (i-i0);
                grid[i,j].p = p[n];
                grid[i,j].v = {vx[n],vy[n],vz[n]};
            }
        }
    } else { //Load from Conservative Variables
        std::vector<double> E = readArray(file, key_E);
        std::vector<double> px = readArray(file, key_px);
        std::vector<double> py = readArray(file, key_py);
        std::vector<double> pz = readArray(file, key_pz);
        for(int i = i0; i<in; i++){
            for(int j=j0; j<jn; j++){
                size_t n = (j-j0)*(in-i0) + (i-i0);
                ConservativeState U;
                U.rho = grid[i,j].rho;
                U.p = {px[n],py[n],pz[n]};
                U.E = E[n];
                #ifdef MHD
                U.B = grid[i,j].B;
                #endif
                grid[i,j] = U;
            }
        }
    }
}
//MARK: 3D
void DRAGONHOARD::loadFromFile(Grid3D& grid, double& t, int& cycle, const std::string& filename){
    std::string path = CONFIG::output_dir + "/" + checkExtension(filename);
    H5::H5File file(path, H5F_ACC_RDONLY);
    
    //Verify compatibility
    const int fmt = readIntAttribute(file, key_fmt);
    if(fmt > 1) throw std::runtime_error("Please update to the latest version of DRAGON to read this file");
    const int dim = readIntAttribute(file, key_dim);
    if(dim != 3) throw std::runtime_error(std::format("Expected 3 dimensions but got {}",dim));
    const bool mhd = readIntAttribute(file, key_mhd);
    #ifndef MHD
    if(mhd) throw std::runtime_error("File was saved using MHD. Please enable MHD in Config.h");
    #endif
    
    //Extract size
    const int nx = readIntAttribute(file, key_nx);
    const int ny = readIntAttribute(file, key_ny);
    const int nz = readIntAttribute(file, key_nz);
    if(nx != grid.getSizeX() || ny != grid.getSizeY() || nz != grid.getSizeZ()) throw std::runtime_error("File size does not match the grid size");
    
    int ng = readIntAttribute(file, key_ng);
    if(ng > grid.getGhosts()) ng = grid.getGhosts();
    const int i0 = -ng, in = nx+ng, j0 = -ng, jn = ny+ng, k0 = -ng, kn = nz+ng;

    //Other metadata
    grid.dx = readDoubleAttribute(file, key_dx);
    grid.dy = readDoubleAttribute(file, key_dy);
    grid.dz = readDoubleAttribute(file, key_dz);
    t = readDoubleAttribute(file, key_time);
    cycle = readIntAttribute(file, key_cyc);
    
    const int write_opt = readIntAttribute(file, key_wrt_opt);

    //Write-option-indepependent values
    std::vector<double> rho = readArray(file, key_rho);
    #ifdef MHD
    std::vector<double> Bx = mhd ? readArray(file, key_Bx) : std::vector<double>();
    std::vector<double> By = mhd ? readArray(file, key_By) : std::vector<double>();
    std::vector<double> Bz = mhd ? readArray(file, key_Bz) : std::vector<double>();
    std::vector<double> Ax = mhd ? readArray(file, key_Ax) : std::vector<double>();
    std::vector<double> Ay = mhd ? readArray(file, key_Ay) : std::vector<double>();
    std::vector<double> Az = mhd ? readArray(file, key_Az) : std::vector<double>();
    #endif
    for(int i = i0; i<in; i++){
        for(int j=j0; j<jn; j++){
            for(int k=k0; k<kn; k++){
                size_t n = ((k-k0)*(jn-j0) + (j-j0))*(in-i0) + (i-i0);
                grid[i,j,k].rho = rho[n];
                #ifdef MHD
                if(mhd) {
                    grid[i,j,k].B = {Bx[n],By[n],Bz[n]};
                } else {
                    grid[i,j,k].B = {0,0,0};
                }
                #endif
            }
        }
    }
    #ifdef MHD
    for(int i = i0; i<=in; i++){
        for(int j=j0; j<=jn; j++){
            for(int k=k0; k<=kn; k++){
                if(mhd) {
                    size_t n = ((k-k0)*(jn+1-j0) + (j-j0))*(in+1-i0) + (i-i0);
                    grid._A()[i,j,k] = {Ax[n],Ay[n],Az[n]};
                } else {
                    grid._A()[i,j,k] = {0,0,0};
                }
            }
        }
    }
    #endif
    
    if(write_opt & HDF5_WRITE_PRIMITIVE) { //Load from Primitive Variables
        std::vector<double> p = readArray(file, key_p);
        std::vector<double> vx = readArray(file, key_vx);
        std::vector<double> vy = readArray(file, key_vy);
        std::vector<double> vz = readArray(file, key_vz);
        for(int i = i0; i<in; i++){
            for(int j=j0; j<jn; j++){
                for(int k=k0; k<kn; k++){
                    size_t n = ((k-k0)*(jn-j0) + (j-j0))*(in-i0) + (i-i0);
                    grid[i,j,k].p = p[n];
                    grid[i,j,k].v = {vx[n],vy[n],vz[n]};
                }
            }
        }
    } else { //Load from Conservative Variables
        std::vector<double> E = readArray(file, key_E);
        std::vector<double> px = readArray(file, key_px);
        std::vector<double> py = readArray(file, key_py);
        std::vector<double> pz = readArray(file, key_pz);
        for(int i = i0; i<in; i++){
            for(int j=j0; j<jn; j++){
                for(int k=k0; k<kn; k++){
                    size_t n = ((k-k0)*(jn-j0) + (j-j0))*(in-i0) + (i-i0);
                    ConservativeState U;
                    U.rho = grid[i,j,k].rho;
                    U.p = {px[n],py[n],pz[n]};
                    U.E = E[n];
                    #ifdef MHD
                    U.B = grid[i,j,k].B;
                    #endif
                    grid[i,j,k] = U;
                }
            }
        }
    }
}



