//
//  HDF5Output.cpp
//  DRAGON/IO
//
//  Created by Bobbie Markwick on 03/07/2026.
//

#include "HDF5Output.hpp"
#include "Config.h"
#include "HDF5_Attrs.hpp"
#include <H5Cpp.h>
#include <vector>
#include <format>
#include <iostream>
#include <stdexcept>

#define HDF5_WRITE_ENERGY (HDF5_WRITE_PRIMITIVE_AND_ENERGY - HDF5_WRITE_PRIMITIVE)
#define HDF5_WRITE_PRIMS (HDF5_WRITE_OPTION & HDF5_WRITE_PRIMITIVE)
#define HDF5_WRITE_CONSV (HDF5_WRITE_OPTION & (HDF5_WRITE_CONSERVATIVE - HDF5_WRITE_ENERGY))
#define HDF5_WRITE_E (HDF5_WRITE_OPTION & HDF5_WRITE_ENERGY)


#if !(HDF5_WRITE_PRIMS || HDF5_WRITE_CONSV)
#error Specified write option does not write a sufficient description of the fluid. Please reconfigure HDF5_WRITE_OPTION in Config.h
#endif


std::string IO::cycle_string(int n){
    std::string zeroes = n<10 ? "0000" : (n<100 ? "000" : (n<1000 ? "00" : (n<10000 ? "0" : "")));
    return zeroes +  std::to_string(n);
}

//MARK: Helpers
namespace{

void writeArray(H5::H5File& file,const std::string& name, const std::vector<double>& data, hsize_t nx){
    hsize_t dims[1] = {nx};
    H5::DataSpace dataspace(1, dims);
    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}
void writeArray(H5::H5File& file,const std::string& name, const std::vector<double>& data, hsize_t nx,hsize_t ny){
    hsize_t dims[2] = {ny, nx};
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}
void writeArray(H5::H5File& file,const std::string& name, const std::vector<double>& data, hsize_t nx,hsize_t ny,hsize_t nz){
    hsize_t dims[3] = {nz, ny, nx};
    H5::DataSpace dataspace(3, dims);
    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}
void writeDoubleAttribute(H5::H5File& file, const std::string& name, double value) {
    H5::DataSpace dataspace(H5S_SCALAR);
    H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    attr.write(H5::PredType::NATIVE_DOUBLE, &value);
}

void writeIntAttribute(H5::H5File& file, const std::string& name, int value) {
    H5::DataSpace dataspace(H5S_SCALAR);
    H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_INT, dataspace);
    attr.write(H5::PredType::NATIVE_INT, &value);
}
}


//MARK: Dispatch
void IO::writeToFile(Grid& grid, double t, int cycle, const std::string& filename){
    Grid3D* grid3D = dynamic_cast<Grid3D*>(&grid);
    if(grid3D){
        writeToFile(*grid3D, t, cycle, filename);
        return;
    }
    Grid2D* grid2D = dynamic_cast<Grid2D*>(&grid);
    if(grid2D){
        writeToFile(*grid2D, t, cycle, filename);
        return;
    }
    Grid1D* grid1D = dynamic_cast<Grid1D*>(&grid);
    if(grid1D){
        writeToFile(*grid1D, t, cycle, filename);
        return;
    }
    
    throw std::runtime_error("IO attempted with unsupported Grid type");
}

//MARK: Writing - 1D

void IO::writeToFile(Grid1D& grid, double t, int cycle, const std::string& filename){
    const int nx = grid.getSize(), ng = grid.getGhosts();
    #ifdef WRITE_GHOSTS_TO_FILE
    const int i0 = -ng, in = nx+ng;
    #else
    const int i0 = 0, in = nx;
    #endif
 
    std::string path = CONFIG::output_dir + "/" + filename + file_ext;
    H5::H5File file(path, H5F_ACC_TRUNC);
    
    //Metadata
    writeIntAttribute(file, key_fmt, 1);
    writeIntAttribute(file, key_wrt_opt, HDF5_WRITE_OPTION);
    writeIntAttribute(file, key_dim, 1);
    #ifdef MHD
    writeIntAttribute(file, key_mhd, 1);
    #else
    writeIntAttribute(file, key_mhd, 0);
    #endif
    
    #ifdef WRITE_GHOSTS_TO_FILE
    writeIntAttribute(file, key_ng, ng);
    #else
    writeIntAttribute(file, key_ng, 0);
    #endif
    writeIntAttribute(file, key_nx, nx);
    writeDoubleAttribute(file, key_dx, grid.dx);
    writeIntAttribute(file, key_cyc, cycle);
    writeDoubleAttribute(file, key_time, t);
    
    //Flatten Grid
    size_t size = (in-i0);
    std::vector<double> rho(size);
    #if HDF5_WRITE_PRIMS
    std::vector<double> vx(size);
    std::vector<double> vy(size);
    std::vector<double> vz(size);
    std::vector<double> p(size);
    #endif
    #if HDF5_WRITE_E
    std::vector<double> E(size);
    #if HDF5_WRITE_CONSV
    std::vector<double> px(size);
    std::vector<double> py(size);
    std::vector<double> pz(size);
    #endif
    #endif
    #ifdef MHD
    std::vector<double> Bx(size);
    std::vector<double> By(size);
    std::vector<double> Bz(size);
    #endif
    for(int i = i0; i<in; i++){
        PrimitiveState w = grid[i];
        size_t n = i-i0;

        rho[n] = w.rho;
        #if HDF5_WRITE_PRIMS
        vx[n] = w.v.x;
        vy[n] = w.v.y;
        vz[n] = w.v.z;
        p[n] = w.p;
        #endif
        #if HDF5_WRITE_CONSV
        ConservativeState u(w);
        px[n] = u.p.x;
        py[n] = u.p.y;
        pz[n] = u.p.z;
        E[n] = u.E;
        #elif HDF5_WRITE_E
        E[n] = w.energy();
        #endif
        #ifdef MHD
        Bx[n] = w.B.x;
        By[n] = w.B.y;
        Bz[n] = w.B.z;
        #endif
    }
    
    //Write Grid
    file.createGroup(key_fluid);
    writeArray(file, key_rho, rho, in-i0);
    #if HDF5_WRITE_PRIMS
    file.createGroup(key_v);
    writeArray(file, key_vx, vx, in-i0);
    writeArray(file, key_vy, vy, in-i0);
    writeArray(file, key_vz, vz, in-i0);
    writeArray(file, key_p, p, in-i0);
    #endif
    #if HDF5_WRITE_E
    writeArray(file, key_E, E, in-i0);
    #if HDF5_WRITE_CONSV
    file.createGroup(key_mom);
    writeArray(file, key_px, px, in-i0);
    writeArray(file, key_py, py, in-i0);
    writeArray(file, key_pz, pz, in-i0);
    #endif
    #endif
    #ifdef MHD
    file.createGroup(key_B);
    writeArray(file, key_Bx, Bx, in-i0);
    writeArray(file, key_By, By, in-i0);
    writeArray(file, key_Bz, Bz, in-i0);
    #endif

}

//MARK: Writing - 2D
void IO::writeToFile(Grid2D& grid, double t, int cycle, const std::string& filename){
    const int nx = grid.getSizeX(), ny = grid.getSizeY(), ng = grid.getGhosts();
    #ifdef WRITE_GHOSTS_TO_FILE
    const int i0 = -ng, in = nx+ng, j0 = -ng, jn = ny+ng;
    #else
    const int i0 = 0, in = nx, j0 = 0, jn = ny;
    #endif

    std::string path = CONFIG::output_dir + "/" + filename + file_ext;
    H5::H5File file(path, H5F_ACC_TRUNC);
    
    //Metadata
    writeIntAttribute(file, key_fmt, 1);
    writeIntAttribute(file, key_wrt_opt, HDF5_WRITE_OPTION);
    writeIntAttribute(file, key_dim, 2);
    #ifdef MHD
    writeIntAttribute(file, key_mhd, 1);
    #else
    writeIntAttribute(file, key_mhd, 0);
    #endif
    
    #ifdef WRITE_GHOSTS_TO_FILE
    writeIntAttribute(file, key_ng, ng);
    #else
    writeIntAttribute(file, key_ng, 0);
    #endif
    writeIntAttribute(file, key_nx, nx);
    writeIntAttribute(file, key_ny, ny);
    writeDoubleAttribute(file, key_dx, grid.dx);
    writeDoubleAttribute(file, key_dy, grid.dy);
    writeIntAttribute(file, key_cyc, cycle);
    writeDoubleAttribute(file, key_time, t);
    
    //Flatten Grid
    size_t size = (in-i0)*(jn-j0);
    std::vector<double> rho(size);
    #if HDF5_WRITE_PRIMS
    std::vector<double> vx(size);
    std::vector<double> vy(size);
    std::vector<double> vz(size);
    std::vector<double> p(size);
    #endif
    #if HDF5_WRITE_E
    std::vector<double> E(size);
    #if HDF5_WRITE_CONSV
    std::vector<double> px(size);
    std::vector<double> py(size);
    std::vector<double> pz(size);
    #endif
    #endif
    #ifdef MHD
    size_t size_A = (in+1-i0)*(jn+1-j0);
    std::vector<double> Bx(size);
    std::vector<double> By(size);
    std::vector<double> Bz(size);
    std::vector<double> Az(size_A);
    #endif
    for(int i = i0; i<in; i++){
        for(int j=j0; j<jn; j++){
            PrimitiveState w = grid[i,j];
            size_t n = (j-j0)*(in-i0) + (i-i0);
            
            rho[n] = w.rho;
            #if HDF5_WRITE_PRIMS
            vx[n] = w.v.x;
            vy[n] = w.v.y;
            vz[n] = w.v.z;
            p[n] = w.p;
            #endif
            #if HDF5_WRITE_CONSV
            ConservativeState u(w);
            px[n] = u.p.x;
            py[n] = u.p.y;
            pz[n] = u.p.z;
            E[n] = u.E;
            #elif HDF5_WRITE_E
            E[n] = w.energy();
            #endif
            #ifdef MHD
            Bx[n] = w.B.x;
            By[n] = w.B.y;
            Bz[n] = w.B.z;
            #endif
        }
    }
    #ifdef MHD
    for(int i = i0; i<=in; i++){
        for(int j=j0; j<=jn; j++){
            size_t n =  (j-j0)*(in-i0+1) + (i-i0);
            Az[n] = grid._A()[i,j].z;
        }
    }
    #endif
    
    //Write Grid data
    file.createGroup(key_fluid);
    writeArray(file, key_rho, rho, in-i0, jn-j0);
    #if HDF5_WRITE_PRIMS
    file.createGroup(key_v);
    writeArray(file, key_vx, vx, in-i0, jn-j0);
    writeArray(file, key_vy, vy, in-i0, jn-j0);
    writeArray(file, key_vz, vz, in-i0, jn-j0);
    writeArray(file, key_p, p, in-i0, jn-j0);
    #endif
    #if HDF5_WRITE_E
    writeArray(file, key_E, E, in-i0, jn-j0);
    #if HDF5_WRITE_CONSV
    file.createGroup(key_mom);
    writeArray(file, key_px, px, in-i0, jn-j0);
    writeArray(file, key_py, py, in-i0, jn-j0);
    writeArray(file, key_pz, pz, in-i0, jn-j0);
    #endif
    #endif
    #ifdef MHD
    file.createGroup(key_B);
    writeArray(file, key_Bx, Bx, in-i0, jn-j0);
    writeArray(file, key_By, By, in-i0, jn-j0);
    writeArray(file, key_Bz, Bz, in-i0, jn-j0);
    file.createGroup(key_A);
    writeArray(file, key_Az, Az, in+1-i0, jn+1-j0);
    #endif

}

//MARK: Writing - 3D
void IO::writeToFile(Grid3D& grid, double t, int cycle, const std::string& filename){
    const int nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ(), ng = grid.getGhosts();
    #ifdef WRITE_GHOSTS_TO_FILE
    const int i0 = -ng, in = nx+ng, j0 = -ng, jn = ny+ng, k0 = -ng, kn = nz+ng;
    #else
    const int i0 = 0, in = nx, j0 = 0, jn = ny, k0 = 0, kn = nz;
    #endif

    std::string path = CONFIG::output_dir + "/" + filename + file_ext;
    H5::H5File file(path, H5F_ACC_TRUNC);
    
    //Metadata
    writeIntAttribute(file, key_fmt, 1);
    writeIntAttribute(file, key_wrt_opt, HDF5_WRITE_OPTION);
    writeIntAttribute(file, key_dim, 3);
    #ifdef MHD
    writeIntAttribute(file, key_mhd, 1);
    #else
    writeIntAttribute(file, key_mhd, 0);
    #endif
    #ifdef WRITE_GHOSTS_TO_FILE
    writeIntAttribute(file, key_ng, ng);
    #else
    writeIntAttribute(file, key_ng, 0);
    #endif
    writeIntAttribute(file, key_nx, nx);
    writeIntAttribute(file, key_ny, ny);
    writeIntAttribute(file, key_nz, nz);
    writeDoubleAttribute(file, key_dx, grid.dx);
    writeDoubleAttribute(file, key_dy, grid.dy);
    writeDoubleAttribute(file, key_dz, grid.dz);
    writeIntAttribute(file, key_cyc, cycle);
    writeDoubleAttribute(file, key_time, t);

    
    //Flatten Grid
    size_t size = (in-i0)*(jn-j0)*(kn-k0);
    std::vector<double> rho(size);
    #if HDF5_WRITE_PRIMS
    std::vector<double> vx(size);
    std::vector<double> vy(size);
    std::vector<double> vz(size);
    std::vector<double> p(size);
    #endif
    #if HDF5_WRITE_E
    std::vector<double> E(size);
    #if HDF5_WRITE_CONSV
    std::vector<double> px(size);
    std::vector<double> py(size);
    std::vector<double> pz(size);
    #endif
    #endif
    #ifdef MHD
    size_t size_A = (in+1-i0)*(jn+1-j0)*(kn+1-k0);
    std::vector<double> Bx(size);
    std::vector<double> By(size);
    std::vector<double> Bz(size);
    std::vector<double> Ax(size_A);
    std::vector<double> Ay(size_A);
    std::vector<double> Az(size_A);
    #endif
    for(int i = i0; i<in; i++){
        for(int j=j0; j<jn; j++){
            for(int k=k0; k<kn; k++){
                PrimitiveState w = grid[i,j,k];
                size_t n = ((k-k0)*(jn-j0) + (j-j0))*(in-i0) + (i-i0);
                
                rho[n] = w.rho;
                #if HDF5_WRITE_PRIMS
                vx[n] = w.v.x;
                vy[n] = w.v.y;
                vz[n] = w.v.z;
                p[n] = w.p;
                #endif
                #if HDF5_WRITE_CONSV
                ConservativeState u(w);
                px[n] = u.p.x;
                py[n] = u.p.y;
                pz[n] = u.p.z;
                E[n] = u.E;
                #elif HDF5_WRITE_E
                E[n] = w.energy();
                #endif
                #ifdef MHD
                Bx[n] = w.B.x;
                By[n] = w.B.y;
                Bz[n] = w.B.z;
                #endif
            }
        }
    }
    #ifdef MHD
    for(int i = i0; i<=in; i++){
        for(int j=j0; j<=jn; j++){
            for(int k=k0; k<=kn; k++){
                vec3 A = grid._A()[i,j,k];
                size_t n = ((k-k0)*(jn-j0+1) + j-j0)*(in-i0+1) + (i-i0);
                Ax[n] = A.x;
                Ay[n] = A.y;
                Az[n] = A.z;
            }
        }
    }
    #endif
    
    //Write Grid
    file.createGroup(key_fluid);
    writeArray(file, key_rho, rho, in-i0, jn-j0, kn-k0);
    #if HDF5_WRITE_PRIMS
    file.createGroup(key_v);
    writeArray(file, key_vx, vx, in-i0, jn-j0, kn-k0);
    writeArray(file, key_vy, vy, in-i0, jn-j0, kn-k0);
    writeArray(file, key_vz, vz, in-i0, jn-j0, kn-k0);
    writeArray(file, key_p, p, in-i0, jn-j0, kn-k0);
    #endif
    #if HDF5_WRITE_E
    writeArray(file, key_E, E, in-i0, jn-j0, kn-k0);
    #if HDF5_WRITE_CONSV
    file.createGroup(key_mom);
    writeArray(file, key_px, px, in-i0, jn-j0, kn-k0);
    writeArray(file, key_py, py, in-i0, jn-j0, kn-k0);
    writeArray(file, key_pz, pz, in-i0, jn-j0, kn-k0);
    #endif
    #endif
    #ifdef MHD
    file.createGroup(key_B);
    writeArray(file, key_Bx, Bx, in-i0, jn-j0, kn-k0);
    writeArray(file, key_By, By, in-i0, jn-j0, kn-k0);
    writeArray(file, key_Bz, Bz, in-i0, jn-j0, kn-k0);
    file.createGroup(key_A);
    writeArray(file, key_Ax, Ax, in+1-i0, jn+1-j0, kn+1-k0);
    writeArray(file, key_Ay, Ay, in+1-i0, jn+1-j0, kn+1-k0);
    writeArray(file, key_Az, Az, in+1-i0, jn+1-j0, kn+1-k0);
    #endif
}



