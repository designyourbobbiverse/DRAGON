//
//  DistGrid.cpp
//  DRAGON/Refinement
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#include "DistGrid.hpp"
#include "DragonWing.hpp"
#include "CFL.hpp"
#include <math.h>
#include <cassert>
#include <iostream>

//MARK: Bin Setup
//Compute the size of the ith children, given that we have nx cells spread across ncx bins
static int computeChildSize(int nx, int ncx, int i){
    int base_size = nx / ncx;
    int remainder = nx % ncx;
    return base_size + (i < remainder ? 1 : 0);
}

static int validGhosts(int g){ //How many ghost cells are needed to do this correctly
#if !defined(DIMENSION_UNSPLIT) //Split schemes need extra ghosts to avoid needing to sync between substeps
    #if defined(MUSCL_Hancock)
    return std::max(g, 4);
    #else
    return std::max(g, 3);
    #endif
#elif  defined(MHD) && defined (CTU)
    return std::max(g, 3);
#else
    return std::max(g, 2); //Unsplit syncs afeter every advance, so unless using MHD+CTU it's safe to use 2 ghosts
#endif
}

DistGrid1D::DistGrid1D(int nx, double dx_, int g, bool root): Grid1D(nx, dx_,validGhosts(g)), size_x(nx*dx_){
    
    
    ncx = root ? CONFIG::core_count : 1;

    if (ncx > 1){ //Children
        children.reserve(ncx);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
        
            auto child = std::make_unique<DistGrid1D>(_nx, dx, getGhosts(), false);
            child->boundary = Boundary::Ignore();
            children.push_back(std::move(child));
        }
    }
}

DistGrid2D::DistGrid2D(int nx, int ny, double dx_, double dy_, int g, bool root): Grid2D(nx, ny, dx_, dy_, validGhosts(g)), size_x(nx*dx_),size_y(ny*dy_) {

    if(!root){ ncx = 1; ncy = 1; return;}
    
    //Goal: nx/ncx = ny/ncy,  ncx*ncy = core_count
    ncx = ceil(sqrt( CONFIG::core_count * double(nx)/double(ny)) ); //ncx = core_count / ncy =  core_count * (nx/ny)/ ncx
    ncy = ceil(sqrt( CONFIG::core_count * double(ny)/double(nx)) ); //ncy = core_count / ncx =  core_count * (ny/nx)/ ncy

    
    if (ncx > 1 || ncy > 1){ //Children
        children.reserve(ncx * ncy);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
            for(int j = 0; j<ncy; j++){
                int _ny = computeChildSize(ny, ncy, j);

                auto child = std::make_unique<DistGrid2D>(_nx, _ny, dx, dy, getGhosts(), false);
                child->boundary = Boundary::Ignore();
                children.push_back(std::move(child));
            }
        }
    }
        
}

DistGrid3D::DistGrid3D(int nx, int ny, int nz, double dx_, double dy_, double dz_, int g, bool root):
    Grid3D(nx, ny, nz, dx_, dy_, dz_, validGhosts(g)), size_x(nx*dx_),size_y(ny*dy_), size_z(nz*dz_){

    if(!root){ ncx = 1; ncy = 1; ncz = 1; return;}
        
    //Goal: nx/ncx = ny/ncy = nz/ncz,  ncx*ncy*ncz = core_count
    double rxy = double(nx)/double(ny), rxz = double(nx)/double(nz), ryz = double(ny)/double(nz);
    ncx = ceil(pow( CONFIG::core_count * rxy*rxz, 0.3333) ); //ncx = core_count / (ncy*ncz) =  core_count * (nx/ny) * (nx/nz)/ ncx^2
    ncy = ceil(pow( CONFIG::core_count * ryz/rxy, 0.3333) ); //ncy = core_count / (ncx*ncz) =  core_count * (ny/nx) * (ny/nz)/ ncy^2
    ncz = ceil(pow( CONFIG::core_count / (rxz*ryz), 0.3333) ); //ncz = core_count / (ncy*ncz) =  core_count * (nz/nx) * (nz/ny)/ ncz^2

    if (ncx > 1 || ncy > 1 || ncz > 1){ //Children
        children.reserve(ncx * ncy * ncz);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
            for(int j = 0; j<ncy; j++){
                int _ny = computeChildSize(ny, ncy, j);
                for(int k = 0; k<ncz; k++){
                    int _nz = computeChildSize(nz, ncz, k);
                    
                    auto child = std::make_unique<DistGrid3D>(_nx, _ny, _nz, dx, dy, dz, getGhosts(), false);
                    child->boundary = Boundary::Ignore();
                    children.push_back(std::move(child));
                }
            }
        }
    }
        
}

//MARK: Parent -> Child

void DistGrid1D::pushToChildren(){
    if(ncx == 1) return; //This is the child
    const int ng = getGhosts();
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        std::unique_ptr<DistGrid1D>& child = children[zi];
        //Copy parent data to child
        for(int i = -ng; i < child->getSize() + ng; i++){
            (*child)[i] = w[i+x_offset];
            
        }
        child->boundary = Boundary::Ignore(); //Don't let the boundary API overwrite this
        child->pushToChildren(); //If children have children, make them sync too
           
        x_offset += child->getSize();
    }
}
void DistGrid2D::pushToChildren(){
    if(ncx*ncy == 1) return;
    const int ng = getGhosts();

    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        const int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<DistGrid2D>& child = children[ncy*zi + zj];
            const int _ny = child->getSizeY();
            //Copy parent fluid to child
            for(int i = -ng; i < _nx + ng; i++){
                for(int j = -ng; j < _ny +ng; j++){
                    (*child)[i,j] = w[i+x_offset, j+y_offset];
                }
            }
            #ifdef MHD//Copy parent magnetic potential to child
            for(int i = -ng; i <= _nx + ng; i++){
                for(int j = -ng; j <= _ny +ng; j++){
                    child->A[i,j] = A[i+x_offset, j+y_offset];
                }
            }
            #endif
            child->boundary = Boundary::Ignore();//Don't let the boundary API overwrite this
            child->pushToChildren();//If children have children, make them sync too
            
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}
void DistGrid3D::pushToChildren(){
    if(ncx*ncy*ncz == 1) return;
    const int ng = getGhosts();
    
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        const int _nx = children[ncz*ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            const int _ny = children[ncz*(ncy*zi + zj)]->getSizeY();
            int z_offset = 0;
            for(int zk = 0; zk<ncz; zk++){
                std::unique_ptr<DistGrid3D>& child = children[ncz*(ncy*zi + zj) + zk];
                const int _nz = child->getSizeZ();
                //Copy parent fluid to child
                for(int i = -ng; i < _nx + ng; i++){
                    for(int j = -ng; j < _ny + ng; j++){
                        for(int k = -ng; k < _nz + ng; k++){
                            (*child)[i,j,k] = w[i+x_offset, j+y_offset, k+z_offset];
                        }
                    }
                }
                #ifdef MHD//Copy parent magnetic potential to child
                for(int i = -ng; i <= _nx + ng; i++){
                    for(int j = -ng; j <= _ny + ng; j++){
                        for(int k = -ng; k <= _nz + ng; k++){
                            child->A[i,j,k] = A[i+x_offset, j+y_offset, k+z_offset];
                        }
                    }
                }
                #endif
                child->boundary = Boundary::Ignore();//Don't let the boundary API overwrite this
                child->pushToChildren();//If children have children, make them sync too
                
                z_offset += _nz;
            }
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}


//MARK: Child -> Parent
void DistGrid1D::loadFromChildren(){
    if(ncx == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        std::unique_ptr<DistGrid1D>& child = children[zi];
        child->loadFromChildren(); //If children have children, make them sync first
        //Copy child data to parent
        for(int i = 0; i < child->getSize(); i++) {
            w[i+x_offset] = (*child)[i];
        }
        x_offset += child->getSize();
    }
}

void DistGrid2D::loadFromChildren(){
    if(ncx*ncy == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<DistGrid2D>& child = children[ncy*zi + zj];
            child->loadFromChildren();//If children have children, make them sync first
            int _ny = child->getSizeY();
            //Copy child fluid to parent
            for(int i = 0; i < _nx; i++){
                for(int j=0; j < _ny; j++){
                    w[i+x_offset, j+y_offset] = (*child)[i,j];
                }
            }
            #ifdef MHD//Copy child magnetic potential to parent
            for(int i = 0; i <= _nx; i++){
                for(int j=0; j <= _ny; j++){
                    A[i+x_offset, j+y_offset] = child->A[i,j];
                }
            }
            #endif
            
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}

void DistGrid3D::loadFromChildren(){
    if(ncx*ncy*ncz == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncz*ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            int _ny = children[ncz*(ncy*zi+zj)]->getSizeY();
            int z_offset = 0;
            for(int zk = 0; zk<ncz; zk++){
                std::unique_ptr<DistGrid3D>& child = children[ncz*(ncy*zi + zj) + zk];
                child->loadFromChildren();//If children have children, make them sync first
                int _nz = child->getSizeZ();
                //Copy child fluid to parent
                for(int i = 0; i < _nx; i++){
                    for(int j=0; j < _ny; j++){
                        for(int k=0; k < _nz; k++){
                            w[i+x_offset, j+y_offset, k+z_offset] = (*child)[i,j,k];
                        }
                    }
                }
                #ifdef MHD //Copy child magnetic potential to parent
                for(int i = 0; i <= _nx; i++){
                    for(int j=0; j <= _ny; j++){
                        for(int k=0; k <= _nz; k++){
                            _A()[i+x_offset, j+y_offset, k+z_offset] = child->_A()[i,j,k];
                        }
                    }
                }
                #endif

                z_offset += _nz;
            }
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}

//MARK: Advance
template <typename T> void DistGrid<T>::step(double dt){
    pushToChildren();
    DRAGONWING::initialize(static_cast<int>(children.size()));
    for(auto& child : children){
        DRAGONWING::launchParallel(child.get(), dt);
    }
    bool success = DRAGONWING::waitForCheckpoint2(); //Wait for children to finish
    if(!success) { //If we failed, try again with half time step
        throw std::runtime_error(DRAGONWING::restartMsg());
    }
    //Copy Back
    loadFromChildren();
}


//MARK: Subclass Dispatch
void DistGrid1D::split_step(double dt){
    if(children.size() <= 1 ){
        Grid1D::split_step(dt);
        DRAGONWING::reportCheckpoint2();
    } else { DistGrid::step(dt); }
}
void DistGrid1D::unsplit_step(double dt){
    if(children.size() <= 1 ){
        Grid1D::unsplit_step(dt);
        DRAGONWING::reportCheckpoint2();
        return;
    } else { DistGrid::step(dt); }
}
void DistGrid2D::split_step(double dt){
    if(children.size() <= 1 ){
        Grid2D::split_step(dt);
        DRAGONWING::reportCheckpoint2();
    } else { DistGrid::step(dt); }
}
void DistGrid2D::unsplit_step(double dt){
    if(children.size() <= 1 ){
        Grid2D::unsplit_step(dt);
        DRAGONWING::reportCheckpoint2();
        return;
    } else { DistGrid::step(dt); }
}

void DistGrid3D::split_step(double dt){
    if(children.size() <= 1 ){
        Grid3D::split_step(dt);
        DRAGONWING::reportCheckpoint2();
    } else {  DistGrid::step(dt); }
}
void DistGrid3D::unsplit_step(double dt){
    if(children.size() <= 1 ){
        Grid3D::unsplit_step(dt);
        DRAGONWING::reportCheckpoint2();
    } else { DistGrid::step(dt); }
}


bool DistGrid1D::on_step_fail(const std::exception& e){
    if(children.size() <= 1 ) return Grid1D::on_step_fail(e);
    return true;
}
bool DistGrid2D::on_step_fail(const std::exception& e){
    if(children.size() <= 1 ) return Grid2D::on_step_fail(e);
    return true;
}
bool DistGrid3D::on_step_fail(const std::exception& e){
    if(children.size() <= 1 ) return Grid3D::on_step_fail(e);
    return true;
}
