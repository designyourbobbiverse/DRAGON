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
#if defined(MUSCL_Hancock) && !defined(DIMENSION_UNSPLIT)
    return std::max(g, 3); //MUSCL Split needs an extra ghost to avoid needing to sync after each substep
#else
    return std::max(g, 2); //Unsplit syncs afeter every advance, so it's safe to use the usual 2
#endif
}

DistGrid1D::DistGrid1D(int nx, double dx_, int g, bool root): size_x(nx*dx_),dx(dx_), ghosts(validGhosts(g)), data(nx, dx_,validGhosts(g)) {
    
    
    ncx = root ? core_count : 1;

    if (ncx > 1){ //Children
        children.reserve(ncx);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
        
            auto child = std::make_unique<DistGrid1D>(_nx, dx, ghosts, false);
            child->data.boundary = Boundary::Ignore();
            children.push_back(std::move(child));
        }
    }
}

DistGrid2D::DistGrid2D(int nx, int ny, double dx_, double dy_, int g, bool root): size_x(nx*dx_),size_y(ny*dy_), dx(dx_), dy(dy_), ghosts(validGhosts(g)), data(nx, ny, dx_, dy_, validGhosts(g)) {

    if(!root){ ncx = 1; ncy = 1; return;}
    
    //Goal: nx/ncx = ny/ncy,  ncx*ncy = core_count
    ncx = ceil(sqrt( core_count * double(nx)/double(ny)) ); //ncx = core_count / ncy =  core_count * (nx/ny)/ ncx
    ncy = ceil(sqrt( core_count * double(ny)/double(nx)) ); //ncy = core_count / ncx =  core_count * (ny/nx)/ ncy

    
    if (ncx > 1 || ncy > 1){ //Children
        children.reserve(ncx * ncy);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
            for(int j = 0; j<ncy; j++){
                int _ny = computeChildSize(ny, ncy, j);

                auto child = std::make_unique<DistGrid2D>(_nx, _ny, dx, dy, ghosts, false);
                child->data.boundary = Boundary::Ignore();
                children.push_back(std::move(child));
            }
        }
    }
        
}

DistGrid3D::DistGrid3D(int nx, int ny, int nz, double dx_, double dy_, double dz_, int g, bool root):
    size_x(nx*dx_),size_y(ny*dy_), size_z(nz*dz_), dx(dx_), dy(dy_), dz(dz_), ghosts(validGhosts(g)),
    data(nx, ny, nz, dx_, dy_, dz_, validGhosts(g)) {

    if(!root){ ncx = 1; ncy = 1; ncz = 1; return;}
    //Goal: nx/ncx = ny/ncy = nz/ncz,  ncx*ncy*ncz = core_count
    double rxy = double(nx)/double(ny), rxz = double(nx)/double(nz), ryz = double(ny)/double(nz);
    ncx = ceil(pow( core_count * rxy*rxz, 0.3333) ); //ncx = core_count / (ncy*ncz) =  core_count * (nx/ny) * (nx/nz)/ ncx^2
    ncy = ceil(pow( core_count * ryz/rxy, 0.3333) ); //ncy = core_count / (ncx*ncz) =  core_count * (ny/nx) * (ny/nz)/ ncy^2
    ncz = ceil(pow( core_count / (rxz*ryz), 0.3333) ); //ncz = core_count / (ncy*ncz) =  core_count * (nz/nx) * (nz/ny)/ ncz^2

    if (ncx > 1 || ncy > 1 || ncz > 1){ //Children
        children.reserve(ncx * ncy * ncz);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
            for(int j = 0; j<ncy; j++){
                int _ny = computeChildSize(ny, ncy, j);
                for(int k = 0; k<ncz; k++){
                    int _nz = computeChildSize(nz, ncz, k);
                    
                    auto child = std::make_unique<DistGrid3D>(_nx, _ny, _nz, dx, dy, dz, ghosts, false);
                    child->data.boundary = Boundary::Ignore();
                    children.push_back(std::move(child));
                }
            }
        }
    }
        
}


//MARK: Grid Access
PrimitiveState& DistGrid1D::operator[](int i){ return data[i];}
const PrimitiveState& DistGrid1D::operator[](int i) const { return data[i]; }

PrimitiveState& DistGrid2D::operator[](int i,int j){ return data[i,j];}
const PrimitiveState& DistGrid2D::operator[](int i,int j) const { return data[i,j]; }

PrimitiveState& DistGrid3D::operator[](int i,int j,int k){ return data[i,j,k];}
const PrimitiveState& DistGrid3D::operator[](int i,int j,int k) const { return data[i,j,k]; }

int DistGrid1D::getSize() const{ return data.getSize(); }
int DistGrid1D::getGhosts() const{ return data.getGhosts(); }

int DistGrid2D::getSizeX() const{ return data.getSizeX(); }
int DistGrid2D::getSizeY() const{ return data.getSizeY(); }
int DistGrid2D::getGhosts() const{ return data.getGhosts(); }

int DistGrid3D::getSizeX() const{ return data.getSizeX(); }
int DistGrid3D::getSizeY() const{ return data.getSizeY(); }
int DistGrid3D::getSizeZ() const{ return data.getSizeZ(); }
int DistGrid3D::getGhosts() const{ return data.getGhosts(); }


//MARK: Parent -> Child

void DistGrid1D::pushToChildren(){
    if(ncx == 1) return; //This is the child
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        std::unique_ptr<DistGrid1D>& child = children[zi];
        //Copy parent data to child
        for(int i = -ghosts; i < child->getSize() + ghosts; i++){
            (*child)[i] = data[i+x_offset];
            
        }
        child->boundary = Boundary::Ignore(); //Don't let the boundary API overwrite this
        child->data.boundary = Boundary::Ignore();
        child->pushToChildren(); //If children have children, make them sync too
           
        x_offset += child->getSize();
    }
}
void DistGrid2D::pushToChildren(){
    if(ncx*ncy == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<DistGrid2D>& child = children[ncy*zi + zj];
            int _ny = child->getSizeY();
            //Copy parent fluid to child
            for(int i = -ghosts; i < _nx + ghosts; i++){
                for(int j = -ghosts; j < _ny +ghosts; j++){
                    (*child)[i,j] = data[i+x_offset, j+y_offset];
                }
            }
            #ifdef MHD//Copy parent magnetic potential to child
            for(int i = -ghosts; i <= _nx + ghosts; i++){
                for(int j = -ghosts; j <= _ny +ghosts; j++){
                    child->A()[i,j] = A()[i+x_offset, j+y_offset];
                }
            }
            #endif
            child->boundary = Boundary::Ignore();//Don't let the boundary API overwrite this
            child->data.boundary = Boundary::Ignore();
            child->pushToChildren();//If children have children, make them sync too
            
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}
void DistGrid3D::pushToChildren(){
    if(ncx*ncy*ncz == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncz*ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            int _ny = children[ncz*(ncy*zi + zj)]->getSizeY();
            int z_offset = 0;
            for(int zk = 0; zk<ncz; zk++){
                std::unique_ptr<DistGrid3D>& child = children[ncz*(ncy*zi + zj) + zk];
                int _nz = child->getSizeZ();
                //Copy parent fluid to child
                for(int i = -ghosts; i < _nx + ghosts; i++){
                    for(int j = -ghosts; j < _ny + ghosts; j++){
                        for(int k = -ghosts; k < _nz + ghosts; k++){
                            (*child)[i,j,k] = data[i+x_offset, j+y_offset, k+z_offset];
                        }
                    }
                }
                #ifdef MHD//Copy parent magnetic potential to child
                for(int i = -ghosts; i <= _nx + ghosts; i++){
                    for(int j = -ghosts; j <= _ny + ghosts; j++){
                        for(int k = -ghosts; k <= _nz + ghosts; k++){
                            child->A()[i,j,k] = A()[i+x_offset, j+y_offset, k+z_offset];
                        }
                    }
                }
                #endif
                child->boundary = Boundary::Ignore();//Don't let the boundary API overwrite this
                child->data.boundary = Boundary::Ignore();
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
            data[i+x_offset] = (*child)[i];
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
                    data[i+x_offset, j+y_offset] = (*child)[i,j];
                }
            }
            #ifdef MHD//Copy child magnetic potential to parent
            for(int i = 0; i <= _nx; i++){
                for(int j=0; j <= _ny; j++){
                    A()[i+x_offset, j+y_offset] = child->A()[i,j];
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
                            data[i+x_offset, j+y_offset, k+z_offset] = (*child)[i,j,k];
                            #ifdef MHD
                            A()[i+x_offset, j+y_offset, k+z_offset] = child->A()[i,j,k];
                            #endif
                        }
                    }
                }
                #ifdef MHD //Copy child magnetic potential to parent
                for(int i = 0; i <= _nx; i++){
                    for(int j=0; j <= _ny; j++){
                        for(int k=0; k <= _nz; k++){
                            A()[i+x_offset, j+y_offset, k+z_offset] = child->A()[i,j,k];
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
void DistGrid1D::advance(double dt, bool check_cfl){
    data.boundary = std::move(boundary);
    if(ncx == 1 ){
        data.advance(dt, check_cfl);
        DRARGONWING::reportCheckpoint2();
        return;
    }
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(data);
        pushToChildren();
        
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(data)) : dt;
        //Execute the Advancement
        bool success = false;
        do{
            DRARGONWING::initialize(ncx);
            for(auto& child : children){
                DRARGONWING::launchParallel(child.get(), t1);
            }
            success = DRARGONWING::waitForCheckpoint2(); //Wait for children to finish
            if(!success) t1 /= 2; //If we failed, try again with half time step
        } while (!success);
        dt -= t1;
        //Copy Back
        loadFromChildren();
    }
}


void DistGrid2D::advance(double dt, bool check_cfl){
    data.boundary = std::move(boundary);
    if(ncx*ncy == 1 ){
        data.advance(dt, check_cfl);
        DRARGONWING::reportCheckpoint2();
        boundary = std::move(data.boundary);
        return;
    }
    data.initialize_B_fields();
    boundary = std::move(data.boundary);
    
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(data);
        data.initialize_B_fields();
        pushToChildren();
        
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(data)) : dt;
        //Execute the Advancement
        bool success = false;
        do{
            DRARGONWING::initialize(ncx*ncy);
            for(auto& child : children){
                DRARGONWING::launchParallel(child.get(), t1);
            }
            success = DRARGONWING::waitForCheckpoint2(); //Wait for children to finish
            if(!success) t1 /= 2; //If we failed, try again with half time step
        } while (!success);
        dt -= t1;
        //Copy Back
        loadFromChildren();
    }
    boundary = std::move(data.boundary);
}


void DistGrid3D::advance(double dt, bool check_cfl){
    data.boundary = std::move(boundary);
    if(ncx*ncy*ncz == 1 ){
        data.advance(dt, check_cfl);
        DRARGONWING::reportCheckpoint2();
        boundary = std::move(data.boundary);
        return;
    }
    data.initialize_B_fields();
    boundary = std::move(data.boundary);

    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(data);
        pushToChildren();
        
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(data)) : dt;
        //Execute the Advancement
        bool success = false;
        do{
            DRARGONWING::initialize(ncx*ncy*ncz);
            for(auto& child : children){
                DRARGONWING::launchParallel(child.get(), t1);
            }
            success = DRARGONWING::waitForCheckpoint2(); //Wait for children to finish
            if(!success) t1 /= 2; //If we failed, try again with half time step
        } while (!success);
        dt -= t1;
        //Copy Back
        loadFromChildren();
    }
}
