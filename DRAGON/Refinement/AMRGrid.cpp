//
//  AMRGrid.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#include "AMRGrid.hpp"
#include "DragonWing.hpp"
#include "CFL.hpp"
#include <math.h>
#include <cassert>

//MARK: Bin Setup
//Compute the size of the ith children, given that we have nx cells spread across ncx bins
int computeChildSize(int nx, int ncx, int i){
    if(nx <= bin_size) return nx; //Single Bin
    int _nx = bin_size;
    if (i == 0 || i+1 == ncx) { //Reduce bin size if leftmost or rightmost
        _nx -=  (ncx*bin_size - nx)/2;
        if (i==0 && nx%2==1) _nx -= 1;
    }
    return _nx;
}

static int validGhosts(int g){ //How many ghost cells are needed to do this correctly
#if defined(MUSCL_Hancock) && !defined(DIMENSION_UNSPLIT)
    return std::max(g, 3); //MUSCL Split needs an extra ghost to avoid needing to sync after each substep
#else
    return std::max(g, 2); //Unsplit syncs afeter every advance, so it's safe to use the usual 2
#endif
}

AMRGrid1D::AMRGrid1D(int nx, double dx_, int g): size_x(nx*dx_),dx(dx_), ghosts(validGhosts(g)), data(nx, dx_,validGhosts(g)) {
    ncx = (nx / bin_size) + (nx % bin_size == 0 ? 0 : 1);
    
    if (ncx > 1){ //Children
        children.reserve(ncx);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
        
            auto child = std::make_unique<AMRGrid1D>(_nx, dx, ghosts);
            children.push_back(std::move(child));
        }
    }
}

AMRGrid2D::AMRGrid2D(int nx, int ny, double dx_, double dy_, int g):
    size_x(nx*dx_),size_y(ny*dy_), dx(dx_), dy(dy_), ghosts(validGhosts(g)),
    data(nx, ny, dx_, dy_, validGhosts(g)) {

    ncx = (nx / bin_size) + (nx % bin_size == 0 ? 0 : 1);
    ncy = (ny / bin_size) + (ny % bin_size == 0 ? 0 : 1);
    
    if (ncx > 1 || ncy > 1){ //Children
        children.reserve(ncx * ncy);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
            for(int j = 0; j<ncy; j++){
                int _ny = computeChildSize(ny, ncy, j);

                auto child = std::make_unique<AMRGrid2D>(_nx, _ny, dx, dy, ghosts);
                children.push_back(std::move(child));
            }
        }
    }
        
}

AMRGrid3D::AMRGrid3D(int nx, int ny, int nz, double dx_, double dy_, double dz_, int g):
    size_x(nx*dx_),size_y(ny*dy_), size_z(nz*dz_), dx(dx_), dy(dy_), dz(dz_), ghosts(validGhosts(g)),
    data(nx, ny, nz, dx_, dy_, dz_, validGhosts(g)) {

    ncx = (nx / bin_size) + (nx % bin_size == 0 ? 0 : 1);
    ncy = (ny / bin_size) + (ny % bin_size == 0 ? 0 : 1);
    ncz = (nz / bin_size) + (nz % bin_size == 0 ? 0 : 1);

    if (ncx > 1 || ncy > 1 || ncz > 1){ //Children
        children.reserve(ncx * ncy * ncz);
        for(int i = 0; i<ncx; i++){
            int _nx = computeChildSize(nx, ncx, i);
            for(int j = 0; j<ncy; j++){
                int _ny = computeChildSize(ny, ncy, j);
                for(int k = 0; k<ncz; k++){
                    int _nz = computeChildSize(nz, ncz, k);
                    
                    auto child = std::make_unique<AMRGrid3D>(_nx, _ny, _nz, dx, dy, dz, ghosts);
                    children.push_back(std::move(child));
                }
            }
        }
    }
        
}


//MARK: Grid Access
PrimitiveState& AMRGrid1D::operator[](int i){ return data[i];}
const PrimitiveState& AMRGrid1D::operator[](int i) const { return data[i]; }

PrimitiveState& AMRGrid2D::operator[](int i,int j){ return data[i,j];}
const PrimitiveState& AMRGrid2D::operator[](int i,int j) const { return data[i,j]; }

PrimitiveState& AMRGrid3D::operator[](int i,int j,int k){ return data[i,j,k];}
const PrimitiveState& AMRGrid3D::operator[](int i,int j,int k) const { return data[i,j,k]; }

int AMRGrid1D::getSize() const{ return data.getSize(); }
int AMRGrid1D::getGhosts() const{ return data.getGhosts(); }

int AMRGrid2D::getSizeX() const{ return data.getSizeX(); }
int AMRGrid2D::getSizeY() const{ return data.getSizeY(); }
int AMRGrid2D::getGhosts() const{ return data.getGhosts(); }

int AMRGrid3D::getSizeX() const{ return data.getSizeX(); }
int AMRGrid3D::getSizeY() const{ return data.getSizeY(); }
int AMRGrid3D::getSizeZ() const{ return data.getSizeZ(); }
int AMRGrid3D::getGhosts() const{ return data.getGhosts(); }


//MARK: Parent -> Child

void AMRGrid1D::pushToChildren(){
    if(ncx == 1) return; //This is the child
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        std::unique_ptr<AMRGrid1D>& child = children[zi];
        //Copy parent data to child
        for(int i = -ghosts; i < child->getSize() + ghosts; i++){
            (*child)[i] = data[i+x_offset];
        }
        child->boundary = Boundary::Ignore(); //Don't let the boundary API overwrite this
        child->pushToChildren(); //If children have children, make them sync too
           
        x_offset += child->getSize();
    }
}
void AMRGrid2D::pushToChildren(){
    if(ncx*ncy == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<AMRGrid2D>& child = children[ncy*zi + zj];
            int _ny = child->getSizeY();
            //Copy parent data to child
            for(int i = -ghosts; i < _nx + ghosts; i++){
                for(int j = -ghosts; j < _ny +ghosts; j++){
                    (*child)[i,j] = data[i+x_offset, j+y_offset];
                }
            }
            child->boundary = Boundary::Ignore(); //Don't let the boundary API overwrite this
            child->pushToChildren(); //If children have children, make them sync too
            
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}
void AMRGrid3D::pushToChildren(){
    if(ncx*ncy*ncz == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncz*ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            int _ny = children[ncz*(ncy*zi + zj)]->getSizeY();
            int z_offset = 0;
            for(int zk = 0; zk<ncz; zk++){
                std::unique_ptr<AMRGrid3D>& child = children[ncz*(ncy*zi + zj) + zk];
                int _nz = child->getSizeZ();
                //Copy parent data to child
                for(int i = -ghosts; i < _nx + ghosts; i++){
                    for(int j = -ghosts; j < _ny + ghosts; j++){
                        for(int k = -ghosts; k < _nz + ghosts; k++){
                            (*child)[i,j,k] = data[i+x_offset, j+y_offset, k+z_offset];
                        }
                    }
                }
                child->boundary = Boundary::Ignore(); //Don't let the boundary API overwrite this
                child->pushToChildren(); //If children have children, make them sync too
                
                z_offset += _nz;
            }
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}


//MARK: Child -> Parent
void AMRGrid1D::loadFromChildren(){
    if(ncx == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        std::unique_ptr<AMRGrid1D>& child = children[zi];
        child->loadFromChildren(); //If children have children, make them sync first
        //Copy child data to parent
        for(int i = 0; i < child->getSize(); i++) {
            data[i+x_offset] = (*child)[i];
        }
        x_offset += child->getSize();
    }
}

void AMRGrid2D::loadFromChildren(){
    if(ncx*ncy == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<AMRGrid2D>& child = children[ncy*zi + zj];
            child->loadFromChildren();//If children have children, make them sync first
            int _ny = child->getSizeY();
            //Copy child data to parent
            for(int i = 0; i < _nx; i++){
                for(int j=0; j < _ny; j++){
                    data[i+x_offset, j+y_offset] = (*child)[i,j];
                }
            }
            
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}

void AMRGrid3D::loadFromChildren(){
    if(ncx*ncy*ncz == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncz*ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            int _ny = children[ncz*(ncy*zi+zj)]->getSizeY();
            int z_offset = 0;
            for(int zi = 0; zi<ncx; zi++){
                std::unique_ptr<AMRGrid3D>& child = children[ncy*zi + zj];
                child->loadFromChildren();//If children have children, make them sync first
                int _nz = child->getSizeZ();
                //Copy child data to parent
                for(int i = 0; i < _nx; i++){
                    for(int j=0; j < _ny; j++){
                        for(int k=0; k < _nz; k++){
                            data[i+x_offset, j+y_offset, k+z_offset] = (*child)[i,j,k];
                        }
                    }
                }
                
                z_offset += _nz;
            }
            y_offset += _ny;
        }
        x_offset += _nx;
    }
}

//MARK: Advance
void AMRGrid1D::advance(double dt, bool check_cfl){
    if(ncx == 1 ){
        data.boundary = std::move(boundary);
        data.advance(dt, check_cfl);
        return;
    }
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(data);
        pushToChildren();
        
        //debugCheckChildGhosts();
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(data)) : dt;
        dt -= t1;
        //Execute the Advancement
        for(auto& child : children){
            DRARGONWING::launchParallel(child.get(), t1);
        }
        DRARGONWING::synchronize();
        
        //Copy Back
        loadFromChildren();
    }
}


void AMRGrid2D::advance(double dt, bool check_cfl){
    if(ncx*ncy == 1 ){
        data.boundary = std::move(boundary);
        data.advance(dt, check_cfl);
        return;
    }
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(data);
        pushToChildren();
        
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(data)) : dt;
        dt -= t1;
        //Execute the Advancement
        for(auto& child : children){
            DRARGONWING::launchParallel(child.get(), t1);
        }
        DRARGONWING::synchronize();//Wait for children to finish
        
        //Copy Back
        loadFromChildren();
    }
}


void AMRGrid3D::advance(double dt, bool check_cfl){
    if(ncx*ncy*ncz == 1 ){
        data.boundary = std::move(boundary);
        data.advance(dt, check_cfl);
        return;
    }
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(data);
        pushToChildren();
        
        //debugCheckChildGhosts();
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(data)) : dt;
        dt -= t1;
        //Execute the Advancement
        for(auto& child : children){
            DRARGONWING::launchParallel(child.get(), t1);
        }
        DRARGONWING::synchronize();
        
        //Copy Back
        loadFromChildren();
    }
}
