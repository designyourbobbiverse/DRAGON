//
//  AMRGrid.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#include "AMRGrid.hpp"
#include "ParentGridFill.hpp"
#include "DragonWing.hpp"
#include "CFL.hpp"
#include <math.h>
#include <cassert>


int computeChildSize(int nx, int ncx, int i){
    if(nx <= bin_size) return nx;
    int _nx = bin_size;
    if (i == 0 || i+1 == ncx) {
        _nx -=  (ncx*bin_size - nx)/2;
        if (i==0 && nx%2==1) _nx -= 1;
    }
    return _nx;
}

static int validGhosts(int g){
#ifdef MUSCL_Hancock
    return std::max(g, 3);
#else
    return std::max(g, 2);
#endif
}


AMRGrid2D::AMRGrid2D(int nx, int ny, double dx_, double dy_, int g):
    size_x(nx*dx_),size_y(ny*dy_), dx(dx_), dy(dy_), ghosts(validGhosts(g)),
    data(nx, ny, dx_, dy_, validGhosts(g)) {

    ncx = (nx / bin_size) + (nx % bin_size == 0 ? 0 : 1);
    ncy = (ny / bin_size) + (ny % bin_size == 0 ? 0 : 1);
    
        
    if (ncx > 1 || ncy > 1){ //Children
        children.reserve(ncx * ncy);
        int x_offset = 0;
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
AMRGrid2D::~AMRGrid2D(){  }
PrimitiveState& AMRGrid2D::operator[](int i,int j){ return data[i,j];}
const PrimitiveState& AMRGrid2D::operator[](int i,int j) const { return data[i,j]; }
int AMRGrid2D::getSizeX() const{ return data.getSizeX(); }
int AMRGrid2D::getSizeY() const{ return data.getSizeY(); }
int AMRGrid2D::getGhosts() const{ return data.getGhosts(); }


void AMRGrid2D::loadFromChildren(){
    if(ncx*ncy == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<AMRGrid2D>& child = children[ncy*zi + zj];
            int _ny = child->getSizeY();
            
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
void AMRGrid2D::pushToChildren(){
    if(ncx*ncy == 1) return;
    int x_offset = 0;
    for(int zi = 0; zi<ncx; zi++){
        int _nx = children[ncy*zi]->getSizeX();
        int y_offset = 0;
        for(int zj = 0; zj<ncy; zj++){
            std::unique_ptr<AMRGrid2D>& child = children[ncy*zi + zj];
            int _ny = child->getSizeY();
            
            for(int i = -ghosts; i < _nx + ghosts; i++){
                for(int j = -ghosts; j < _ny +ghosts; j++){
                    (*child)[i,j] = data[i+x_offset, j+y_offset];
                }
            }
            child->boundary = Boundary::Ignore();
            child->pushToChildren();
            
            y_offset += _ny;
        }
        x_offset += _nx;
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


