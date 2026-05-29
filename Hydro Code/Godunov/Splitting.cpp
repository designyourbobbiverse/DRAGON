//
//  Untitled.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 16/05/2026.
//

#include "Godunov.hpp"
#include "Boundary.hpp"
#include <cassert>


//MARK: 2D Array

Grid2D::Grid2D(int nx_, int ny_, double dx_, double dy_, int g_):
    nx(nx_), ny(ny_), dx(dx_), dy(dy_), ghosts(g_) {
    w = new PrimitiveState[(nx+2*ghosts)*(ny+2*ghosts)];
}
PrimitiveState& Grid2D::operator[](int i, int j) {
#ifdef TESTMODE
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
#endif
    int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
    return w[m];
}
const PrimitiveState& Grid2D::operator[](int i, int j) const {
#ifdef TESTMODE
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
#endif
    int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
    return w[m];
}
int Grid2D::getSizeX(){ return nx; }
int Grid2D::getSizeY(){ return ny; }
int Grid2D::getGhosts(){ return ghosts; }
Grid2D::~Grid2D(){ delete[] w; }



//MARK: 2D Sweep

void Grid2D::advance(double dt){
    advanceX(dt/2);
    advanceY(dt);
    advanceX(dt/2);
}

void Grid2D::advanceX(double dt){
    if(boundary) boundary->apply(*this);
    
    Grid1D W(nx, dx, ghosts), _B1(nx, dx, ghosts), _B2(nx, dx, ghosts);
    for(int j=-ghosts; j<ny+ghosts; j++){
        for(int i=-ghosts; i<nx+ghosts; i++) W[i] = (*this)[i,j];

        W.god_sweep(dt, _B1, _B2);
        
        for(int i=-ghosts; i<nx+ghosts; i++) (*this)[i,j] = W[i];
    }
}
void Grid2D::advanceY(double dt){
    if(boundary) boundary->apply(*this);

    Grid1D W(ny, dy, ghosts), _B1(ny, dy, ghosts), _B2(ny,dy, ghosts);
    for(int i=-ghosts; i<nx+ghosts; i++){
        for(int j=-ghosts; j<ny+ghosts; j++) {
            W[j] = (*this)[i,j];
            W[j].swapXY();
        }
        W.god_sweep(dt, _B1, _B2);
        for(int j=-ghosts; j<ny+ghosts; j++) {
            W[j].swapXY();
            (*this)[i,j] = W[j];
        }
    }
}




//MARK: 3D Array

Grid3D::Grid3D(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_, int g_):
    nx(nx_), ny(ny_), nz(nz_), dx(dx_), dy(dy_), dz(dz_), ghosts(g_) {
    w = new PrimitiveState[(nx+2*ghosts)*(ny+2*ghosts)*(nz+2*ghosts)];
}
PrimitiveState& Grid3D::operator[](int i, int j, int k) {
#ifdef TESTMODE
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
    assert(k + ghosts >= 0 && k < nz+ghosts);
#endif
    int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
    return w[m];
}
const PrimitiveState& Grid3D::operator[](int i, int j, int k) const {
#ifdef TESTMODE
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
    assert(k + ghosts >= 0 && k < nz+ghosts);
#endif
    int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
    return w[m];
}
int Grid3D::getSizeX(){ return nx; }
int Grid3D::getSizeY(){ return ny; }
int Grid3D::getSizeZ(){ return nz; }
int Grid3D::getGhosts(){ return ghosts; }
Grid3D::~Grid3D(){ delete[] w; }


//MARK: 3D Sweep

void Grid3D::advance(double dt){
    advanceX(dt/2);
    advanceY(dt/2);
    advanceZ(dt);
    advanceY(dt/2);
    advanceX(dt/2);
}

void Grid3D::advanceX(double dt){
    if(boundary) boundary->apply(*this);

    Grid1D W(nx,dx,ghosts), _B1(nx,dx,ghosts), _B2(nx,dx,ghosts);
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int j=-ghosts; j<ny+ghosts; j++){
            for(int i=-ghosts; i<nx+ghosts; i++) W[i] = (*this)[i,j,k];
            
            W.god_sweep(dt, _B1, _B2);
            
            for(int i=-ghosts; i<nx+ghosts; i++) (*this)[i,j,k] = W[i];
        }
    }
}
void Grid3D::advanceY(double dt){
    if(boundary) boundary->apply(*this);

    Grid1D W(ny,dy,ghosts), _B1(ny,dy,ghosts), _B2(ny,dy,ghosts);
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int i=-ghosts; i<nx+ghosts; i++) {
            for(int j=-ghosts; j<ny+ghosts; j++) {
                W[j] = (*this)[i,j,k];
                W[j].swapXY();
            }
            W.god_sweep(dt, _B1, _B2);
            for(int j=-ghosts; j<ny+ghosts; j++) {
                W[j].swapXY();
                (*this)[i,j,k] = W[j];
            }
       }
    }
}
void Grid3D::advanceZ(double dt){
    if(boundary) boundary->apply(*this);

    Grid1D W(nz,dz,ghosts), _B1(nz,dz,ghosts), _B2(nz,dz,ghosts);
    for(int i=-ghosts; i<nx+ghosts; i++) {
        for(int j=-ghosts; j<ny+ghosts; j++) {
            for(int k=-ghosts; k<nz+ghosts; k++){
                W[k] = (*this)[i,j,k];
                W[k].swapXZ();
            }
            W.god_sweep(dt, _B1, _B2);
            for(int k=-ghosts; k<nz+ghosts; k++){
                W[k].swapXZ();
                (*this)[i,j,k] = W[k];
            }
       }
    }
}

