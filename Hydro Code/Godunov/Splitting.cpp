//
//  Untitled.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 16/05/2026.
//

#include "Godunov.hpp"
#include <cassert>


//MARK: 2D Array

Grid2D::Grid2D(int nx, int ny, int ghosts){
    this->nx = nx; this->ny = ny; this->ghosts = ghosts;
    w = new PrimitiveState[(nx+2*ghosts)*(ny+2*ghosts)];
}
PrimitiveState& Grid2D::operator[](int i, int j) {
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
    int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
    return w[m];
}
const PrimitiveState& Grid2D::operator[](int i, int j) const {
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
    int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
    return w[m];
}

Grid2D::~Grid2D(){ delete[] w; }



//MARK: 2D Sweep

void Grid2D::advance(double dt){
    advanceX(dt/2);
    advanceY(dt);
    advanceX(dt/2);
}

void Grid2D::advanceX(double dt){
    Grid1D W(nx,ghosts), _B1(nx,ghosts), _B2(nx,ghosts);
    for(int j=-ghosts; j<ny+ghosts; j++){
        for(int i=-ghosts; i<nx+ghosts; i++) W[i] = (*this)[i,j];

        W.advance(dt, dx, _B1, _B2);
        
        for(int i=-ghosts; i<nx+ghosts; i++) (*this)[i,j] = W[i];
    }
}
void Grid2D::advanceY(double dt){
    Grid1D W(ny,ghosts), _B1(ny,ghosts), _B2(ny,ghosts);
    for(int i=-ghosts; i<nx+ghosts; i++){
        for(int j=-ghosts; j<ny+ghosts; j++) {
            W[j] = (*this)[i,j];
            W[j].swapXY();
        }
        W.advance(dt, dy, _B1, _B2);
        for(int j=-ghosts; j<ny+ghosts; j++) {
            W[j].swapXY();
            (*this)[i,j] = W[j];
        }
    }
}




//MARK: 3D Array

Grid3D::Grid3D(int nx, int ny, int nz, int ghosts){
    this->nx = nx; this->ny = ny; this->nz = nz; this->ghosts = ghosts;
    w = new PrimitiveState[(nx+2*ghosts)*(ny+2*ghosts)*(nz+2*ghosts)];
}
PrimitiveState& Grid3D::operator[](int i, int j, int k) {
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
    assert(k + ghosts >= 0 && k < nz+ghosts);
    int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
    return w[m];
}
const PrimitiveState& Grid3D::operator[](int i, int j, int k) const {
    assert(i + ghosts >= 0 && i < nx+ghosts);
    assert(j + ghosts >= 0 && j < ny+ghosts);
    assert(k + ghosts >= 0 && k < nz+ghosts);
    int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
    return w[m];
}

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
    Grid1D W(nx,ghosts), _B1(nx,ghosts), _B2(nx,ghosts);
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int j=-ghosts; j<ny+ghosts; j++){
            for(int i=-ghosts; i<nx+ghosts; i++) W[i] = (*this)[i,j,k];
            
            W.advance(dt, dx, _B1, _B2);
            
            for(int i=-ghosts; i<nx+ghosts; i++) (*this)[i,j,k] = W[i];
        }
    }
}
void Grid3D::advanceY(double dt){
    Grid1D W(ny,ghosts), _B1(ny,ghosts), _B2(ny,ghosts);
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int i=-ghosts; i<nx+ghosts; i++) {
            for(int j=-ghosts; j<ny+ghosts; j++) {
                W[j] = (*this)[i,j,k];
                W[j].swapXY();
            }
            W.advance(dt, dy, _B1, _B2);
            for(int j=-ghosts; j<ny+ghosts; j++) {
                W[j].swapXY();
                (*this)[i,j,k] = W[j];
            }
       }
    }
}
void Grid3D::advanceZ(double dt){
    Grid1D W(nz,ghosts), _B1(nz,ghosts), _B2(nz,ghosts);
    for(int i=-ghosts; i<nx+ghosts; i++) {
        for(int j=-ghosts; j<ny+ghosts; j++) {
            for(int k=-ghosts; k<nz+ghosts; k++){
                W[k] = (*this)[i,j,k];
                W[k].swapXZ();
            }
            W.advance(dt, dz, _B1, _B2);
            for(int k=-ghosts; k<nz+ghosts; k++){
                W[k].swapXZ();
                (*this)[i,j,k] = W[k];
            }
       }
    }
}

