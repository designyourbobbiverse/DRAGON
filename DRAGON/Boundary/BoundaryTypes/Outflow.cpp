//
//  BoundaryTypes/Outflow.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

//MARK: Constructors
Boundary::Outflow::Outflow(int faces_, bool corners, bool gated_):  GhostFill(faces_, corners), gated(gated_) {}
Boundary::Outflow::Outflow(std::string s, bool corners, bool gated): Outflow(face_mask(s),corners,gated) {}
Outflow Boundary::Outflow::Gated(int faces, bool corner_ghosts){ return Outflow(faces, corner_ghosts,true); }
Outflow Boundary::Outflow::Gated(std::string s, bool corners){ return Outflow(s,corners,true);}

//MARK: 1D
//Set each relevant ghost to the nearest physical cell, possibly subject to inflow gating
void Boundary::Outflow::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = grid[0];
            if(gated && grid[-g].v.x > 0) grid[-g].v.x = 0; //Kill inflows if applicable
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = grid[nx-1];
            if(gated && grid[nx-1+g].v.x < 0) grid[nx-1+g].v.x = 0; //Kill inflows if applicable
        }
    }
}

//MARK: 2D
//Set each relevant ghost to the nearest physical cell, possibly subject to inflow gating
void Boundary::Outflow::apply(Grid2D& grid) {
    //Calculate the bounds ahead of time
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);

    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[0,j];
                if(gated && grid[-g,j].v.x > 0) grid[-g,j].v.x = 0; //Kill inflows if applicable
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = grid[nx-1,j];
                if(gated && grid[nx-1+g,j].v.x < 0) grid[nx-1+g,j].v.x = 0; //Kill inflows if applicable
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,0];
                if(gated && grid[i,-g].v.y > 0) grid[i,-g].v.y = 0; //Kill inflows if applicable
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = grid[i,ny-1];
                if(gated && grid[i,ny-1+g].v.y < 0) grid[i,ny-1+g].v.y = 0; //Kill inflows if applicable
            }
        }
    }
//MARK: 2D MHD
    #ifdef MHD //Copy Transverse Fields
    auto& _A = grid._A();
    // A has one more physical point per dimension than w.
    if (faces & X_negative) {
        for (int j = j0; j <= jn; j++) {
            for (int g = 1; g <= ng; g++) {
                _A[-g,j] = 2*_A[1-g,j] - _A[2-g,j];
                _A[-g,j].x = _A[1-g,j].x;
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0; j <= jn; j++) {
            for (int g = 1; g <= ng; g++) {
                _A[nx+g,j] = 2*_A[nx+g-1,j] - _A[nx+g-2,j];
                _A[nx+g,j].x = _A[nx+g-1,j].x;
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0; i <= in; i++) {
            for (int g = 1; g <= ng; g++) {
                _A[i,-g] = 2*_A[i,1-g] - _A[i,2-g];
                _A[i,-g].y = _A[i,1-g].y;
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0; i <= in; i++) {
            for (int g = 1; g <= ng; g++) {
                _A[i,ny+g] = 2*_A[i,ny+g-1] - _A[i,ny+g-2];
                _A[i,ny+g].y = _A[i,ny+g-1].y;
            }
        }
    }
    #endif
}
//MARK: 3D
//Set each relevant ghost to the nearest physical cell, possibly subject to inflow gating
void Boundary::Outflow::apply(Grid3D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);
    
    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[0,j,k];
                    if(gated && grid[-g,j,k].v.x > 0) grid[-g,j,k].v.x = 0; //Kill inflows if applicable
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = grid[nx-1,j,k];
                    if(gated && grid[nx-1+g,j,k].v.x < 0) grid[nx-1+g,j,k].v.x = 0; //Kill inflows if applicable
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,0,k];
                    if(gated && grid[i,-g,k].v.y > 0) grid[i,-g,k].v.y = 0; //Kill inflows if applicable
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = grid[i,ny-1,k];
                    if(gated && grid[i,ny-1+g,k].v.y < 0) grid[i,ny-1+g,k].v.y = 0; //Kill inflows if applicable
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,0];
                    if(gated && grid[i,j,-g].v.z > 0) grid[i,j,-g].v.z = 0; //Kill inflows if applicable
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = grid[i,j,nz-1];
                    if(gated && grid[i,j,nz-1+g].v.z < 0) grid[i,j,nz-1+g].v.z = 0; //Kill inflows if applicable
                }
            }
        }
    }
    
//MARK: 3D MHD
#ifdef MHD //Copy Transverse Fields
    auto& _A = grid._A();
    // A is vertex-centred and has one more physical point per dimension than w.
    if (faces & X_negative) {
        for (int j = j0; j <= jn; j++) {
            for (int k = k0; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _A[-g,j,k] = 2*_A[1-g,j,k] - _A[2-g,j,k];
                    _A[-g,j,k].x = _A[1-g,j,k].x;
                }
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0; j <= jn; j++) {
            for (int k = k0; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _A[nx+g,j,k] = 2*_A[nx+g-1,j,k] - _A[nx+g-2,j,k];
                    _A[nx+g,j,k].x = _A[nx+g-1,j,k].x;
                }
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0; i <= in; i++) {
            for (int k = k0; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _A[i,-g,k] = 2*_A[i,1-g,k] - _A[i,2-g,k];
                    _A[i,-g,k].y = _A[i,1-g,k].y;
                }
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0; i <= in; i++) {
            for (int k = k0; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _A[i,ny+g,k] = 2*_A[i,ny+g-1,k] - _A[i,ny+g-2,k];
                    _A[i,ny+g,k].y = _A[i,ny+g-1,k].y;
                }
            }
        }
    }
    if (faces & Z_negative) {
        for (int i = i0; i <= in; i++) {
            for (int j = j0; j <= jn; j++) {
                for (int g = 1; g <= ng; g++) {
                    _A[i,j,-g] = 2*_A[i,j,1-g] - _A[i,j,2-g];
                    _A[i,j,-g].z = _A[i,j,1-g].z;
                }
            }
        }
    }
    if (faces & Z_positive) {
        for (int i = i0; i <= in; i++) {
            for (int j = j0; j <= jn; j++) {
                for (int g = 1; g <= ng; g++) {
                    _A[i,j,nz+g] = 2*_A[i,j,nz+g-1] - _A[i,j,nz+g-2];
                    _A[i,j,nz+g].z = _A[i,j,nz+g-1].z;
                }
            }
        }
    }
#endif
}
