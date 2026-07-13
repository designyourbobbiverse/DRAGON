//
//  BoundaryTypes/Reflective.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

#ifdef MHD
Boundary::Reflective::Reflective(int faces_, bool corners, bool conductive):  GhostFill(faces_, corners), conductive(conductive) {}
Boundary::Reflective::Reflective(std::string s, bool corners, bool conductive) : Reflective(face_mask(s),corners,conductive) {}
#else
Boundary::Reflective::Reflective(int faces_, bool corners):  GhostFill(faces_, corners) {}
Boundary::Reflective::Reflective(std::string s, bool corners) : Reflective(face_mask(s),corners) {}
#endif

//MARK: 1D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = grid[g-1];
            grid[-g].v.x *= -1;//Mirror Normal Velocity
            #ifdef MHD
            if(conductive) { //Mirror Transverse Magnetic Fields
                grid[-g].B *= -1;
                grid[-g].B.x *= -1;
            } else {
                grid[-g].B = grid[0].B;
            }
            #endif
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = grid[nx-g];
            grid[nx-1+g].v.x *= -1;//Mirror Normal Velocity
            #ifdef MHD
            if(conductive) { //Mirror Transverse Magnetic Fields
                grid[nx-1+g].B *= -1;
                grid[nx-1+g].B.x *= -1;
            } else {
                grid[nx-1+g].B = grid[nx-1].B;
            }
            #endif
        }
    }
}

//MARK: 2D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid2D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);


    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[g-1,j];
                grid[-g,j].v.x *= -1;//Mirror Normal Velocity
                #ifdef MHD
                if(conductive) {
                    grid[-g,j].B *= -1;
                    grid[-g,j].B.x *= -1;
                } else {
                    grid[-g,j].B =  grid[0,j].B;
                }
                #endif
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = grid[nx-g,j];
                grid[nx-1+g,j].v.x *= -1;//Mirror Normal Velocity
                #ifdef MHD
                if(conductive) {
                    grid[nx-1+g,j].B *= -1;
                    grid[nx-1+g,j].B.x *= -1;
                } else {
                    grid[nx-1+g,j].B =  grid[0,j].B;
                }
                #endif
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,g-1];
                grid[i,-g].v.y *= -1;//Mirror Normal Velocity
                #ifdef MHD
                if(conductive) {
                    grid[i,-g].B *= -1;
                    grid[i,-g].B.y *= -1;
                } else {
                    grid[i,-g].B =  grid[i,0].B;
                }
                #endif
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = grid[i,ny-g];
                grid[i,ny-1+g].v.y *= -1;//Mirror Normal Velocity
                #ifdef MHD //2D is weird and computes Bz on the body, but Bx/By via the edge potential
                if(conductive) {
                    grid[i,ny-1+g].B *= -1;
                    grid[i,ny-1+g].B.y *= -1;
                } else {
                    grid[i,ny-1+g].B =  grid[i,0].B;
                }
                #endif
            }
        }
    }
//MARK: 2D MHD
    #ifdef MHD
    auto& _A = grid._A();
    // A has one more physical point per dimension than w.
    if (faces & X_negative){
        for(int j = j0 ; j <= jn; j++){
            for(int g = 1; g <= ng; g++){
                if (conductive) {//Mirror Transverse Magnetic Fields
                    _A[-g,j] = _A[g,j]; //A reflects over zero, w reflects over -1/2
                } else {
                    _A[-g,j].z = 2*_A[1-g,j].z - _A[2-g,j].z;
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j <= jn; j++){
            for(int g = 1; g <= ng; g++){
                if (conductive) {//Mirror Transverse Magnetic Fields
                    _A[nx+g,j] = _A[nx-g,j];//A reflects over zero, w reflects over -1/2
                } else {
                    _A[nx+g,j].z = 2*_A[nx+g-1,j].z - _A[nx+g-2,j].z;
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i <= in; i++){
            for(int g = 1; g <= ng; g++){
                if (conductive) {
                    _A[i,-g] = _A[i,g];//A reflects over zero, w reflects over -1/2
                } else {
                    _A[i,-g].z = 2*_A[i,1-g].z - _A[i,2-g].z;
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i <= in; i++){
            for(int g = 1; g <= ng; g++){
                if (conductive) {
                    _A[i,ny+g] = _A[i,ny-g];//A reflects over zero, w reflects over -1/2
                } else {
                    _A[i,ny+g].z = 2*_A[i,ny+g-1].z - _A[i,ny+g-2].z;
                }
            }
        }
    }
    #endif
    
}
//MARK: 3D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid3D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);


    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[g-1,j,k];
                    grid[-g,j,k].v.x *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    if(conductive) {
                        grid[-g,j,k].B *= -1;
                        grid[-g,j,k].B.x *= -1;
                    } else {
                        grid[-g,j,k].B =  grid[0,j,k].B;
                    }
                    #endif
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k].v.x *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    if(conductive) {
                        grid[nx-1+g,j,k].B *= -1;
                        grid[nx-1+g,j,k].B.x *= -1;
                    } else {
                        grid[nx-1+g,j,k].B =  grid[nx-1,j,k].B;
                    }
                    #endif
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,g-1,k];
                    grid[i,-g,k].v.y *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    if(conductive) {
                        grid[i,-g,k].B *= -1;
                        grid[i,-g,k].B.y *= -1;
                    } else {
                        grid[i,-g,k].B =  grid[i,0,k].B;
                    }
                    #endif
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k].v.y *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    if(conductive) {
                        grid[i,ny-1+g,k].B *= -1;
                        grid[i,ny-1+g,k].B.y *= -1;
                    } else {
                        grid[i,ny-1+g,k].B =  grid[i,ny-1,k].B;
                    }
                    #endif
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,g-1];
                    grid[i,j,-g].v.z *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    if(conductive) {
                        grid[i,j,-g].B *= -1;
                        grid[i,j,-g].B.z *= -1;
                    } else {
                        grid[i,j,-g].B =  grid[i,j,0].B;
                    }
                    #endif
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g].v.z *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    if(conductive) {
                        grid[i,j,nz-1+g].B *= -1;
                        grid[i,j,nz-1+g].B.z *= -1;
                    } else {
                        grid[i,j,nz-1+g].B =  grid[i,j,nz-1].B;
                    }
                    #endif
                }
            }
        }
    }
//MARK: 3D MHD
    #ifdef MHD
    auto& _A = grid._A();
    // A has one more physical point per dimension than w.

    if (faces & X_negative){
        for(int j = j0 ; j <= jn; j++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        _A[-g,j,k] = _A[g,j,k]; //A reflects over zero, w reflects over -1/2
                        _A[-g,j,k].x = -_A[g-1,j,k].x; //Except normal A does refelct over  -1/2
                    } else {
                        _A[-g,j,k] = 2*_A[1-g,j,k] - _A[2-g,j,k];
                        _A[-g,j,k].x = _A[1-g,j,k].x;
                    }
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j <= jn; j++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        _A[nx+g,j,k] = _A[nx-g,j,k]; //A reflects over zero, w reflects over -1/2
                        _A[nx-1+g,j,k].x = -_A[nx-g,j,k].x; //Except normal A does reflect over -1/2
                    } else {
                        _A[nx+g,j,k] = 2*_A[nx+g-1,j,k] - _A[nx+g-2,j,k];
                        _A[nx+g,j,k].x = _A[nx+g-1,j,k].x;
                    }
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i <= in; i++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        _A[i,-g,k] = _A[i,g,k]; //A reflects over zero, w reflects over -1/2
                        _A[i,-g,k].y = -_A[i,g-1,k].y; //Except normal A does reflect over -1/2
                    } else {
                        _A[i,-g,k] = 2*_A[i,1-g,k] - _A[i,2-g,k];
                        _A[i,-g,k].y = _A[i,1-g,k].y;
                    }
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i <= in; i++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        _A[i,ny+g,k] = _A[i,ny-g,k]; //A reflects over zero, w reflects over -1/2
                        _A[i,ny-1+g,k].y = -_A[i,ny-g,k].y; //Except normal A does reflect over -1/2
                    } else {
                        _A[i,ny+g,k] = 2*_A[i,ny+g-1,k] - _A[i,ny+g-2,k];
                        _A[i,ny+g,k].y = _A[i,ny+g-1,k].y;
                    }
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i <= in; i++){
            for(int j = j0 ; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        _A[i,j,-g] = _A[i,j,g]; //A reflects over zero, w reflects over -1/2
                        _A[i,j,-g].z = -_A[i,j,g-1].z; //Except normal A does reflect over -1/2
                    } else {
                        _A[i,j,-g] = 2*_A[i,j,1-g] - _A[i,j,2-g];
                        _A[i,j,-g].z = _A[i,j,1-g].z;
                    }
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i <= in; i++){
            for(int j = j0 ; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        _A[i,j,nz+g] = _A[i,j,nz-g]; //A reflects over zero, w reflects over -1/2
                        _A[i,j,nz-1+g].z = -_A[i,j,nz-g].z; //Except normal A does reflect over -1/2
                    } else {
                        _A[i,j,nz+g] = 2*_A[i,j,nz+g-1] - _A[i,j,nz+g-2];
                        _A[i,j,nz+g].z = _A[i,j,nz+g-1].z;
                    }
                }
            }
        }
    }
    #endif
}
