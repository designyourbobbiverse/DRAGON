//
//  Boundary/Reflective.cpp
//  DRAGON
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
                if (conductive) {//Mirror Transverse Magnetic Fields
                    grid.getA()[-g,j] = grid.getA()[g,j]; //A reflects over zero, w reflects over -1/2
                    grid[-g,j].B.z *= -1; //2D is weird and computes Bz on the body, but Bx/By via the edge potential
                } else {
                    grid.getA()[-g,j].z = 2*grid.getA()[1-g,j].z - grid.getA()[2-g,j].z;
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
                if (conductive) {//Mirror Transverse Magnetic Fields
                    grid.getA()[nx+g,j] = grid.getA()[nx-g,j];//A reflects over zero, w reflects over -1/2
                    grid[nx-1+g,j].B.z *= -1; //2D is weird and computes Bz on the body, but Bx/By via the edge potential
                } else {
                    grid.getA()[nx+g,j].z = 2*grid.getA()[nx+g-1,j].z - grid.getA()[nx+g-2,j].z;
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
                if (conductive) {
                    grid.getA()[i,-g] = grid.getA()[i,g];//A reflects over zero, w reflects over -1/2
                    grid[i,-g].B.z *= -1; //2D is weird and computes Bz on the body, but Bx/By via the edge potential
                } else {
                    grid.getA()[i,-g].z = 2*grid.getA()[i,1-g].z - grid.getA()[i,2-g].z;
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
                #ifdef MHD//Mirror Transverse Magnetic Fields
                if (conductive) {
                    grid.getA()[i,ny+g] = grid.getA()[i,ny-g];//A reflects over zero, w reflects over -1/2
                    grid[i,ny-1+g].B.z *= -1; //2D is weird and computes Bz on the body, but Bx/By via the edge potential
                } else {
                    grid.getA()[i,ny+g].z = 2*grid.getA()[i,ny+g-1].z - grid.getA()[i,ny+g-2].z;
                }
                #endif
            }
        }
    }
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
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        grid.getA()[-g,j,k] = grid.getA()[g,j,k]; //A reflects over zero, w reflects over -1/2
                        grid.getA()[-g,j,k].x = grid.getA()[g-1,j,k].x; //Except normal A does refelct over  -1/2
                    } else {
                        grid.getA()[-g,j,k] = 2*grid.getA()[1-g,j,k] - grid.getA()[2-g,j,k];
                        grid.getA()[-g,j,k].x = grid.getA()[1-g,j,k].x;
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
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        grid.getA()[nx+g,j,k] = grid.getA()[nx-g,j,k]; //A reflects over zero, w reflects over -1/2
                        grid.getA()[nx-1+g,j,k].x = -grid.getA()[nx-g,j,k].x; //Except normal A does reflect over -1/2
                    } else {
                        grid.getA()[nx+g,j,k] = 2*grid.getA()[nx+g-1,j,k] - grid.getA()[nx+g-2,j,k];
                        grid.getA()[nx+g,j,k].x = grid.getA()[nx+g-1,j,k].x;
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
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        grid.getA()[i,-g,k] = grid.getA()[i,g,k]; //A reflects over zero, w reflects over -1/2
                        grid.getA()[i,-g,k].y = -grid.getA()[i,g-1,k].y; //Except normal A does reflect over -1/2
                    } else {
                        grid.getA()[i,-g,k] = 2*grid.getA()[i,1-g,k] - grid.getA()[i,2-g,k];
                        grid.getA()[i,-g,k].y = grid.getA()[i,1-g,k].y;
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
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        grid.getA()[i,ny+g,k] = grid.getA()[i,ny-g,k]; //A reflects over zero, w reflects over -1/2
                        grid.getA()[i,ny-1+g,k].x = -grid.getA()[i,ny-g,k].y; //Except normal A does reflect over -1/2
                    } else {
                        grid.getA()[i,ny+g,k] = 2*grid.getA()[i,ny+g-1,k] - grid.getA()[i,ny+g-2,k];
                        grid.getA()[i,ny+g,k].y = grid.getA()[i,ny+g-1,k].y;
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
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        grid.getA()[i,j,-g] = grid.getA()[i,j,g]; //A reflects over zero, w reflects over -1/2
                        grid.getA()[i,j,-g].z = -grid.getA()[i,j,g-1].z; //Except normal A does reflect over -1/2
                    } else {
                        grid.getA()[i,j,-g] = 2*grid.getA()[i,j,1-g] - grid.getA()[i,j,2-g];
                        grid.getA()[i,j,-g].z = grid.getA()[i,j,1-g].z;
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
                    if (conductive) { //Mirror Transverse Magnetic Fields
                        grid.getA()[i,j,nz+g] = grid.getA()[i,j,nz-g]; //A reflects over zero, w reflects over -1/2
                        grid.getA()[i,j,nz-1+g].z = -grid.getA()[i,j,nz-g].z; //Except normal A does reflect over -1/2
                    } else {
                        grid.getA()[i,j,nz+g] = 2*grid.getA()[i,j,nz+g-1] - grid.getA()[i,j,nz+g-2];
                        grid.getA()[i,j,nz+g].z = grid.getA()[i,j,nz+g-1].z;
                    }
                    #endif
                }
            }
        }
    }
}
