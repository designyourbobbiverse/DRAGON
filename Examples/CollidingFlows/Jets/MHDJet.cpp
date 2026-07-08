//
//  BoundaryTypes/Jet.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/06/2026.
//

#include "Jets.hpp"
#include "Grid.hpp"
#include <cmath>
#include "Constants.h"

using namespace Boundary;
#ifdef MHD
//MARK: Constructors
MHDJet::MHDJet(double rho, double v, double p, double beta, double rm, double rj, int face): Jet(rho,v,p,rj,face), beta(beta), rm(rm) { }
MHDJet::MHDJet(double rho, double v, double p, double beta, double rm, double rj, std::string face): Jet(rho,v,p,rj,face), beta(beta), rm(rm){ }

//MARK: 1D
void MHDJet::apply(Grid1D& grid) {
    //Construct a jet with reduced thermal pressure
    double pj = p * fabs(beta)/(1.0+ fabs(beta));
    std::swap(pj,p);
    Jet::apply(grid);
    std::swap(pj,p);
    
    int ng = grid.getGhosts();
    double B = sqrt(pj / fabs(beta)) * (beta/fabs(beta)) * sq8pi;

    if (jetface & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g].B.z = B;
        }
    }
    if (jetface & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g].B.z = B;
        }
    }
}

//MARK: 2D
//All applicable ghost cells are set to equal this->state
void MHDJet::apply(Grid2D& grid) {
    //Construct a jet with reduced thermal pressure
    double pj = p * fabs(beta)/(1.0+ fabs(beta));
    std::swap(pj,p);
    Jet::apply(grid);
    std::swap(pj,p);

    
    double dx = grid.dx, dy = grid.dy;
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = std::max(0.0,floor(nx*0.5 - rj/dx)), in = fmin(nx,ceil(nx*0.5 + rj/dx));
    int j0 = std::max(0.0,floor(ny*0.5 - rj/dy)), jn = fmin(ny,ceil(ny*0.5 + rj/dy));


    double B = sqrt(pj / fabs(beta)) * (beta/fabs(beta)) * sq8pi;

    if (jetface & X_negative){
        for(int j = j0; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j].B.z = B * (j*2 > ny ? 1 : -1);
            }
        }
    }
    if (jetface & X_positive){
        for(int j = j0; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j].B.z = B * (j*2 > ny ? 1 : -1);
            }
        }
    }
    if (jetface & Y_negative){
        for(int i = i0; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g].B.z = B * (i*2 > nx ? -1 : 1);
            }
        }
    }
    if (jetface & Y_positive){
        for(int i = i0; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g].B.z = B * (i*2 > nx ? -1 : 1);
            }
        }
    }
    
}


//MARK: 3D
//All applicable ghost cells are set to equal this->state
void MHDJet::apply(Grid3D& grid) {
    //Reuse existing logic to set rho, v, and p
    double pj = p / (1.0 + (rm*rm)/(fabs(beta)*rj*rj)); //Lower thermal pressure to account for B
    std::swap(pj,p);
    Jet::apply(grid);
    std::swap(pj,p); //Restore ambient pressure
    
    
    //Calculate the bounds ahead of time
    double dx = grid.dx, dy = grid.dy, dz = grid.dz;
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = -ng, in = nx + ng;
    int j0 = i0, jn = ny + ng;
    int k0 = i0, kn = nz + ng;
    
    auto& A = grid._A();
    
    double pmag = (pj/fabs(beta));
    double Bm = sqrt(pmag) * (beta < 0 ? -1 : 1) * sq8pi;
    double _rm = 1/rm;
    double _rm2 = _rm*_rm;
    double Am = Bm*rm*log(rj/rm);
    
    if (jetface & X_negative){
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                //Pressure on Body-centered cells
                double y = (j+0.5-ny*0.5)*dy, z = (k+0.5-nz*0.5)*dz;
                double r = sqrt(y*y+z*z);
                if (r < rm){
                    for(int g = 1; g <= ng; g++){
                        grid[-g,j,k].p += pmag * (1 - (r*r)*_rm2);
                    }
                }
                //A field is edge-centered
                y -= 0.5*dy; z -= 0.5*dz;
                r = sqrt(y*y+z*z);
                //Start by copying existing data
                for(int g = 1; g <= ng; g++) A[-g,j,k] = A[0,j,k];
                //Compute the A field
                double Ax;
                if(r <= rj) {
                    if(r < rm)  Ax = Am + (Bm*0.5)*(rm -r*r*_rm);
                    else Ax = Bm*rm*log(rj/r);
                } else Ax = 0;
                //Apply
                for(int g = 1; g <= ng; g++){ //A[0].x is interior, so exclude g=0
                    A[-g,j,k].x = Ax;
                }
            }
        }
    }
    if (jetface & X_positive){
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                //Pressure on Body-centered cells
                double y = (j+0.5-ny*0.5)*dy, z = (k+0.5-nz*0.5)*dz;
                double r = sqrt(y*y+z*z);
                if (r < rm){
                    for(int g = 1; g <= ng; g++){
                        grid[nx-1+g,j,k].p += pmag * (1 - (r*r)*_rm2);
                    }
                }
                //A field is edge-centered
                y -= 0.5*dy; z -= 0.5*dz;
                r = sqrt(y*y+z*z);
                //Start by copying existing data
                for(int g = 1; g <= ng; g++) A[nx+g,j,k] = A[nx,j,k];
                //Compute the A field
                double Ax;
                if(r <= rj) {
                    if(r < rm)  Ax = Am + (Bm*0.5)*(rm -r*r*_rm);
                    else Ax = Bm*rm*log(rj/r);
                } else Ax = 0;
                //Apply
                for(int g = 0; g <= ng; g++){ //A[nx].x is exterior, so include g=0
                    A[nx+g,j,k].x = Ax;
                }
            }
        }
    }
    if (jetface & Y_negative){
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                //Pressure on Body-centered cells
                double x = (i+0.5-nx*0.5)*dx, z = (k+0.5-nz*0.5)*dz;
                double r = sqrt(x*x+z*z);
                if (r < rm){
                    for(int g = 1; g <= ng; g++){
                        grid[i,-g,k].p += pmag * (1 - (r*r)*_rm2);
                    }
                }
                //A field is edge-centered
                x -= 0.5*dx; z -= 0.5*dz;
                r = sqrt(x*x+z*z);
                //Start by copying existing data
                for(int g = 1; g <= ng; g++) A[i,-g,k] = A[i,0,k];
                //Compute the A field
                double Ay;
                if(r <= rj) {
                    if(r < rm)  Ay = Am + (Bm*0.5)*(rm -r*r*_rm);
                    else Ay = Bm*rm*log(rj/r);
                }else Ay = 0;
                //Apply
                for(int g = 1; g <= ng; g++){ //A[0].y is interior, so exclude g=0
                    A[i,-g,k].y = Ay;
                }
            }
        }
    }
    if (jetface & Y_positive){
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                //Pressure on Body-centered cells
                double x = (i+0.5-nx*0.5)*dx, z = (k+0.5-nz*0.5)*dz;
                double r = sqrt(x*x+z*z);
                if (r < rm){
                    for(int g = 1; g <= ng; g++){
                        grid[i,ny-1+g,k].p += pmag * (1 - (r*r)*_rm2);
                    }
                }
                //A field is edge-centered
                x -= 0.5*dx; z -= 0.5*dz;
                r = sqrt(x*x+z*z);
                //Start by copying existing data
                for(int g = 1; g <= ng; g++) A[i,ny+g,k] = A[i,ny,k];
                //Compute the A field
                double Ay;
                if(r <= rj) {
                    if(r < rm)  Ay = Am + (Bm*0.5)*(rm -r*r*_rm);
                    else Ay = Bm*rm*log(rj/r);
                }else Ay = 0;
                //Apply
                for(int g = 0; g <= ng; g++){ //A[ny].y is exterior, so include g=0
                    A[i,ny+g,k].y = Ay;
                }
            }
        }
    }
    if (jetface & Z_negative){
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                //Pressure on Body-centered cells
                double x = (i+0.5-nx*0.5)*dx, y = (j+0.5-ny*0.5)*dy;
                double r = sqrt(x*x+y*y);
                if (r < rm){
                    for(int g = 1; g <= ng; g++){
                        grid[i,j,-g].p += pmag * (1 - (r*r)*_rm2);
                    }
                }
                //A field is edge-centered
                x -= 0.5*dx; y -= 0.5*dy;
                r = sqrt(x*x+y*y);
                //Start by copying existing data
                for(int g = 1; g <= ng; g++) A[i,j,-g] = A[i,j,0];
                //Compute the A field
                double Az;
                if(r <= rj) {
                    if(r < rm)  Az = Am + (Bm*0.5)*(rm -r*r*_rm);
                    else Az = Bm*rm*log(rj/r);
                }  else Az = 0;
                //Apply
                for(int g = 1; g <= ng; g++){ //A[0].z is interior, so exclude g=0
                    A[i,j,-g].z = Az;
                }

            }
        }
    }
    if (jetface & Z_positive){
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                //Pressure on Body-centered cells
                double x = (i+0.5-nx*0.5)*dx, y = (j+0.5-ny*0.5)*dy;
                double r = sqrt(x*x+y*y);
                if (r < rm){
                    for(int g = 1; g <= ng; g++){
                        grid[i,j,nz-1+g].p += pmag * (1 - (r*r)*_rm2);
                    }
                }
                //A field is edge-centered
                x -= 0.5*dx; y -= 0.5*dy;
                r = sqrt(x*x+y*y);
                //Start by copying existing data
                for(int g = 1; g <= ng; g++) A[i,j,nz+g] = A[i,j,nz];
                //Compute the A field
                double Az;
                if(r <= rj) {
                    if(r < rm)  Az = Am + (Bm*0.5)*(rm -r*r*_rm);
                    else Az = Bm*rm*log(rj/r);
                } else Az = 0;
                //Apply
                for(int g = 0; g <= ng; g++){ //A[nz].z is exterior, so include g=0
                    A[i,j,nz+g].z = Az;
                }
            }
        }
    }
    
}
#endif
