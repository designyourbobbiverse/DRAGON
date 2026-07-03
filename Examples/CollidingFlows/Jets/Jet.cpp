//
//  BoundaryTypes/Jet.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/06/2026.
//

#include "Jets.hpp"
#include "Grid.hpp"
#include <cmath>

using namespace Boundary;

//MARK: Constructors
Jet::Jet(double rho, double v, double p, double r, int face): rho(rho), v(v), p(p), rj(r), jetface(face), GhostFill(0, true){ }
Jet::Jet(double rho, double v, double p, double r, std::string face): Jet(rho,v,p,r, face_mask(face)){ }

//MARK: 1D
//All applicable ghost cells are set to equal this->state
void Jet::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    
    PrimitiveState state;
    state.rho = rho;
    state.p = p;
    
    if (jetface & X_negative){
        state.v = {v,0,0};
        for(int g = 1; g <= ng; g++){
            grid[-g] = state;
        }
    }
    if (jetface & X_positive){
        state.v = {-v,0,0};
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = state;
        }
    }
}

//MARK: 2D
//All applicable ghost cells are set to equal this->state
void Jet::apply(Grid2D& grid) {
    //Calculate the bounds ahead of time
    double dx = grid.dx, dy = grid.dy;
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = fmax(0,floor(nx*0.5 - rj/dx)), in = fmin(nx,ceil(nx*0.5 + rj/dx));
    int j0 = fmax(0,floor(ny*0.5 - rj/dy)), jn = fmin(ny,ceil(ny*0.5 + rj/dy));

    PrimitiveState state;
    state.rho = rho;
    state.p = p;

    if (jetface & X_negative){
        state.v = {v,0,0};
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = state;
            }
        }
    }
    if (jetface & X_positive){
        state.v = {-v,0,0};
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = state;
            }
        }
    }
    if (jetface & Y_negative){
        state.v = {0,v,0};
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = state;
            }
        }
    }
    if (jetface & Y_positive){
        state.v = {0,-v,0};
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = state;
            }
        }
    }
}
//MARK: 3D
//All applicable ghost cells are set to equal this->state
void Jet::apply(Grid3D& grid) {
    //Calculate the bounds ahead of time
    double dx = grid.dx, dy = grid.dy, dz = grid.dz;
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = fmax(0,floor(nx*0.5 - rj/dx)), in = fmin(nx,ceil(nx*0.5 + rj/dx));
    int j0 = fmax(0,floor(ny*0.5 - rj/dy)), jn = fmin(ny,ceil(ny*0.5 + rj/dy));
    int k0 = fmax(0,floor(nz*0.5 - rj/dz)), kn = fmin(nz,ceil(nz*0.5 + rj/dz));

    PrimitiveState state;
    state.rho = rho;
    state.p = p;
    
    if (jetface & X_negative){
        state.v = {v,0,0};
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                double y = (j+0.5-ny*0.5)*dy, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(y*y+z*z)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = state;
                }
            }
        }
    }
    if (jetface & X_positive){
        state.v = {-v,0,0};
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                double y = (j+0.5-ny*0.5)*dy, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(y*y+z*z)>rj) continue;

                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = state;
                }
            }
        }
    }
    if (jetface & Y_negative){
        state.v = {0,v,0};
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                double x = (i+0.5-nx*0.5)*dx, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(x*x+z*z)>rj) continue;
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = state;
                }
            }
        }
    }
    if (jetface & Y_positive){
        state.v = {0,-v,0};
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                double x = (i+0.5-nx*0.5)*dx, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(x*x+z*z)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = state;
                }
            }
        }
    }
    if (jetface & Z_negative){
        state.v = {0,0,v};
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                double x = (i+0.5-nx*0.5)*dx, y = (j+0.5-ny*0.5)*dy;
                if(sqrt(x*x+y*y)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = state;
                }
            }
        }
    }
    if (jetface & Z_positive){
        state.v = {0,0,-v};
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                double x = (i+0.5-nx*0.5)*dx, y = (j+0.5-ny*0.5)*dy;
                if(sqrt(x*x+y*y)>rj) continue;

                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = state;
                }
            }
        }
    }
}
