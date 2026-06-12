//
//  GridGroup.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#include "GridGroup.hpp"
#include "CFL.hpp"

using namespace Boundary;

//MARK: Setup
GridGroup2D::GridGroup2D(int nx, int ny, double dx, double dy, int ghosts) {
#ifdef MUSCL_Hancock
    if(ghosts < 2) ghosts = 2;
#else
    if(ghosts < 1) ghosts = 1;
#endif
    for(int n = 0; n<16; n++){
        zones[n] = std::make_unique<Grid2D>(nx/4, ny/4, dx, dy, ghosts);
    }
    reloadLinks();
}
void GridGroup2D::reloadLinks(){
    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            (*zones[4*i+j]).boundary = Boundary::Reflective();
        }
    }
    for(int i = 0; i<4; i++){
        for(int j = 0; j<3; j++){
            (*zones[4*i+j]).boundary += GridLink(zones[4*i+j+1].get(), Y_positive);
            (*zones[4*i+j+1]).boundary += GridLink(zones[4*i+j].get(), Y_negative);
        }
    }
    for(int i = 0; i<3; i++){
        for(int j = 0; j<4; j++){
            (*zones[4*i+j]).boundary += GridLink(zones[4*i+j+4].get(), X_positive);
            (*zones[4*i+j+4]).boundary += GridLink(zones[4*i+j].get(), X_negative);
        }
    }
}


//MARK: Grid Access
PrimitiveState& GridGroup2D::operator[](int i, int j) {
    int znx =  (*zones[0]).getSizeX(), zny =  (*zones[0]).getSizeX();
    int _i = i % znx, zi = i/znx;
    int _j = j % zny, zj = j/zny;
    return (*zones[4*zi + zj])[_i,_j];
}
const PrimitiveState& GridGroup2D::operator[](int i, int j) const {
    int znx =  (*zones[0]).getSizeX(), zny =  (*zones[0]).getSizeX();
    int _i = i % znx, zi = i/znx;
    int _j = j % zny, zj = j/zny;
    return (*zones[4*zi + zj])[_i,_j];
}
int GridGroup2D::getSizeX() const { return (*zones[0]).getSizeX()*4; }
int GridGroup2D::getSizeY() const { return  (*zones[0]).getSizeY()*4; }
int GridGroup2D::getGhosts() const { return (*zones[0]).getGhosts(); }

//MARK: Advance
void GridGroup2D::advance(double dt){
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        reloadLinks();
        //CFL Time Constraint
        double t1 = dt;
        for(int n=0; n<16; n++){
            t1 = std::min(t1,CFL::cfl_time(*zones[n]));
        }
        dt -= t1;
        //Execute the Advancement
        for(int n=0; n<16; n++){
            zones[n]->advance(t1,false);
        }
    }    
}
