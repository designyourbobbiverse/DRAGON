//
//  Problem.cpp
//  DRAGON/Examples/1D/Dai-Woodward
//
//  Created by Bobbie Markwick on 7/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
using namespace DRAGON;

typedef DistGrid1D MyGrid;

constexpr double rho_L = 1.08;
constexpr double rho_R = 1.0;
const double vx = 1.2;
const double vy = 0.01;
const double vz = 0.5;
const double Bx = 2.0;
const double ByL = 3.6;
const double ByR = 4.0;
const double Bz = 2.0;
constexpr double p_L = 0.95;
constexpr double p_R = 1.0;
const int nx = 1024*8;

Grid& Problem::makeProblem(){
    auto grid = new MyGrid(nx, 1.0/nx);
    grid->boundary = Boundary::Outflow();
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //dx/2 corresponds to the [0] cell. y and z will always be zero.
    PrimitiveState w;
    w.rho = x < 0.5 ? rho_L : rho_R;
    w.p = x < 0.5 ? p_L : p_R;
    w.v = x < 0.5 ? vec3{vx,vy,vz} : vec3{0,0,0};
    w.B = x < 0.5 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //This function is ignored in 1D MHD
    return {0,0,0};
}


void Problem::completeProblemInit(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    //Here you can do any initialization not covered by initialFluidState and initialMagneticPotential

}


void Problem::beforeCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the Godunov scheme. It isn't called for the initial frame
    
}


void Problem::afterCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called after the Godunov scheme but before the file write. It isn't called for the initial frame
    
}

void Problem::problemComplete(Grid& problem, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);

    //This is called only after the final time is reached.
}

