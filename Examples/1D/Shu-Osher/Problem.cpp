//
//  Problem_Shu_Osher.cpp
//  DRAGON/Examples/1D/Shu-Osher
//
//  Created by Bobbie Markwick on 7/07/2026.
//  Shu and Osher (1989). https://doi.org/10.1016/0021-9991(89)90222-2
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

typedef DistGrid1D MyGrid;//Choose the dimension of your grid here

constexpr double rho_L = 3.857143;
constexpr double rho_R = 1.0;
constexpr double p_L = 10.333333;
constexpr double p_R = 1.0;
constexpr double vxL = 2.629369;
constexpr int nx = 1280*8;
constexpr double dx = 10.0/nx;

Grid& Problem::makeProblem(){
    auto grid = new MyGrid(nx, dx);
    grid->boundary = Boundary::Outflow();
    return *grid;
}

PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //dx/2 corresponds to the [0] cell. y and z will always be zero.
    PrimitiveState w;
    w.rho = x < 1.0 ? rho_L : rho_R;
    w.p = x < 1.0 ? p_L : p_R;
    w.v = x < 1.0 ? vec3{vxL,0,0} : vec3{0,0,0};
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //This function is ignored in pure hydro
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
