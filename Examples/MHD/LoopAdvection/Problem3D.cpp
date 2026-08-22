//
//  Problem3D.cpp
//  DRAGON/Examples/MHD/LoopAdvection
//
//  Created by Bobbie Markwick on 20/07/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

#include <cmath>           //For std::sqrt
#include <algorithm>       //For std::max
#include "DragonHoard.hpp" //For loading initial data to compare to final
#include <iostream>        //For error output

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1.0;
constexpr double p_amb = 5.0;
constexpr double vx = 2.0;
constexpr double vy = 1.0;
constexpr double vz = 1.0;

constexpr bool gaussian = false;
constexpr double r0 = gaussian ? 0.1 : 0.2;
constexpr double B0 = 1e-6;

constexpr int n = 128;
constexpr double dx = 2.0/n;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n/2,n/2, dx, dx, dx);
    grid->boundary = Boundary::Periodic();
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2, dy/2, dz/2) corresponds to the [0,0,0] cell
    PrimitiveState w;
    w.rho =  rho0;
    w.p = p_amb;
    w.v = {vx,vy,vz};
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //Initialize the vector potential at point (x,y,z)
        // (0,0,0) corresponds to the [0,0,0] cell.
    //As such, we need to adjust them to be center-relative
    x -= 1.0; y -= 0.5;
    //Magnetic Fields will be initialized from this potential to ensure div B = 0
    double r = std::sqrt(x*x + y*y);
    double Az = gaussian ? B0*exp(-0.5*std::pow(r/r0,2)) : B0*std::max(r0-r,0.0);
    return {0,0,Az};
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
    
    MyGrid& expected =  *dynamic_cast<MyGrid*>(&makeProblem());
    double ig; int nore;
    DRAGONHOARD::loadFromFile(expected, ig,nore, DRAGONHOARD::output_base_name + "_" + DRAGONHOARD::cycle_string(0));
    
    vec3 L1;
    for(int i=0; i<n;i++){
        for(int j=0; j<n/2; j++){
            for(int k=0; k<n/2; k++){
                vec3 err = grid[i,j,k].B - expected[i,j,k].B;
                L1 += vec3{std::abs(err.x), std::abs(err.y), std::abs(err.z)} ;
            }
        }
    }
    
    L1 /= (n*n*n/4);
    std::cout<<"L1 error: "<< std::sqrt(L1*L1) <<"\n";
}
