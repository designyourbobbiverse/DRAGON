//
//  Problem3D.cpp
//  DRAGON/Examples/SphericalBlast
//
//  Created by Bobbie Markwick on 8/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
using namespace DRAGON;

#include <cmath>        //For std::sqrt
#include "Constants.h"  //For gamma

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho = 1.0;
constexpr double p_amb = 1e-5;
constexpr double E_blast = 1.0;
constexpr int n = 256;
constexpr double r0 = 12.0/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Outflow();
    
    //Compute blast area
    int blast_cells = 0;
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            double x = (i + 0.5)/n  - 0.5;
            double y = (j + 0.5)/n - 0.5;
            double z = (k + 0.5)/n - 0.5;
            double r = std::sqrt(x*x + y*y + z*z);
            if(r < r0) blast_cells++;
        }
    }
    p_blast *= (n*n*n) / blast_cells;
    
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    PrimitiveState w;
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2,dz/2) corresponds to the [0,0,0] cell

    w.rho = rho;
    w.p = p_amb;
    
    double r = std::sqrt(x*x + y*y + z*z);
    if(r < r0) w.p += p_blast;
    
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    vec3 A;
    //This function is ignored in pure Hydro
    return A;
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
