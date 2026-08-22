//
//  Problem2D.cpp
//  DRAGON/Examples/Hydro/Diagonal
//
//  Created by Bobbie Markwick on 09/07/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

#include <cmath> //For std::sin
#include <iostream> //For std::cout at completion

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1.0;
constexpr double rho1 = 0.1;
constexpr double p_amb = 1.0;
constexpr int n = 512;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //dx/2, dy/2 corresponds to the [0,0] cell. z will always be zero.
    PrimitiveState w;
    w.rho =  rho0 + rho1 * std::sin(2*M_PI*(x+y));
    w.p = p_amb;
    w.v = {1,1,0};
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //This function is ignored in Pure Hydro
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
    double L1 = 0;
    double L2 = 0;
    double Linf = 0;
    
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            double x = (i + 0.5)/n;
            double y = (j + 0.5)/n;
            double rho_exact = rho0 + rho1 * std::sin(2.0 * M_PI * (x + y - 2 * t));

            double err = std::abs(grid[i,j].rho - rho_exact);
            if(err > Linf) Linf = err;
            L1 += err ;
            L2 += err*err;
            
        }
    }
    L2 = std::sqrt(L2 / (n*n));
    
    std::cout<<"L1 error: "<<L1 / (n*n)<<"\n";
    std::cout<<"L2 error: "<<L2<<"\n";
    std::cout<<"Max error: "<<Linf<<"\n";

}
