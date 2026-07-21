//
//  Problem2D.cpp
//  DRAGON/Examples/LoopAdvection
//
//  Created by Bobbie Markwick on 20/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>
#include <iostream>

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1.0;
constexpr double p_amb = 5.0;
constexpr double vx = 2.0;
constexpr double vy = 0.5;
constexpr double r0 = 0.2;
constexpr double B0 = 1e-6;

constexpr int n = 256;
constexpr double dx = 2.0/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n/2, dx, dx);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);


    //Set up the ambient grid and calculate the size of the blast
    for(int i=0; i<=grid.getSizeX();i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            double x = (i) * dx - 1.0;
            double y = (j) * dx - 0.5;
            double r = sqrt(x*x + y*y);

            grid[i,j].rho = rho0;
            grid[i,j].v = {vx,vy,0};
            grid[i,j].p = p_amb;
            
            double Az = B0*fmax(r0-r, 0);
            grid._A()[i,j] = {0,0, Az};
        }
    }
    grid.initialize_B_fields();
}


void Problem::beforeCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the Godunov scheme. It isn't called for the initial frame
    problemComplete(problem, t);
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
    
    double cx = fmod(vx * t + 1.0, 2.0);
    double cy = fmod(vy * t + 0.5, 1.0);

    
    for(int i=0; i<n;i++){
        for(int j=0; j<n/2; j++){
            double x = (i+0.5) * dx - cx;
            double y = (j+0.5) * dx - cy;
            double r = sqrt(x*x + y*y);
            double B2_exact = r <= r0 ? B0*B0 : 0;

            double err = fabs(grid[i,j].B*grid[i,j].B - B2_exact) / (B0*B0);
            if(err > Linf) Linf = err;
            L1 += err ;
            L2 += err*err;
            
        }
    }
    L2 = sqrt(L2 / (n*n));
    
    std::cout<<"L1 error: "<<L1 / (n*n)<<"\n";
    std::cout<<"L2 error: "<<L2<<"\n";
    std::cout<<"L-infinity error: "<<Linf<<"\n";
}
