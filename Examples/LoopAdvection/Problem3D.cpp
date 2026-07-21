//
//  Problem3D.cpp
//  DRAGON/Examples/LoopAdvection
//
//  Created by Bobbie Markwick on 20/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>
#include <iostream>

constexpr bool gaussian = false;

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1.0;
constexpr double p_amb = 5.0;
constexpr double vx = 2.0;
constexpr double vy = 1.0;
constexpr double vz = 1.0;
constexpr double r0 = gaussian ? 0.1 : 0.2;
constexpr double B0 = 1e-6;


constexpr int n = 256;
constexpr double dx = 2.0/n;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n/2,n/2, dx, dx, dx);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);


    //Set up the ambient grid and calculate the size of the blast
    for(int i=0; i<=grid.getSizeX();i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            for(int k=0; k<=grid.getSizeZ(); k++){
                double x = (i) * dx - 1.0;
                double y = (j) * dx - 0.5;
                double r = sqrt(x*x + y*y);
                
                grid[i,j,k].rho = rho0;
                grid[i,j,k].v = {vx,vy,vz};
                grid[i,j,k].p = p_amb;
                
                double dist = gaussian ? exp(-0.5*pow(r/r0,2)) : fmax(r0-r,0);
                grid._A()[i,j,k] = {0,0, B0 *  dist};
            }
        }
    }
    grid.initialize_B_fields();
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
    initializeProblem(expected);

    //This is called only after the final time is reached.
    double L1 = 0;
    double L2 = 0;
    double Linf = 0;

    
    for(int i=0; i<n;i++){
        for(int j=0; j<n/2; j++){
            for(int k=0; k<n/2; k++){
                double B2_expect = expected[i,j,k].B * expected[i,j,k].B;
                
                double err = fabs(grid[i,j,k].B*grid[i,j,k].B - B2_expect) / (B0*B0);
                if(err > Linf) Linf = err;
                L1 += err ;
                L2 += err*err;
            }
            
        }
    }
    L2 = sqrt(L2 / (n*n*n));
    
    std::cout<<"L1 error: "<<L1 / (n*n*n)<<"\n";
    std::cout<<"L2 error: "<<L2<<"\n";
    std::cout<<"L-infinity error: "<<Linf<<"\n";
}
