//
//  Problem3D.cpp
//  DRAGON/Examples/Diagonal
//
//  Created by Bobbie Markwick on 09/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>
#include <iostream>

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1.0;
constexpr double rho1 = 0.1;
constexpr double p_amb = 1.0;
constexpr int n = 128;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, n, 1.0/n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);


    //Set up the ambient grid and calculate the size of the blast
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                double x = (i + 0.5)/n;
                double y = (j + 0.5)/n;
                double z = (k + 0.5)/n;

                grid[i,j,k].rho = rho0 + rho1 * sin(2*M_PI*(x+y+z));
                grid[i,j,k].v = {1,1,1};
                grid[i,j,k].p = p_amb;
            }
        }
    }
    
    
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
    //This is called only after the final time is reached.
    double L1 = 0;
    double L2 = 0;
    double Linf = 0;
    
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                double x = (i + 0.5)/n;
                double y = (j + 0.5)/n;
                double z = (k + 0.5)/n;
                double rho_exact = 1.0 + 0.1 * sin(2.0 * M_PI * (x + y + z - 3 * t));
                
                double err = fabs(grid[i,j,k].rho - rho_exact);
                if(err > Linf) Linf = err;
                L1 += err ;
                L2 += err*err;
            }
        }
    }
    L2 = sqrt(L2 / (n*n*n));
    
    std::cout<<"L1 error: "<<L1 / (n*n*n)<<"\n";
    std::cout<<"L2 error: "<<L2<<"\n";
    std::cout<<"Max error: "<<Linf<<"\n";
}
