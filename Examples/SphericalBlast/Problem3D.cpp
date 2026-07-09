//
//  Problem3D.cpp
//  DRAGON/Examples/SphericalBlast
//
//  Created by Bobbie Markwick on 8/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>
#include "Constants.h"

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho = 1.0;
constexpr double p_amb = 1e-5;
constexpr double E_blast = 1.0;
constexpr int n = 256;
constexpr double r0 = 12.0/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n,n, 1.0/n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Outflow();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);


    //Set up the ambient grid and calculate the size of the blast
    double blast_vol = 0;
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                double x = (i + 0.5)/n  - 0.5;
                double y = (j + 0.5)/n - 0.5;
                double z = (k + 0.5)/n - 0.5;
                double r = sqrt(x*x + y*y + z*z);
                
                grid[i,j,k].rho = rho;
                grid[i,j,k].v = {0,0,0};
                grid[i,j,k].p = p_amb;
                
                if(r < r0) blast_vol += 1.0/(n*n*n);
            }
        }
    }
    //Inject the blast
    const double p_blast = (_gamma - 1.0) * E_blast / blast_vol;

    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                double x = (i + 0.5)/n  - 0.5;
                double y = (j + 0.5)/n - 0.5;
                double z = (k + 0.5)/n - 0.5;
                double r = sqrt(x*x + y*y + z*z);
                
                if(r < r0) grid[i,j,k].p += p_blast ;
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
}
