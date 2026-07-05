//
//  Problem.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Problem.hpp"
#include "Refinement/DistGrid.hpp"
#include "Jets.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho_jet = 3;
constexpr double rho_amb = rho_jet/6.0;
constexpr double p_amb = rho_amb * 1.5e-5 ;
double Bz = 0.0;//sqrt(8*M_PI*p_amb)/10;

Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    int res = 2;
    double dx = 0.1/res, dy = 0.1/res, dz = 0.1/res;
    auto _grid = new MyGrid(64*res, 64*res,64*res, dx, dy, dz);
    MyGrid& grid = *_grid;
    grid.boundary = MHDJet(rho_jet, 7e-2, p_amb, 10.0, 1.0, 2.4, "Z");
    
    return grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);

    //load your initial conditions for t=0
    for(int i=0; i<=grid.getSizeX(); i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            for(int k=0; k<=grid.getSizeZ(); k++){
                grid[i,j,k].rho = rho_amb;
                grid[i,j,k].p = p_amb;
                grid[i,j,k].v = {0,0,0};
                grid._A()[i,j,k] = 0.5*Bz * vec3{-j*grid.dy,i*grid.dx,0};
            }
        }
    }
    grid.initialize_B_fields();
    
}


void Problem::cycleComplete(Grid& problem, int cycle){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the file write, and is called for the initial frame as well as subsequent updates
    
}
