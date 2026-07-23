//
//  Problem2D.cpp
//  DRAGON/Examples/Orszag–Tang
//
//  Created by Bobbie Markwick on 22/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include "Constants.h"
#include <cmath>
#include <iostream>

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 25.0 / (36*M_PI);
constexpr double p_amb = 5.0 / (12*M_PI);
constexpr double B0 = 1.0;

constexpr int n = 512;
constexpr double dx = 2*M_PI/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, dx, dx);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);


    //Set up the ambient grid and calculate the size of the blast
    for(int i=0; i<=grid.getSizeX();i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            double x = (i+0.5) * dx;
            double y = (j+0.5) * dx;

            grid[i,j].rho = rho0;
            grid[i,j].v = {-sin(y),sin(x),0};
            grid[i,j].p = p_amb;
            
            x -= 0.5 * dx;
            y -= 0.5 * dx;
            
            double Az = cos(y) + 0.25*cos(2*x);
            grid._A()[i,j] = {0,0, B0*Az};
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
}
