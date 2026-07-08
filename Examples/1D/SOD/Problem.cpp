//
//  Problem.cpp
//  DRAGON/Examples/1D/SOD
//
//  Created by Bobbie Markwick on 7/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"

typedef Grid1D MyGrid;//Choose the dimension of your grid here

constexpr double rho_L = 1.0;
constexpr double rho_R = 0.125;
constexpr double p_L = 1.0;
constexpr double p_R = 0.1;
const int nx = 1024;

Grid& Problem::makeProblem(){
    auto grid = new MyGrid(nx, 1.0/nx);
    grid->boundary = Boundary::Outflow();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    for(int i = 0; i<nx; i++){
        grid[i].rho = i < nx/2 ? rho_L : rho_R;
        grid[i].p = i < nx/2 ? p_L : p_R;
        grid[i].v = {0,0,0};
    }
}


void Problem::cycleComplete(Grid& problem, int cycle){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the file write, and is called for the initial frame as well as subsequent updates
    
}
