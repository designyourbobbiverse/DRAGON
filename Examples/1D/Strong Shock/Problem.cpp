//
//  Problem_SOD.cpp
//  DRAGON/Example/1D/Strong Shock
//
//  Created by Bobbie Markwick on 7/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"

typedef DistGrid1D MyGrid;//Choose the dimension of your grid here

constexpr double rho_L = 1.0;
constexpr double rho_R = 1.0;
constexpr double p_L = 1000.0;
constexpr double p_R = 0.01;
const int nx = 1024*8;

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
