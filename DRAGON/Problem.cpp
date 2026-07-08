//
//  Problem.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"

typedef Grid1D MyGrid;//Choose the dimension of your grid here

Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto _grid = new MyGrid(64, 1.0);
    return *_grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    //This is where you should initialize the initial data of the grid
    
}


void Problem::cycleComplete(Grid& problem, int cycle){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the file write, and is called for the initial frame as well as subsequent updates
    
}
