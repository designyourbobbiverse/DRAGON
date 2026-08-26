//
//  Problem.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

typedef Grid1D MyGrid;//Replace Grid1D with one of the following:
    //Grid1D: single dimension, single thread
    //Grid2D: two dimensions, single thread
    //Grid3D: three dimensions, single thread
    //DistGrid1D: single dimension, multiple threads
    //DistGrid2D: two dimensions, multiple threads
    //DistGrid3D: three dimensions, multiple threads

Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto _grid = new MyGrid(64, 1.0);
    return *_grid;
}

PrimitiveState Problem::initialFluidState(double x, double y, double z){
    PrimitiveState w;
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2,dz/2) corresponds to the [0,0,0] cell
        //In 1D and 2D, z will always be zero, as will y in 1D

    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    vec3 A;
    //Initialize the vector potential at point (x,y,z)
        // (0,0,0)corresponds to the [0,0,0] cell. In 2D, z will always be zero
    //In 2D & 3D, Magnetic Fields will be initialized from this potential to ensure div B = 0
    //This function is ignored in pure Hydro and in 1D (1D B should be given in initialFluidState)
    
    return A;
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
}
