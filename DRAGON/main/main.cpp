//  DRAGON/main/main.cpp
//  Created by Bobbie Markwick
//

#include "Problem.hpp" //The user's problem

#include "Config.h"
#include "DragonHoard.hpp"  //File output
#include <iostream>         //Console messages
#include <chrono>           //Runtime measurement
#include "CT.hpp"           //Magnetic Field initialization

#ifndef TESTMODE



//MARK: Problem initialization
static void init1D(Grid1D& grid){
    const int nx = grid.getSize();
    const double dx = grid.dx;
    for(int i=0; i<nx; i++){
        grid[i] = Problem::initialFluidState((i+0.5)*dx, 0, 0);
    }
    //1D MHD includes B on the cells, no need for special CT array
}
static void init2D(Grid2D& grid){
    const int nx = grid.getSizeX(), ny = grid.getSizeY();
    const double dx = grid.dx, dy = grid.dy;
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            grid[i,j] = Problem::initialFluidState((i+0.5)*dx, (j+0.5)*dy, 0);
        }
    }
    #ifdef MHD
    ExtendedArray2D<vec3> A(nx+1, ny+1, 0);
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            grid._B()[i,j] = {0,0,0};
            A[i,j] = Problem::initialMagneticPotential(i*dx, j*dy, 0);
        }
    }
    CT::Faraday(A, grid._B(), -1/dx, -1/dy, 0);
    grid.initialize_B_fields();
    #endif
}
static void init3D(Grid3D& grid){
    const int nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    const double dx = grid.dx, dy = grid.dy, dz = grid.dz;
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                grid[i,j,k] = Problem::initialFluidState((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz);
            }
        }
    }
    #ifdef MHD
    ExtendedArray3D<vec3> A(nx+1, ny+1, nz+1, 0);
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<=nz; k++){
                grid._B()[i,j,k] = {0,0,0};
                A[i,j,k] = Problem::initialMagneticPotential(i*dx, j*dy, k*dz);
            }
        }
    }
    CT::Faraday(A, grid._B(), -1/dx, -1/dy, -1/dz, 0);
    grid.initialize_B_fields();
    #endif
}
static void load(Grid& problem, double& time, int& cycle){
    #ifdef RESTART_FROM_FILE
    std::string file = DRAGONHOARD::restartFileName();
    if(file.size() > 0) {
        try{
            DRAGONHOARD::loadFromFile(problem, time, cycle, file);
            return;
        } catch(std::exception& e){
            std::cerr << e.what()<<"\nInitializing from scratch\n";
        }
    }
    #endif
    Grid1D* grid1D = dynamic_cast<Grid1D*>(&problem);
    Grid2D* grid2D = dynamic_cast<Grid2D*>(&problem);
    Grid3D* grid3D = dynamic_cast<Grid3D*>(&problem);
    if(grid1D) init1D(*grid1D);
    else if(grid2D) init2D(*grid2D);
    else if(grid3D) init3D(*grid3D);
    Problem::completeProblemInit(problem);

    DRAGONHOARD::writeToFile(problem, 0, 0, DRAGONHOARD::output_base_name + "_" + DRAGONHOARD::cycle_string(0));

}



//MARK: File I/O
static void verify_dir(){
    try{
        DRAGONHOARD::verifyOutputDirectory();
    } catch (std::exception& e){
        std::cout<<e.what()<<std::endl;
        throw e;
    }
}
 
//MARK: Cycle output
static void cycle_output(std::string cycleStr, double clock_time){
    //Cycle Number
    std::cout<<"Frame "<<cycleStr<<" computed, ";
    //Time
    int h = floor(clock_time/3600.0), m = floor((clock_time-h*3600.0)/60.0);
    double s = round((clock_time - h*3600 - m*60)*100)/100.0;
    std::cout << "Time: "<< h << "h "<< m <<"m " << s <<"s \n";
}


//MARK: main
int main(int argc, const char * argv[]) {
    verify_dir();
    
    auto start = std::chrono::system_clock::now();
    Grid& problem = Problem::makeProblem();
    double time = 0.0;
    int cycle = 0;
    
    load(problem, time, cycle);
    //Monitor Output
    std::string cycleStr = DRAGONHOARD::cycle_string(cycle);
    double clock_time = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
    cycle_output(cycleStr, clock_time);
    
    
    while(time < CONFIG::final_time){
        //Let the problem code to do any special processing
        Problem::beforeCycle(problem, cycle, time);
        //Advance the (Magneto)Hydrodynamics
        problem.advance(CONFIG::dt);
        time += CONFIG::dt;
        cycle++;
        //Let the problem code to do any special processing
        Problem::afterCycle(problem, cycle, time);
        
        //Monitor Output
        std::string cycleStr = DRAGONHOARD::cycle_string(cycle);
        double clock_time = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
        cycle_output(cycleStr, clock_time);
        
        //Write to File
        DRAGONHOARD::writeToFile(problem, time, cycle, DRAGONHOARD::output_base_name + "_" + cycleStr);
    }
    
    Problem::problemComplete(problem, time);
    
    return 0;
}


#endif
