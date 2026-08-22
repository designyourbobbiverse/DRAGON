//
//  Problem3D.cpp
//  DRAGON/Examples/MHD/Waves
//
//  Created by Bobbie Markwick on 09/07/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

#include "Constants.h" //For gamma, sq4pi
#include <cmath> //For std::sin
#include "DragonHoard.hpp" //For loading initial data to compare to final
#include <iostream>        //For error output

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here


typedef enum {
    Fast,
    Alfven,
    Slow,
    Entropy
} Mode;

constexpr Mode mode = Fast;
constexpr bool minus_mode = false;
constexpr int n = 64;

constexpr double epsilon = 1e-6;

const double sin_alpha = 2.0/3.0;
const double sin_beta = 2.0/std::sqrt(5);
const double cos_alpha = std::sqrt(1-sin_alpha*sin_alpha);
const double cos_beta = std::sqrt(1-sin_beta*sin_beta);
const vec3 e_x1 = vec3{cos_alpha*cos_beta, cos_alpha*sin_beta, sin_alpha};
const vec3 e_x2 = vec3{-sin_beta, cos_beta, 0};
const vec3 e_x3 = vec3{-sin_alpha*cos_beta, -sin_alpha*sin_beta, cos_alpha};


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(2*n,n, n, 1.5/n, 1.5/n, 1.5/n);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2, dy/2, dz/2) corresponds to the [0,0,0] cell
    PrimitiveState w;
    
    //Background State
    w.rho = 1.0;
    w.p = 1.0/gamma;
    w.v = (mode == Entropy ? 1.0 : 0.0) * e_x1;
    w.B = (e_x1 + 1.5 * e_x2) * sq4pi;  //Initialize magnetic field so we can perturb in conservative values
    
    //Perturbation Eigenvector
    ConservativeState R;
    switch(mode){
    case Fast:
        R.rho = 1.0 / std::sqrt(5);
        R.mom = (minus_mode ? -1 : 1) * R.rho * (2.0 * e_x1 - 1.0 * e_x2) ;
        R.B = 2 * R.rho * e_x2 * sq4pi;
        R.E = 4.5 * R.rho;
        break;
    case Alfven:
        R.mom = (minus_mode ? 1 : -1) * e_x3;
        R.B = e_x3 * sq4pi;
        break;
    case Slow:
        R.rho = 2.0 / std::sqrt(5);
        R.mom = (minus_mode ? -1 : 1) * R.rho * (0.5 * e_x1 + 1.0 * e_x2);
        R.B = -0.5 * R.rho * e_x2 * sq4pi;
        R.E = 0.75 * R.rho;
        break;
    case Entropy:
        R.rho = 1;
        R.mom = e_x1;
        R.E = 0.5;
        break;
    }
    //Apply Perturbation
    double x1 = x * cos_alpha*cos_beta + y * cos_alpha*sin_beta + z * sin_alpha;
    w += (epsilon * R) * std::cos(2*M_PI*x1);
    
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //Initialize the vector potential at point (x,y,z)
        // (0,0,0) corresponds to the [0,0,0] cell.
    //The Magnetic Fields will be initialized from this potential to ensure div B = 0
    double x1 = x * cos_alpha*cos_beta + y * cos_alpha*sin_beta + z * sin_alpha;
    double x2 = y* cos_beta - x * sin_beta;

    //Perturbation
    vec3 perturb  = {0,0,0};
    switch(mode){
        case Fast: perturb.y = 2.0 * epsilon / std::sqrt(5); break;
        case Alfven: perturb.z = epsilon; break;
        case Slow: perturb.y = -epsilon / std::sqrt(5); break;
    }
    
    //Vector Potential
    double sinx = std::sin(2*M_PI*x1) / (2*M_PI);
    double Ax2 = perturb.z * sinx;
    double Ax3 = x2 - 1.5 * x1 - perturb.y * sinx;
    return (Ax2*e_x2 + Ax3*e_x3) * sq4pi;
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


vec3 abs(vec3 v){ return {std::abs(v.x), std::abs(v.y), std::abs(v.z)}; }
void Problem::problemComplete(Grid& problem, double t){
    //This is called only after the final time is reached.
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    MyGrid& expected =  *dynamic_cast<MyGrid*>(&makeProblem());
    double ig; int nore;
    DRAGONHOARD::loadFromFile(expected, ig,nore, DRAGONHOARD::output_base_name + "_" + DRAGONHOARD::cycle_string(0));
    
    ConservativeState L1;
    
    for(int i=0; i<2*n;i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                ConservativeState err = grid[i,j,k] - expected[i,j,k];
                err.rho = std::abs(err.rho);
                err.mom = abs(err.mom);
                err.B = abs(err.B);
                err.E = std::abs(err.E);
                
                L1 += err ;
            }
        }
    }
    
    L1 /= (2*n*n*n);
    double total_err = std::sqrt(L1.rho*L1.rho + L1.mom*L1.mom + L1.B*L1.B + L1.E*L1.E);

    std::cout<<"L1 error: "<< total_err <<"\n";
}
