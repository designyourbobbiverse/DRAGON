//  DRAGON/main/main.cpp
//  Created by Bobbie Markwick
//

#include "Problem.hpp" //The user's problem

#include "DragonHoard.hpp"  //File output
#include <iostream>         //Console messages
#include <chrono>           //Runtime measurement
using namespace DRAGON;

#ifndef TESTMODE
static void verify_dir();
static void printCycleComplete(std::string cycleStr, double clock_time);


int main(int argc, const char * argv[]) {
    verify_dir();
    
    auto start = std::chrono::system_clock::now();
    Grid& problem = Problem::makeProblem();
    double time = 0.0;
    int cycle = 0;
    Problem::load(problem, time, cycle);
    
    //Monitor Output
    std::string cycleStr = DRAGONHOARD::cycle_string(cycle);
    double clock_time = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
    printCycleComplete(cycleStr, clock_time);
    
    
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
        printCycleComplete(cycleStr, clock_time);
        
        //Write to File
        DRAGONHOARD::writeToFile(problem, time, cycle, DRAGONHOARD::output_base_name + "_" + cycleStr);
    }
    
    Problem::problemComplete(problem, time);
    
    return 0;
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
static void printCycleComplete(std::string cycleStr, double clock_time){
    //Cycle Number
    std::cout<<"Frame "<<cycleStr<<" computed, ";
    //Time
    int h = floor(clock_time/3600.0), m = floor((clock_time-h*3600.0)/60.0);
    double s = round((clock_time - h*3600 - m*60)*100)/100.0;
    std::cout << "Time: "<< h << "h "<< m <<"m " << s <<"s \n";
}



#endif
