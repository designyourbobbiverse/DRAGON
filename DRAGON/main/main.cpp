//  DRAGON/main/main.cpp
//  Created by Bobbie Markwick
//

#include "Problem.hpp"
#include "Config.h"
#include "DragonHoard.hpp"
#include <iostream>
#include <chrono>

#ifndef TESTMODE


static void verify_dir(){
    try{
        DRAGONHOARD::verifyOutputDirectory();
    } catch (std::exception& e){
        std::cout<<e.what()<<std::endl;
        throw e;
    }
}
static void load(Grid& problem, double& time, int& cycle){
    #ifdef RESTART_FROM_FILE
    std::string file = DRAGONHOARD::restartFileName();
    if(file.size() > 0) {
        DRAGONHOARD::loadFromFile(problem, time, cycle, file);
    } else {
        Problem::initializeProblem(problem);
        DRAGONHOARD::writeToFile(problem, 0, 0, CONFIG::output_base_name + "_" + DRAGONHOARD::cycle_string(0));
    }
    #else
    Problem::initializeProblem(problem);
    IO::writeToFile(problem, 0, 0, CONFIG::output_base_name + "_" + IO::cycle_string(0));
    #endif
}
 

static void cycle_output(std::string cycleStr, double clock_time){
    //Cycle Number
    std::cout<<"Frame "<<cycleStr<<" computed, ";
    //Time
    int h = floor(clock_time/3600.0), m = floor((clock_time-h*3600.0)/60.0);
    double s = round((clock_time - h*3600 - m*60)*100)/100.0;
    std::cout << "Time: "<< h << "h "<< m <<"m " << s <<"s \n";
}


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
        DRAGONHOARD::writeToFile(problem, time, cycle, CONFIG::output_base_name + "_" + cycleStr);
    }
    
    Problem::problemComplete(problem, time);
    
    return 0;
}


#endif
