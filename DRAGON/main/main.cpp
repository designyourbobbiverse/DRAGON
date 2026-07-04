//  DRAGON/main/main.cpp
//  Created by Bobbie Markwick
//

#include <chrono>
#include <iostream>
#include "Problem.hpp"
#include "HDF5Output.hpp"

#ifdef TESTMODE
#error Test mode was enabled, but main.cpp is compiled. Check your build settings
#endif


int main(int argc, const char * argv[]) {
    
    auto start = std::chrono::system_clock::now();

    Grid& problem = Problem::makeProblem();

    for(int n = 0; n<1000; n++){
        if(n>0) problem.advance(1);
        //Let the problem code to do any special processing
        Problem::cycleComplete(problem, n);
        
        //Ouptut Cycle Number
        auto numStr = std::string(n<1000 ? "0" : "") + std::string(n<100 ? "0" : "") + std::string(n<10 ? "0" : "") +  std::to_string(n);
        std::cout<<"Frame "<<numStr<<" computed, ";
        //Time
        std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
        double time = elapsed_seconds.count();
        int h = floor(time/3600.0), m = floor((time-h*3600.0)/60.0);
        double s = round((time - h*3600 - m*60)*100)/100.0;
        std::cout << "Time: "<< h << "h "<< m <<"m " << s <<"s \n";
        
        
        //Write to File
       // IO::writeToFile(problem, n, n, "Frame_"+numStr+".h5");
    }
    
    return 0;
}
