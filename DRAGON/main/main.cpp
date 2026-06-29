//  DRAGON/main/main.cpp
//  Created by Bobbie Markwick
//

#include <chrono>
#include "Problem.hpp"

#ifdef TESTMODE
#error Test mode was enabled, check your build settings
#endif


int main(int argc, const char * argv[]) {
    
    auto start = std::chrono::system_clock::now();

    Grid& problem = Problem::makeProblem();

    for(int n = 0; n<1000; n++){
        if(n>0) problem.advance(1);
        
        auto now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = now - start;
        Problem::cycleComplete(problem, n, elapsed_seconds.count());
    }
    
    return 0;
}
