//
//  WING_Concurrency_Tests.cpp
//  DRAGON/Testing/Infrastructure/WING
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "DragonWing.hpp"

#include "Hydro/Grid.hpp" //To make a dummy grid

#include <vector>   //Execution order tests
#include <mutex>    //To make the execution order tests threadsafe

#include <stdexcept>      //To test the error handling
#include <iostream>       //To print the restart message
#include <unistd.h> //For std::sleep
using namespace DRAGON_Test;

//Things not tested:
//Any sort of concurrency stress testing
//Behaviour of Multiple calls to restart message (this one is more likely to get user overwrites, so I'd rather not force them to worry about breaking it).

thread_local int thread_num = 2112; //To verify that TestGrid is called on a separate thread

static std::mutex m{};
static std::vector<int> stages{};


class SingleTestGrid: public Grid {
public:
    bool didCall = false; //To verify that TestGrid is called
    int callThread = 0; //To verify that TestGrid is called on a separate thread
    double dt = 0;  //To verify that TestGrid is called with correct dt
    
    bool shouldThrow = false; //Set to true to test launchParallel's catching capability
    
    void split_step(double dt) override{ unsplit_step(dt); }
    void unsplit_step(double dt) override {
        didCall = true;
        callThread = thread_num;
        this->dt = dt;
        
        if(shouldThrow){
            throw std::runtime_error("Passed");
        } else {
            DRAGONWING::reportCheckpoint2();
        }
        
        
    }
    
};
class OrderTestGrid: public Grid {
public:
    bool shouldThrow = false;

    
    void split_step(double dt) override{ unsplit_step(dt); }
    void unsplit_step(double dt) override {
        if(!DRAGONWING::waitForRelease()) return;
        m.lock();
        stages.push_back(1);
        m.unlock();
        usleep(2);
        if(shouldThrow) throw std::runtime_error("Oopsie");
        m.lock();
        stages.push_back(-1);
        m.unlock();
        DRAGONWING::reportCheckpoint1();
        if(!DRAGONWING::waitForCheckpoint1()) return;
        m.lock();
        stages.push_back(2);
        m.unlock();
        usleep(1);
        m.lock();
        stages.push_back(-2);
        m.unlock();
        DRAGONWING::reportCheckpoint2();
        
    }
    
};


//MARK: Launch Parallel
void DRAGON_Test::verify_WING_launchParallel_calls_advance(){
    SingleTestGrid grid{};
    thread_num++;//Make this thread look different from the others
    
    auto tp = DRAGONWING::ThreadPool(1);
    tp.launchParallel(&grid, 3.14);
    bool succ = tp.waitForCompletion();
    
    assert(grid.didCall);
    assert(grid.callThread != thread_num);
    assert(grid.dt == 3.14);
    assert(succ);
}
void DRAGON_Test::verify_WING_launchParallel_catches(){
    SingleTestGrid grid{};
    grid.shouldThrow = true;
    
    auto tp = DRAGONWING::ThreadPool(1);
    tp.launchParallel(&grid, 3.14);
    bool succ = tp.waitForCompletion();
    
    assert(!succ);
    std::cout << tp.restartMsg();
}


//MARK: Restart Message
void DRAGON_Test::verify_WING_restart_message(){
    auto tp = DRAGONWING::ThreadPool(1);
    assert(tp.restartMsg().length() == 0);
    DRAGONWING::requestRestart("Hallucination"); //This shouldn't do anything
    assert(tp.restartMsg().length() == 0);
    tp.requestRestart("Passed");
    auto str = tp.restartMsg();
    assert(tp.restartMsg().length() == 0);//Clears
    std::cout<<str;
}

//MARK: Checkpoints
void DRAGON_Test::verify_WING_checkpoints(){
    int N = 100;
    auto tp = DRAGONWING::ThreadPool(N);
    auto* grid = new OrderTestGrid[N]();
    for(int i=0; i<N; i++){
        tp.launchParallel(&grid[i], 3.14);
    }
    bool succ = tp.waitForCompletion();
    assert(succ);
    delete [] grid;
    
    
    m.lock();
    int phase1 = 0;
    bool phase2 = false;
    int net = 0;
    int didPhase2 = 0;
    for(auto p: stages){
        net += p;
        if (p==2) didPhase2++;
        //Make sure everyone finishes phase 1 before anyone starts phase 2
        if(phase2) {
            assert(p != 1 && p != -1);
        } else if(p == 2 || p == -2) {
            phase2 = true;
        } else { //In phase 1, make sure we don't have too many
            phase1 += p;
            if (DRAGONWING::Config::phase_1_max_threads <= 0) continue;
            assert(phase1 <= DRAGONWING::Config::phase_1_max_threads);
        }
    }
    assert(didPhase2 == N); //Verify that everyone entered phase 2
    assert(net == 0); //Verify that everyone in phase 2 exited phase 2
    stages.clear();
    m.unlock();
}

void DRAGON_Test::verify_WING_checkpoints_restart(){
    int N = 100;
    auto tp = DRAGONWING::ThreadPool(N);
    auto* grid = new OrderTestGrid[N]();
    grid[69].shouldThrow = true;
    for(int i=0; i<N; i++){
        tp.launchParallel(&grid[i], 3.14);
    }
    bool succ = tp.waitForCompletion();
    assert(!succ);
    delete [] grid;
    
    
    m.lock();
    int phase1 = 0;
    int net = 0;
    for(auto p: stages){
        net += p;
        assert(p != 2 && p != -2); //Nobody should have done pahse 2
        phase1 += p; //An extra check on phase_1_max_threads doesn't hurt
        if (DRAGONWING::Config::phase_1_max_threads <= 0) continue;
        assert(phase1 <= DRAGONWING::Config::phase_1_max_threads);
    }
    //Everyone who started phase 1 finished phase 1 before the threadpool returned
    assert(net == 1); //Except the one that threw.
    
    stages.clear();
    m.unlock();
}
