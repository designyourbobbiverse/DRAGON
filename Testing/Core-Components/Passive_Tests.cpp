//
//  Passive_Tests.cpp
//  DRAGON/Testing/Core-Components
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Passives/Passives.hpp"

#include "Config.h"
#include <iostream> //For std::cout

using namespace DRAGON_Test;

void DRAGON_Test::verify_passive_scalars(bool output){
    if (output) std::cout << "Passive Scalars:\n";
    
    if (output) std::cout << "- Key Management:\n";
    if (output) std::cout << "\t- lookup throws if bad: ";
    verify_passive_scalar_bad_lookup();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Adding Keys: ";
    verify_passive_scalar_add();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Removing Keys: ";
    verify_passive_scalar_remove();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Key add/remove no-op: ";
    verify_passive_scalar_no_op();
    if (output) std::cout << "Passed\n";
    
    
    if (output) std::cout << "- 1D:\n";
    if (output) std::cout << "\t- Adding Scalars: ";
    verify_passive_scalar_add_1D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Removing Scalars: ";
    verify_passive_scalar_remove_1D();
    if (output) std::cout << "Passed\n";

    if (output) std::cout << "- 2D:\n";
    if (output) std::cout << "\t- Adding Scalars: ";
    verify_passive_scalar_add_2D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Removing Scalars: ";
    verify_passive_scalar_remove_2D();
    if (output) std::cout << "Passed\n";

    if (output) std::cout << "- 3D:\n";
    if (output) std::cout << "\t- Adding Scalars: ";
    verify_passive_scalar_add_3D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Removing Scalars: ";
    verify_passive_scalar_remove_3D();
    if (output) std::cout << "Passed\n";


    if (output) std::cout << "All Passive Scalar Tests Passed.\n\n";
}


//MARK: Key Management
void DRAGON_Test::verify_passive_scalar_bad_lookup(){
    PassiveArray1D q{1,0};
    PassiveSet p = {1, 3};
    try{
        double val = q.lookup(p, "Key");
        assert(val != val);
    } catch (...) {} //Test passes on throw
}

void DRAGON_Test::verify_passive_scalar_add(){
    PassiveArray1D q{1,0};
    
    PassiveSet p = {1, 2, 3};

    q.add("One");
    assert(q.count() == 1);
    assert(q.lookup(p, "One") == 1);
    
    q.add("Two");
    assert(q.lookup(p, "Two") == 2);

    q.add("Three");
    assert(q.count() == 3);
    assert(q.lookup(p, "Three") == 3);

}
void DRAGON_Test::verify_passive_scalar_remove(){
    PassiveArray1D q{1,0};
    q.add("One");
    q.add("Two");
    q.add("Three");
    assert(q.count() == 3);
    
    PassiveSet p = {1, 3};
    assert(q.lookup(p, "Two") == 3);
    
    q.remove("Two");
    assert(q.count() == 2);
    
    assert(q.lookup(p, "One") == 1);
    assert(q.lookup(p, "Three") == 3);
    try{
        double val = q.lookup(p, "Two");
        assert(val != val);
    } catch (...) {} //Test passes on throw
}
void DRAGON_Test::verify_passive_scalar_no_op(){
    PassiveArray1D q{1,0};
    q.add("Essence");
    q.add("Essence");
    assert(q.count() == 1);
    q.remove("Silence");
    assert(q.count() == 1);
}


//MARK: Adding/Removing in Arrays

void DRAGON_Test::verify_passive_scalar_add_1D(){
    PassiveArray1D q{3,2};
    
    q.add("One"); //Verify that "One" was added
    for(int i=-2; i<3+2; i++){
        auto& val = q[i,"One"];
        assert(val == 0.0);
        val = 1.0;
    }
    q.add("Two"); //Verify that adding "Two" doesn't alter the values of "One"
    for(int i=-2; i<3+2; i++){
        auto p = q[i];
        assert(p.size() == 2);
        auto v1 = q.lookup(p, "One"), v2 = q.lookup(p, "Two");
        assert(v1 == 1.0);
        assert(v2 == 0.0);
    }
    
}
void DRAGON_Test::verify_passive_scalar_remove_1D(){
    PassiveArray1D q{3,2};
    
    q.add("One"); //Test that the correct key is being removed
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++) {
        q[i,"One"] = 1.0;
        q[i,"Two"] = 2.0;
        q[i,"Three"] = 3.0;
    }
    q.remove("Two");
    for(int i=-2; i<3+2; i++){
        auto p = q[i];
        assert(p.size() == 2);
        assert(q.lookup(p, "One") == 1.0);
        assert(q.lookup(p, "Three") == 3.0);
    }
}
void DRAGON_Test::verify_passive_scalar_add_2D(){
    PassiveArray2D q{3,4,2};

    q.add("One");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            auto& val = q[i,j,"One"];
            assert(val == 0.0);
            val = 1.0;
        }
    }
    q.add("Two");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            auto p = q[i,j];
            assert(p.size() == 2);
            auto v1 = q.lookup(p, "One"), v2 = q.lookup(p, "Two");
            assert(v1 == 1.0);
            assert(v2 == 0.0);
        }
    }
}
void DRAGON_Test::verify_passive_scalar_remove_2D(){
    PassiveArray2D q{3,4,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            q[i,j,"One"]   = 1.0;
            q[i,j,"Two"]   = 2.0;
            q[i,j,"Three"] = 3.0;
        }
    }
    q.remove("Two");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            auto p = q[i,j];
            assert(p.size() == 2);
            assert(q.lookup(p, "One")   == 1.0);
            assert(q.lookup(p, "Three") == 3.0);
        }
    }
}
void DRAGON_Test::verify_passive_scalar_add_3D(){
    PassiveArray3D q{3,4,5,2};

    q.add("One");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                auto& val = q[i,j,k,"One"];
                assert(val == 0.0);
                val = 1.0;
            }
        }
    }
    q.add("Two");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                auto p = q[i,j,k];
                assert(p.size() == 2);
                auto v1 = q.lookup(p, "One"), v2 = q.lookup(p, "Two");
                assert(v1 == 1.0);
                assert(v2 == 0.0);
            }
        }
    }
}
void DRAGON_Test::verify_passive_scalar_remove_3D(){
    PassiveArray3D q{3,4,5,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                q[i,j,k,"One"]   = 1.0;
                q[i,j,k,"Two"]   = 2.0;
                q[i,j,k,"Three"] = 3.0;
            }
        }
    }
    q.remove("Two");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                auto p = q[i,j,k];
                assert(p.size() == 2);
                assert(q.lookup(p, "One")   == 1.0);
                assert(q.lookup(p, "Three") == 3.0);
            }
        }
    }
}
