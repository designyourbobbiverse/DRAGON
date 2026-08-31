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
    if (output) std::cout << "\t- Zero Advection: ";
    verify_passive_scalar_zero_flux_1D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Uniform Advection: ";
    verify_passive_scalar_uniform_1D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Conservation Test: ";
    verify_passive_scalar_conserve_1D();
    if (output) std::cout << "Passed\n";

    if (output) std::cout << "- 2D:\n";
    if (output) std::cout << "\t- Adding Scalars: ";
    verify_passive_scalar_add_2D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Removing Scalars: ";
    verify_passive_scalar_remove_2D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Zero Advection: ";
    verify_passive_scalar_zero_flux_2D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Uniform Advection: ";
    verify_passive_scalar_uniform_2D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Conservation Test: ";
    verify_passive_scalar_conserve_2D();
    if (output) std::cout << "Passed\n";

    if (output) std::cout << "- 3D:\n";
    if (output) std::cout << "\t- Adding Scalars: ";
    verify_passive_scalar_add_3D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Removing Scalars: ";
    verify_passive_scalar_remove_3D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Zero Advection: ";
    verify_passive_scalar_zero_flux_3D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Uniform Advection: ";
    verify_passive_scalar_uniform_3D();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Conservation Test: ";
    verify_passive_scalar_conserve_3D();
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

//MARK: Zero Flux Advection Test
void DRAGON_Test::verify_passive_scalar_zero_flux_1D(){
    PassiveArray1D q{3,2};
    ExtendedArray1D<ConservativeState> F_X{3,2};
    FluidArray1D w{3,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        q[i,"One"]   = (rand()%10000000) * 1e-6;
        q[i,"Two"]   = (rand()%10000000) * 1e-6;
        q[i,"Three"] = (rand()%10000000) * 1e-6;
        
        w[i].rho = 1.0;
        //Make fluxes small, but might be slightly positive or negative
        F_X[i].rho = 1e-20 * (rand()%3 - 1);
    }
    auto q_ = q.advected(F_X, w, w, 1.0);
    
    for(int i=-2; i<3+2; i++){
        auto p = q[i], p_ = (*q_)[i];
        for(int n=0; n<3;n++){
            assert(p[n] = p_[n]);
        }
    }
}
void DRAGON_Test::verify_passive_scalar_zero_flux_2D(){
    PassiveArray2D q{3,4,2};
    FluxArray2D F_X{3,4,2}, F_Y{3,4,2};
    FluidArray2D w{3,4,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            q[i,j,"One"]   = (rand()%10000000) * 1e-6;
            q[i,j,"Two"]   = (rand()%10000000) * 1e-6;
            q[i,j,"Three"] = (rand()%10000000) * 1e-6;
            
            w[i,j].rho = 1.0;
            //Make fluxes small, but might be slightly positive or negative
            F_X[i,j].rho = 1e-20 * (rand()%3 - 1);
            F_Y[i,j].rho = 1e-20 * (rand()%3 - 1);
        }
    }
    auto q_ = q.advected(F_X, F_Y, w, w, 1.0, 1.0);
    
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            auto p = q[i,j], p_ = (*q_)[i,j];
            for(int n=0; n<3;n++){
                assert(p[n] = p_[n]);
            }
        }
    }
}
void DRAGON_Test::verify_passive_scalar_zero_flux_3D(){
    PassiveArray3D q{3,4,5,2};
    FluxArray3D F_X{3,4,5,2}, F_Y{3,4,5,2}, F_Z{3,4,5,2};
    FluidArray3D w{3,4,5,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                q[i,j,k,"One"]   = (rand()%10000000) * 1e-6;
                q[i,j,k,"Two"]   = (rand()%10000000) * 1e-6;
                q[i,j,k,"Three"] = (rand()%10000000) * 1e-6;
                
                w[i,j,k].rho = 1.0;
                //Make fluxes small, but might be slightly positive or negative
                F_X[i,j,k].rho = 1e-20 * (rand()%3 - 1);
                F_Y[i,j,k].rho = 1e-20 * (rand()%3 - 1);
                F_Z[i,j,k].rho = 1e-20 * (rand()%3 - 1);
            }
        }
    }
    auto q_ = q.advected(F_X, F_Y, F_Z, w, w, 1.0, 1.0, 1.0);
    
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                auto p = q[i,j,k], p_ = (*q_)[i,j,k];
                for(int n=0; n<3;n++){
                    assert(p[n] = p_[n]);
                }
            }
        }
    }
}

//MARK: Uniform Advection Test
void DRAGON_Test::verify_passive_scalar_uniform_1D(){
    PassiveArray1D q{3,2};
    ExtendedArray1D<ConservativeState> F_X{3,2};
    FluidArray1D w{3,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        q[i,"One"]   = 1.0;
        q[i,"Two"]   = 1.0;
        q[i,"Three"] = 1.0;
        
        w[i].rho = 1.0;

        F_X[i].rho = 1.0;
    }
    auto q_ = q.advected(F_X, w, w, 1.0);
    //Uniform test isn't valid on the outermost ghost
    for(int i=-1; i<3+1; i++){
        auto p = q[i], p_ = (*q_)[i];
        for(int n=0; n<3;n++){
            assert(p[n] = p_[n]);
        }
    }
}
void DRAGON_Test::verify_passive_scalar_uniform_2D(){
    PassiveArray2D q{3,4,2};
    FluxArray2D F_X{3,4,2}, F_Y{3,4,2};
    FluidArray2D w{3,4,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            q[i,j,"One"]   = 1.0;
            q[i,j,"Two"]   = 2.0;
            q[i,j,"Three"] = 3.0;
            
            w[i,j].rho = 1.0;

            F_X[i,j].rho = 0.5;
            F_Y[i,j].rho = 0.25;
        }
    }
    auto q_ = q.advected(F_X, F_Y, w, w, 1.0, 1.0);
    //Uniform test isn't valid on the outermost ghost
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            auto p = q[i,j], p_ = (*q_)[i,j];
            for(int n=0; n<3;n++){
                assert(p[n] = p_[n]);
            }
        }
    }
}
void DRAGON_Test::verify_passive_scalar_uniform_3D(){
    PassiveArray3D q{3,4,5,2};
    FluxArray3D F_X{3,4,5,2}, F_Y{3,4,5,2}, F_Z{3,4,5,2};
    FluidArray3D w{3,4,5,2};

    q.add("One");
    q.add("Two");
    q.add("Three");
    for(int i=-2; i<3+2; i++){
        for(int j=-2; j<4+2; j++){
            for(int k=-2; k<5+2; k++){
                q[i,j,k,"One"]   = 1.0;
                q[i,j,k,"Two"]   = 2.0;
                q[i,j,k,"Three"] = 3.0;
                
                w[i,j,k].rho = 1.0;

                F_X[i,j,k].rho = 0.5;
                F_Y[i,j,k].rho = 0.25;
                F_Z[i,j,k].rho = 0.125;
            }
        }
    }
    auto q_ = q.advected(F_X, F_Y, F_Z, w, w, 1.0, 1.0, 1.0);
    //Uniform test isn't valid on the outermost ghost
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            for(int k=-1; k<5+1; k++){
                auto p = q[i,j,k], p_ = (*q_)[i,j,k];
                for(int n=0; n<3;n++){
                    assert(p[n] = p_[n]);
                }
            }
        }
    }
}
//MARK: Conservation Test
void DRAGON_Test::verify_passive_scalar_conserve_1D(){
    PassiveArray1D q{3,1};
    ExtendedArray1D<ConservativeState> F_X{4,1}, F_Y{3,1};
    FluidArray1D w{3,1}, w2{3,1};

    q.add("One");
    q.add("Two");
    for(int i=-1; i<3+1; i++){
        q[i,"One"]   =  (i== 1) ? 1 : 0;
        q[i,"Two"]   = (i== 2) ? 1 : 0;
        
        w[i].rho = (rand()%10000) * 1e-2;
        w2[i] = w[i];
        
        F_X[i].rho = 1e-6 * (rand()%1000);
    }
    for(int i=-1; i<3+1; i++){
        w2[i].rho += F_X[i+1].rho - F_X[i].rho;
    }
    auto q_ = q.advected(F_X, w, w2, 1.0);

    double m1 = 0, m2 = 0, m1_ = 0, m2_ = 0;
    for(int i=-1; i<3+1; i++){
        auto p = q[i], p_ = (*q_)[i];
        
        m1 += p[0] * w[i].rho;
        m1_ += p_[0] * w2[i].rho;
        m2 += p[1] * w[i].rho;
        m2_ += p_[1] * w2[i].rho;
    }
    //Assert conservation
    assert(m1 == m1_);
    assert(m2 == m2_);
}
void DRAGON_Test::verify_passive_scalar_conserve_2D(){
    PassiveArray2D q{3,4,1};
    FluxArray2D F_X{4,4,1}, F_Y{3,5,1};
    FluidArray2D w{3,4,1}, w2{3,4,1};

    q.add("One");
    q.add("Two");
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            q[i,j,"One"]   =  (i== 1 && j == 1) ? 1 : 0;
            q[i,j,"Two"]   = (i== 2 && j == 2) ? 1 : 0;
            
            w[i,j].rho = (rand()%10000) * 1e-2;
            w2[i,j] = w[i,j];
            
            F_X[i,j].rho = 1e-6 * (rand()%1000);
            F_Y[i,j].rho = 1e-6 * (rand()%1000);
        }
    }
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            w2[i,j].rho += 1.0*(F_X[i+1,j].rho - F_X[i,j].rho);
            w2[i,j].rho += 0.5*(F_Y[i,j+1].rho - F_Y[i,j].rho);
        }
    }
    auto q_ = q.advected(F_X, F_Y, w, w2, 1.0, 0.5);

    double m1 = 0, m2 = 0, m1_ = 0, m2_ = 0;
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            auto p = q[i,j], p_ = (*q_)[i,j];
            
            m1 += p[0] * w[i,j].rho;
            m1_ += p_[0] * w2[i,j].rho;
            m2 += p[1] * w[i,j].rho;
            m2_ += p_[1] * w2[i,j].rho;
        }
    }
    //Assert conservation
    assert(m1 == m1_);
    assert(m2 == m2_);
}
void DRAGON_Test::verify_passive_scalar_conserve_3D(){
    PassiveArray3D q{3,4,5,1};
    FluxArray3D F_X{4,4,5,1}, F_Y{3,5,5,1}, F_Z{3,4,6,1};
    FluidArray3D w{3,4,5,1}, w2{3,4,5,1};

    q.add("One");
    q.add("Two");
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            for(int k=-1; k<5+1; k++){
                q[i,j,k,"One"]   =  (i== 1 && j == 1 && k == 1) ? 1 : 0;
                q[i,j,k,"Two"]   = (i== 2 && j == 2 && k == 2) ? 1 : 0;
                
                w[i,j,k].rho = (rand()%10000) * 1e-2;
                w2[i,j,k] = w[i,j,k];
                
                F_X[i,j,k].rho = 1e-6 * (rand()%1000);
                F_Y[i,j,k].rho = 1e-6 * (rand()%1000);
                F_Z[i,j,k].rho = 1e-6 * (rand()%1000);
            }
        }
    }
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            for(int k=-1; k<5+1; k++){
                w2[i,j,k].rho += 1.0*(F_X[i+1,j,k].rho - F_X[i,j,k].rho);
                w2[i,j,k].rho += 0.5*(F_Y[i,j+1,k].rho - F_Y[i,j,k].rho);
                w2[i,j,k].rho += 0.2*(F_Z[i,j,k+1].rho - F_Z[i,j,k].rho);
            }
        }
    }
    auto q_ = q.advected(F_X, F_Y, F_Z, w, w2, 1.0, 0.5, 0.2);

    double m1 = 0, m2 = 0, m1_ = 0, m2_ = 0;
    for(int i=-1; i<3+1; i++){
        for(int j=-1; j<4+1; j++){
            for(int k=-1; k<5+1; k++){
                auto p = q[i,j,k], p_ = (*q_)[i,j,k];
                
                m1 += p[0] * w[i,j,k].rho;
                m1_ += p_[0] * w2[i,j,k].rho;
                m2 += p[1] * w[i,j,k].rho;
                m2_ += p_[1] * w2[i,j,k].rho;

            }
        }
    }
    //Assert conservation
    assert(m1 == m1_);
    assert(m2 == m2_);
}

