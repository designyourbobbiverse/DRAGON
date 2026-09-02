//
//  WING_Memory_Tests.cpp
//  DRAGON/Testing/Infrastructure/WING
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "DragonWing.hpp"
using namespace DRAGON_Test;

//Things not tested:
//REUSE_AUX_GRIDS disabled
//Concurrency Stress Testing


//MARK: Memory Allocation
void DRAGON_Test::verify_WING_alloc_prim_dimensions(){
    const int N = 3, nx = 4, ny = 5, nz = 6, g = 2;
    //1D
    auto arrs1D = DRAGONWING::requestPrimitiveArrays(N, nx, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs1D[i] != nullptr);
        ExtendedArray1D<PrimitiveState>& arr = *arrs1D[i];
        //Array has the right dimensions
        assert(arr.getSize() == nx);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs1D[i] != arrs1D[j]);
        }
    }
    //2D
    auto arrs2D = DRAGONWING::requestPrimitiveArrays(N, nx, ny, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs2D[i] != nullptr);
        ExtendedArray2D<PrimitiveState>& arr = *arrs2D[i];
        //Array has the right dimensions
        assert(arr.getSizeX() == nx);
        assert(arr.getSizeY() == ny);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs2D[i] != arrs2D[j]);
        }
    }
    //3D
    auto arrs3D = DRAGONWING::requestPrimitiveArrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs3D[i] != nullptr);
        ExtendedArray3D<PrimitiveState>& arr = *arrs3D[i];
        //Array has the right dimensions
        assert(arr.getSizeX() == nx);
        assert(arr.getSizeY() == ny);
        assert(arr.getSizeZ() == nz);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs3D[i] != arrs3D[j]);
        }
    }
}

void DRAGON_Test::verify_WING_alloc_flux_dimensions(){
    const int N = 2, nx = 4, ny = 5, nz = 6, g = 2;
    //1D
    auto arrs1D = DRAGONWING::requestFluxArrays(N, nx, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs1D[i] != nullptr);
        ExtendedArray1D<ConservativeState>& arr = *arrs1D[i];
        //Array has the right dimensions
        assert(arr.getSize() == nx);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs1D[i] != arrs1D[j]);
        }
    }
    //2D
    auto arrs2D = DRAGONWING::requestFluxArrays(N, nx, ny, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs2D[i] != nullptr);
        ExtendedArray2D<ConservativeState>& arr = *arrs2D[i];
        //Array has the right dimensions
        assert(arr.getSizeX() == nx);
        assert(arr.getSizeY() == ny);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs2D[i] != arrs2D[j]);
        }
    }
    //3D
    auto arrs3D = DRAGONWING::requestFluxArrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs3D[i] != nullptr);
        ExtendedArray3D<ConservativeState>& arr = *arrs3D[i];
        //Array has the right dimensions
        assert(arr.getSizeX() == nx);
        assert(arr.getSizeY() == ny);
        assert(arr.getSizeZ() == nz);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs3D[i] != arrs3D[j]);
        }
    }
}
void DRAGON_Test::verify_WING_alloc_vec3_dimensions(){
    const int N = 4, nx = 7, ny = 6, nz = 5, g = 3;
    //1D
    auto arrs1D = DRAGONWING::requestVec3Arrays(N, nx, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs1D[i] != nullptr);
        ExtendedArray1D<vec3>& arr = *arrs1D[i];
        //Array has the right dimensions
        assert(arr.getSize() == nx);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs1D[i] != arrs1D[j]);
        }
    }
    //2D
    auto arrs2D = DRAGONWING::requestVec3Arrays(N, nx, ny, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs2D[i] != nullptr);
        ExtendedArray2D<vec3>& arr = *arrs2D[i];
        //Array has the right dimensions
        assert(arr.getSizeX() == nx);
        assert(arr.getSizeY() == ny);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs2D[i] != arrs2D[j]);
        }
    }
    //3D
    auto arrs3D = DRAGONWING::requestVec3Arrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        //Array exists
        assert(arrs3D[i] != nullptr);
        ExtendedArray3D<vec3>& arr = *arrs3D[i];
        //Array has the right dimensions
        assert(arr.getSizeX() == nx);
        assert(arr.getSizeY() == ny);
        assert(arr.getSizeZ() == nz);
        assert(arr.getGhosts() == g);
        //Arrays are unique
        for(int j=0; j<N; j++){
            if(i!=j) assert(arrs3D[i] != arrs3D[j]);
        }
    }
}
//MARK: Array Reuse
void DRAGON_Test::verify_WING_array_reuse(){
     //Test is meaningless if the user has disabled array reuse
    const int N = 3, nx = 4, ny = 5, nz = 6, g = 2;
    auto arrs1D = DRAGONWING::requestFluxArrays(N, nx, g);
    auto arrs2D = DRAGONWING::requestVec3Arrays(N, nx, ny, g);
    auto arrs3D = DRAGONWING::requestPrimitiveArrays(N, nx, ny, nz, g);
    //Get references to the arrays
    ExtendedArray1D<ConservativeState>* a1D[N];
    for(int i=0; i<N; i++) a1D[i] = arrs1D[i];
    ExtendedArray2D<vec3>* a2D[N];
    for(int i=0; i<N; i++) a2D[i] = arrs2D[i];
    ExtendedArray3D<PrimitiveState>* a3D[N];
    for(int i=0; i<N; i++) a3D[i] = arrs3D[i];

    //Check that a second batch doesn't get the same arrays
    auto arrs1D_2 = DRAGONWING::requestFluxArrays(N, nx+1, g);
    auto arrs2D_2 = DRAGONWING::requestVec3Arrays(N, nx, ny, g);
    auto arrs3D_2 = DRAGONWING::requestPrimitiveArrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            assert(a1D[i] != arrs1D_2[j]);
            assert(a2D[i] != arrs2D_2[j]);
            assert(a3D[i] != arrs3D_2[j]);
        }
    }
    //Release of one batch doesn't release the others
    arrs2D.release();
    auto arrs1D_3 = DRAGONWING::requestFluxArrays(N, nx, g);
    auto arrs2D_3 = DRAGONWING::requestVec3Arrays(N, nx, ny, g);
    auto arrs3D_3 = DRAGONWING::requestPrimitiveArrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        bool reused = false;
        for(int j=0; j<N; j++){
            assert(a1D[i] != arrs1D_3[j]);
            reused |= (a2D[i] == arrs2D_3[j]);
            assert(a3D[i] != arrs3D_3[j]);
        }
        assert(reused);
    }
    //Requesting array of different dimension is irrelevant
    arrs1D.release();
    arrs3D.release();
    auto arrs1D_4 = DRAGONWING::requestFluxArrays(2*N, nx+1, g);
    auto arrs2D_4 = DRAGONWING::requestVec3Arrays(2*N, nx+1, ny, g);
    auto arrs3D_4 = DRAGONWING::requestPrimitiveArrays(2*N, nx+1, ny, nz, g);
    for(int i=0; i<N; i++){
        for(int j=0; j<2*N; j++){
            assert(a1D[i] != arrs1D_4[j]);
            assert(a2D[i] != arrs2D_4[j]);
            assert(a3D[i] != arrs3D_4[j]);
        }
    }
    //Release & Reuse of the others
    auto arrs1D_5 = DRAGONWING::requestFluxArrays(2*N, nx, g);
    auto arrs2D_5 = DRAGONWING::requestVec3Arrays(2*N, nx, ny, g);
    auto arrs3D_5 = DRAGONWING::requestPrimitiveArrays(2*N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        bool reused1D = false, reused3D = false;
        for(int j=0; j<2*N; j++){
            reused1D |= (a1D[i] == arrs1D_5[j]);
            assert(a2D[i] != arrs2D_5[j]);
            reused3D |= (a3D[i] == arrs3D_5[j]);
        }
        assert(reused1D && reused3D);
    }
}


//MARK: ArrayGuard
void DRAGON_Test::verify_WING_arrayguard_get(){
    const int N = 3, nx = 4, ny = 5, nz = 6, g = 2;
    auto arrs1D = DRAGONWING::requestPrimitiveArrays(N, nx, g);
    auto arrs2D = DRAGONWING::requestFluxArrays(N, nx, ny, g);
    auto arrs3D = DRAGONWING::requestVec3Arrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        assert(arrs1D[i] == arrs1D.get(i));
        assert(arrs2D[i] == arrs2D.get(i));
        assert(arrs3D[i] == arrs3D.get(i));
    }
}
void DRAGON_Test::verify_WING_scope_release(){
    const int N = 2, nx = 4, ny = 5, nz = 6, g = 2;
    ExtendedArray3D<PrimitiveState>* arr;
    {
        auto arrs = DRAGONWING::requestPrimitiveArrays(N, nx,ny,nz, g);
        arr = arrs[0];
    }
    auto arrs2 = DRAGONWING::requestPrimitiveArrays(N, nx,ny,nz, g);
    assert(arrs2[0] == arr);
}

void DRAGON_Test::verify_WING_double_release(){
    const int N = 3, nx = 4, ny = 5, nz = 6, g = 2;
    auto arrs = DRAGONWING::requestPrimitiveArrays(N, nx,ny,nz, g);
    auto* arr = arrs[0];
    arrs.release();
    assert(arr->getSizeX()>0);//This verifies that arr hasn't been deleted yet
    arrs.release(); //Shouldn't actually do anything
    assert(arr->getSizeY()>0);//This verifies that arr hasn't been deleted yet
}

//MARK: Purge
//MARK: Array Reuse
void DRAGON_Test::verify_WING_array_purge_survival(){
     //Test is meaningless if the user has disabled array reuse
    const int N = 3, nx = 4, ny = 5, nz = 6, g = 2;
    auto arrs1D = DRAGONWING::requestFluxArrays(N, nx, g);
    auto arrs2D = DRAGONWING::requestVec3Arrays(N, nx, ny, g);
    auto arrs3D = DRAGONWING::requestPrimitiveArrays(N, nx, ny, nz, g);
    //Get references to the arrays
    ExtendedArray1D<ConservativeState>* a1D[N];
    for(int i=0; i<N; i++) a1D[i] = arrs1D[i];
    ExtendedArray2D<vec3>* a2D[N];
    for(int i=0; i<N; i++) a2D[i] = arrs2D[i];
    ExtendedArray3D<PrimitiveState>* a3D[N];
    for(int i=0; i<N; i++) a3D[i] = arrs3D[i];

    
    //Existing arrays survive the purge
    DRAGONWING::purgeAllBuffers();
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            assert(a1D[i]->getGhosts() > 0);
            assert(a2D[i]->getGhosts() > 0);
            assert(a3D[i]->getGhosts() > 0);
        }
    }
    //Fresh arrays are allocated
    auto arrs1D_3 = DRAGONWING::requestFluxArrays(N, nx, g);
    auto arrs2D_3 = DRAGONWING::requestVec3Arrays(N, nx, ny, g);
    auto arrs3D_3 = DRAGONWING::requestPrimitiveArrays(N, nx, ny, nz, g);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            assert(a1D[i] != arrs1D_3[j]);
            assert(a2D[i] != arrs2D_3[j]);
            assert(a3D[i] != arrs3D_3[j]);
        }
    }
}


