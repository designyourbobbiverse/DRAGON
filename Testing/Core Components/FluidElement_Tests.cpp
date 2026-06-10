//
//  FluidElement_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Testing.hpp"
#include <math.h>
#include <iostream>

using namespace DRAGON_Test;


void DRAGON_Test::verify_fluid_element(bool output) {
    if(output) std::cout << "Fluid Arithmetic: ";
    verify_constructors();
    verify_conversion();
    verify_flux();
    verify_flux_add();
    verify_enthalpy();
    verify_swaps_P();
    verify_swaps_C();
    verify_add();
    verify_sub();
    verify_mult();
    verify_div();
    if(output) std::cout << "All tests passed.\n\n";
}


PrimitiveState DRAGON_Test::make_state(double rho, double vx, double vy, double vz, double p) {
    PrimitiveState W;
    W.rho = rho;
    W.vx = vx;
    W.vy = vy;
    W.vz = vz;
    W.p = p;
    return W;
}


//MARK: Close enough
bool DRAGON_Test::approx(double a, double b, double rel, double abs) {
    return fabs(a - b) <= abs + rel * fmax(fabs(a), fabs(b));
}
void DRAGON_Test::expect_close(const ConservativeState& a, const ConservativeState& b, double rel, double abs) {
    assert(approx(a.rho, b.rho, rel, abs));
    assert(approx(a.px,  b.px,  rel, abs));
    assert(approx(a.py,  b.py,  rel, abs));
    assert(approx(a.pz,  b.pz,  rel, abs));
    assert(approx(a.E,   b.E,   rel, abs));
}
void DRAGON_Test::expect_close(const PrimitiveState& a, const PrimitiveState& b, double rel, double abs) {
    assert(approx(a.rho, b.rho, rel, abs));
    assert(approx(a.vx,  b.vx,  rel, abs));
    assert(approx(a.vy,  b.vy,  rel, abs));
    assert(approx(a.vz,  b.vz,  rel, abs));
    assert(approx(a.p,   b.p,   rel, abs));
}


//MARK: Data Structure Verification
void DRAGON_Test::verify_constructors(){
    PrimitiveState W;
    assert(W.rho == 0);
    assert(W.vx == 0);
    assert(W.vy == 0);
    assert(W.vz == 0);
    assert(W.p == 0);

    ConservativeState U;
    assert(U.rho == 0);
    assert(U.px == 0);
    assert(U.py == 0);
    assert(U.pz == 0);
    assert(U.E == 0);
}

void DRAGON_Test::verify_conversion(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 0.0, 10.0);

    ConservativeState U(W);

    assert(approx(U.rho, 2.0));
    assert(approx(U.px, 6.0));
    assert(approx(U.py, 8.0));
    assert(approx(U.pz, 0.0));
    assert(approx(U.E, 40.0));
    assert(approx(U.pressure(), 10.0));
    
    PrimitiveState W2(U);
    expect_close(W, W2);
}

void DRAGON_Test::verify_flux(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 5.0, 10.0);

    ConservativeState U(W);
    ConservativeState F = U.flux(W.vx);

    assert(approx(F.rho, 6.0));
    assert(approx(F.px, 28.0));  // rho vx^2 + p = 2*9 + 10
    assert(approx(F.py, 24.0));  // rho vx vy = 2*3*4
    assert(approx(F.pz, 30.0));  // rho vx vz = 2*3*5
    assert(approx(F.E, W.vx * (U.E + W.p)));
}
void DRAGON_Test::verify_enthalpy(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 0.0, 10.0);
    
    assert(approx(W.enthalpy(), 25.0));
}

//MARK: Dimension Swaps
//Primitive
void DRAGON_Test::verify_swaps_P(){
    PrimitiveState W = make_state(2.0, 1.0, 2.0, 3.0, 10.0);
    
    PrimitiveState XY = W.swapXY();
    assert(approx(XY.vx, 2));
    assert(approx(XY.vy, 1));
    assert(approx(XY.vz, 3));
    
    PrimitiveState XZ = W.swapXZ();
    assert(approx(XZ.vx, 3));
    assert(approx(XZ.vy, 2));
    assert(approx(XZ.vz, 1));
    
    PrimitiveState YZ = W.swapYZ();
    assert(approx(YZ.vx, 1));
    assert(approx(YZ.vy, 3));
    assert(approx(YZ.vz, 2));
    
    //Original Unchanged
    assert(approx(W.vx, 1));
    assert(approx(W.vy, 2));
    assert(approx(W.vz, 3));
}

//Conservative
void DRAGON_Test::verify_swaps_C(){
    ConservativeState U;
    U.px = 1;
    U.py = 2;
    U.pz = 3;

    ConservativeState _XY = U.swapXY();
    assert(approx(_XY.px, 2));
    assert(approx(_XY.py, 1));
    assert(approx(_XY.pz, 3));
    
    ConservativeState _XZ = U.swapXZ();
    assert(approx(_XZ.px, 3));
    assert(approx(_XZ.py, 2));
    assert(approx(_XZ.pz, 1));

    ConservativeState _YZ = U.swapYZ();
    assert(approx(_YZ.px, 1));
    assert(approx(_YZ.py, 3));
    assert(approx(_YZ.pz, 2));
    
    //Original Unchanged
    assert(approx(U.px, 1));
    assert(approx(U.py, 2));
    assert(approx(U.pz, 3));
}



//MARK: Arithemtic

void DRAGON_Test::verify_add(){
    ConservativeState A;
    A.rho = 1; A.px = 2; A.py = 3; A.pz = 4; A.E = 5;
    ConservativeState _A = A;

    ConservativeState B;
    B.rho = 10; B.px = 20; B.py = 30; B.pz = 40; B.E = 50;
    ConservativeState _B = B;
    
    ConservativeState C = A + B;
    assert(approx(C.rho, 11));
    assert(approx(C.px, 22));
    assert(approx(C.py, 33));
    assert(approx(C.pz, 44));
    assert(approx(C.E, 55));

    expect_close(A, _A);   // original unchanged
    expect_close(B, _B);  // original unchanged
    
    A += B;
    assert(approx(A.rho, 11));
    assert(approx(A.px, 22));
    assert(approx(A.py, 33));
    assert(approx(A.pz, 44));
    assert(approx(A.E, 55));
}
void DRAGON_Test::verify_sub(){
    ConservativeState A;
    A.rho = 1; A.px = 2; A.py = 3; A.pz = 4; A.E = 5;
    ConservativeState _A = A;

    ConservativeState B;
    B.rho = 10; B.px = 20; B.py = 30; B.pz = 40; B.E = 50;
    ConservativeState _B = B;
    
    ConservativeState C = B - A;
    assert(approx(C.rho, 9));
    assert(approx(C.px, 18));
    assert(approx(C.py, 27));
    assert(approx(C.pz, 36));
    assert(approx(C.E, 45));

    expect_close(A, _A);   // original unchanged
    expect_close(B, _B);  // original unchanged
    
    B -= A;
    assert(approx(B.rho, 9));
    assert(approx(B.px, 18));
    assert(approx(B.py, 27));
    assert(approx(B.pz, 36));
    assert(approx(B.E, 45));
}
void DRAGON_Test::verify_mult(){
    ConservativeState A;
    A.rho = 1; A.px = 2; A.py = 3; A.pz = 4; A.E = 5;
    ConservativeState _A = A;

    double k = 2;
    
    ConservativeState C = k * A;
    assert(approx(C.rho, 2));
    assert(approx(C.px, 4));
    assert(approx(C.py, 6));
    assert(approx(C.pz, 8));
    assert(approx(C.E, 10));
    
    C = A * k;
    assert(approx(C.rho, 2));
    assert(approx(C.px, 4));
    assert(approx(C.py, 6));
    assert(approx(C.pz, 8));
    assert(approx(C.E, 10));
    
    expect_close(A, _A);   // original unchanged
    
    A *= k;
    assert(approx(A.rho, 2));
    assert(approx(A.px, 4));
    assert(approx(A.py, 6));
    assert(approx(A.pz, 8));
    assert(approx(A.E, 10));
    
    assert(approx(k,2));  // original unchanged
}

void DRAGON_Test::verify_div(){
    ConservativeState A;
    A.rho = 1; A.px = 2; A.py = 3; A.pz = 4; A.E = 5;
    ConservativeState _A = A;

    double k = 2;
    
    ConservativeState C = A / k;
    assert(approx(C.rho, 0.5));
    assert(approx(C.px, 1.0));
    assert(approx(C.py, 1.5));
    assert(approx(C.pz, 2.0));
    assert(approx(C.E, 2.5));
    
    expect_close(A, _A);   // original unchanged
    
    A /= k;
    assert(approx(A.rho, 0.5));
    assert(approx(A.px, 1.0));
    assert(approx(A.py, 1.5));
    assert(approx(A.pz, 2.0));
    assert(approx(A.E, 2.5));

    assert(approx(k,2));  // original unchanged

}
void DRAGON_Test::verify_flux_add() {
    PrimitiveState W= make_state(1.0, 2.0, 3.0, 4.0, 10.0);
    ConservativeState U0(W);

    ConservativeState dU;
    dU.rho = 0.1;
    dU.px  = 0.2;
    dU.py  = 0.3;
    dU.pz  = 0.4;
    dU.E   = 0.5;

    W += dU;

    ConservativeState U1(W);
    ConservativeState expected = U0 + dU;

    assert(approx(U1.rho, expected.rho, 1e-12, 1e-12));
    assert(approx(U1.px,  expected.px,  1e-12, 1e-12));
    assert(approx(U1.py,  expected.py,  1e-12, 1e-12));
    assert(approx(U1.pz,  expected.pz,  1e-12, 1e-12));
    assert(approx(U1.E,   expected.E,   1e-12, 1e-12));
}
