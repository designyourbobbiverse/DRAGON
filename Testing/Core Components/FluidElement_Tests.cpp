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
    if(output) std::cout << "Fluid Arithmetic: \n";
    if(output) std::cout << "- Element Construction: ";
    verify_constructors();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Form Conversion: ";
    verify_conversion();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- enthalpy: ";
    verify_enthalpy();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Axis Swaps: ";
    verify_swaps_P();
    verify_swaps_C();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Addition: ";
    verify_add();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Subtraction: ";
    verify_sub();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Multiplication: ";
    verify_mult();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Division: ";
    verify_div();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Flux Calculation: ";
    verify_flux();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Flux Addition: ";
    verify_flux_add();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "All Fluid Arithmetic Tests Passed.\n\n";
}


PrimitiveState DRAGON_Test::make_state(double rho, double vx, double vy, double vz, double p) {
    PrimitiveState W;
    W.rho = rho;
    W.v.x = vx;
    W.v.y = vy;
    W.v.z = vz;
    W.p = p;
    return W;
}


//MARK: Close enough
bool DRAGON_Test::approx(double a, double b, double rel, double abs) {
    return fabs(a - b) <= abs + rel * fmax(fabs(a), fabs(b));
}
void DRAGON_Test::expect_close(const ConservativeState& a, const ConservativeState& b, double rel, double abs) {
    assert(approx(a.rho, b.rho, rel, abs));
    assert(approx(a.p.x,  b.p.x,  rel, abs));
    assert(approx(a.p.y,  b.p.y,  rel, abs));
    assert(approx(a.p.z,  b.p.z,  rel, abs));
    assert(approx(a.E,   b.E,   rel, abs));
}
void DRAGON_Test::expect_close(const PrimitiveState& a, const PrimitiveState& b, double rel, double abs) {
    assert(approx(a.rho, b.rho, rel, abs));
    assert(approx(a.v.x,  b.v.x,  rel, abs));
    assert(approx(a.v.y,  b.v.y,  rel, abs));
    assert(approx(a.v.z,  b.v.z,  rel, abs));
    assert(approx(a.p,   b.p,   rel, abs));
}


//MARK: Data Structure Verification
void DRAGON_Test::verify_constructors(){
    PrimitiveState W;
    assert(W.rho == 0);
    assert(W.v.x == 0);
    assert(W.v.y == 0);
    assert(W.v.z == 0);
    assert(W.p == 0);

    ConservativeState U;
    assert(U.rho == 0);
    assert(U.p.x == 0);
    assert(U.p.y == 0);
    assert(U.p.z == 0);
    assert(U.E == 0);
}

void DRAGON_Test::verify_conversion(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 0.0, 10.0);

    ConservativeState U(W);

    assert(approx(U.rho, 2.0));
    assert(approx(U.p.x, 6.0));
    assert(approx(U.p.y, 8.0));
    assert(approx(U.p.z, 0.0));
    assert(approx(U.E, 40.0));
    assert(approx(U.pressure(), 10.0));
    
    PrimitiveState W2(U);
    expect_close(W, W2);
}

void DRAGON_Test::verify_flux(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 5.0, 10.0);

    ConservativeState U(W);
    ConservativeState F = U.flux(W.v);

    assert(approx(F.rho, 6.0));
    assert(approx(F.p.x, 28.0));  // rho v.x^2 + p = 2*9 + 10
    assert(approx(F.p.y, 24.0));  // rho v.x v.y = 2*3*4
    assert(approx(F.p.z, 30.0));  // rho v.x v.z = 2*3*5
    assert(approx(F.E, W.v.x * (U.E + W.p)));
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
    assert(approx(XY.v.x, 2));
    assert(approx(XY.v.y, 1));
    assert(approx(XY.v.z, 3));
    
    PrimitiveState XZ = W.swapXZ();
    assert(approx(XZ.v.x, 3));
    assert(approx(XZ.v.y, 2));
    assert(approx(XZ.v.z, 1));
    
    PrimitiveState YZ = W.swapYZ();
    assert(approx(YZ.v.x, 1));
    assert(approx(YZ.v.y, 3));
    assert(approx(YZ.v.z, 2));
    
    //Original Unchanged
    assert(approx(W.v.x, 1));
    assert(approx(W.v.y, 2));
    assert(approx(W.v.z, 3));
}

//Conservative
void DRAGON_Test::verify_swaps_C(){
    ConservativeState U;
    U.p.x = 1;
    U.p.y = 2;
    U.p.z = 3;

    ConservativeState _XY = U.swapXY();
    assert(approx(_XY.p.x, 2));
    assert(approx(_XY.p.y, 1));
    assert(approx(_XY.p.z, 3));
    
    ConservativeState _XZ = U.swapXZ();
    assert(approx(_XZ.p.x, 3));
    assert(approx(_XZ.p.y, 2));
    assert(approx(_XZ.p.z, 1));

    ConservativeState _YZ = U.swapYZ();
    assert(approx(_YZ.p.x, 1));
    assert(approx(_YZ.p.y, 3));
    assert(approx(_YZ.p.z, 2));
    
    //Original Unchanged
    assert(approx(U.p.x, 1));
    assert(approx(U.p.y, 2));
    assert(approx(U.p.z, 3));
}



//MARK: Arithemtic

void DRAGON_Test::verify_add(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
    ConservativeState _A = A;

    ConservativeState B;
    B.rho = 10; B.p.x = 20; B.p.y = 30; B.p.z = 40; B.E = 50;
    ConservativeState _B = B;
    
    ConservativeState C = A + B;
    assert(approx(C.rho, 11));
    assert(approx(C.p.x, 22));
    assert(approx(C.p.y, 33));
    assert(approx(C.p.z, 44));
    assert(approx(C.E, 55));

    expect_close(A, _A);   // original unchanged
    expect_close(B, _B);  // original unchanged
    
    A += B;
    assert(approx(A.rho, 11));
    assert(approx(A.p.x, 22));
    assert(approx(A.p.y, 33));
    assert(approx(A.p.z, 44));
    assert(approx(A.E, 55));
}
void DRAGON_Test::verify_sub(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
    ConservativeState _A = A;

    ConservativeState B;
    B.rho = 10; B.p.x = 20; B.p.y = 30; B.p.z = 40; B.E = 50;
    ConservativeState _B = B;
    
    ConservativeState C = B - A;
    assert(approx(C.rho, 9));
    assert(approx(C.p.x, 18));
    assert(approx(C.p.y, 27));
    assert(approx(C.p.z, 36));
    assert(approx(C.E, 45));

    expect_close(A, _A);   // original unchanged
    expect_close(B, _B);  // original unchanged
    
    B -= A;
    assert(approx(B.rho, 9));
    assert(approx(B.p.x, 18));
    assert(approx(B.p.y, 27));
    assert(approx(B.p.z, 36));
    assert(approx(B.E, 45));
}
void DRAGON_Test::verify_mult(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
    ConservativeState _A = A;

    double k = 2;
    
    ConservativeState C = k * A;
    assert(approx(C.rho, 2));
    assert(approx(C.p.x, 4));
    assert(approx(C.p.y, 6));
    assert(approx(C.p.z, 8));
    assert(approx(C.E, 10));
    
    C = A * k;
    assert(approx(C.rho, 2));
    assert(approx(C.p.x, 4));
    assert(approx(C.p.y, 6));
    assert(approx(C.p.z, 8));
    assert(approx(C.E, 10));
    
    expect_close(A, _A);   // original unchanged
    
    A *= k;
    assert(approx(A.rho, 2));
    assert(approx(A.p.x, 4));
    assert(approx(A.p.y, 6));
    assert(approx(A.p.z, 8));
    assert(approx(A.E, 10));
    
    assert(approx(k,2));  // original unchanged
}

void DRAGON_Test::verify_div(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
    ConservativeState _A = A;

    double k = 2;
    
    ConservativeState C = A / k;
    assert(approx(C.rho, 0.5));
    assert(approx(C.p.x, 1.0));
    assert(approx(C.p.y, 1.5));
    assert(approx(C.p.z, 2.0));
    assert(approx(C.E, 2.5));
    
    expect_close(A, _A);   // original unchanged
    
    A /= k;
    assert(approx(A.rho, 0.5));
    assert(approx(A.p.x, 1.0));
    assert(approx(A.p.y, 1.5));
    assert(approx(A.p.z, 2.0));
    assert(approx(A.E, 2.5));

    assert(approx(k,2));  // original unchanged

}
void DRAGON_Test::verify_flux_add() {
    PrimitiveState W= make_state(1.0, 2.0, 3.0, 4.0, 10.0);
    ConservativeState U0(W);

    ConservativeState dU;
    dU.rho = 0.1;
    dU.p.x  = 0.2;
    dU.p.y  = 0.3;
    dU.p.z  = 0.4;
    dU.E   = 0.5;

    W += dU;

    ConservativeState U1(W);
    ConservativeState expected = U0 + dU;

    assert(approx(U1.rho, expected.rho, 1e-12, 1e-12));
    assert(approx(U1.p.x,  expected.p.x,  1e-12, 1e-12));
    assert(approx(U1.p.y,  expected.p.y,  1e-12, 1e-12));
    assert(approx(U1.p.z,  expected.p.z,  1e-12, 1e-12));
    assert(approx(U1.E,   expected.E,   1e-12, 1e-12));
}
