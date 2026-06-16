//
//  FluidElement_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Testing.hpp"
#include "Constants.h"
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
    if(output) std::cout << "- Enthalpy: ";
    verify_enthalpy();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Wave Speeds: ";
    verify_wavespeeds();
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
#ifdef MHD
    assert(approx(a.B.x,  b.B.x,  rel, abs));
    assert(approx(a.B.y,  b.B.y,  rel, abs));
    assert(approx(a.B.z,  b.B.z,  rel, abs));
#endif
}
void DRAGON_Test::expect_close(const PrimitiveState& a, const PrimitiveState& b, double rel, double abs) {
    assert(approx(a.rho, b.rho, rel, abs));
    assert(approx(a.v.x,  b.v.x,  rel, abs));
    assert(approx(a.v.y,  b.v.y,  rel, abs));
    assert(approx(a.v.z,  b.v.z,  rel, abs));
    assert(approx(a.p,   b.p,   rel, abs));
#ifdef MHD
    assert(approx(a.B.x,  b.B.x,  rel, abs));
    assert(approx(a.B.y,  b.B.y,  rel, abs));
    assert(approx(a.B.z,  b.B.z,  rel, abs));
#endif
}


//MARK: Data Structure Verification
void DRAGON_Test::verify_constructors(){
    PrimitiveState W;
    assert(W.rho == 0);
    assert(W.v.x == 0);
    assert(W.v.y == 0);
    assert(W.v.z == 0);
    assert(W.p == 0);
#ifdef MHD
    assert(W.B.x == 0);
    assert(W.B.y == 0);
    assert(W.B.z == 0);
#endif

    ConservativeState U;
    assert(U.rho == 0);
    assert(U.p.x == 0);
    assert(U.p.y == 0);
    assert(U.p.z == 0);
    assert(U.E == 0);
#ifdef MHD
    assert(U.B.x == 0);
    assert(U.B.y == 0);
    assert(U.B.z == 0);
#endif
}

//MARK: Physics Verification
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
    
#ifdef MHD
    W.B = {1.0,2.0,3.0};

    ConservativeState U2(W);

    assert(approx(U2.rho, 2.0));
    assert(approx(U2.p.x, 6.0));
    assert(approx(U2.p.y, 8.0));
    assert(approx(U2.p.z, 0.0));
    assert(approx(U2.E, 40.0 + (14.0)/(8*M_PI) ));
    assert(approx(U2.pressure(), 10.0));
    
    PrimitiveState W3(U2);
    expect_close(W, W3);
#endif
}

void DRAGON_Test::verify_flux(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 5.0, 10.0);

    ConservativeState U(W);
    ConservativeState F = U.flux();

    assert(approx(F.rho, 6.0));
    assert(approx(F.p.x, 28.0));  // rho v.x^2 + p = 2*9 + 10
    assert(approx(F.p.y, 24.0));  // rho v.x v.y = 2*3*4
    assert(approx(F.p.z, 30.0));  // rho v.x v.z = 2*3*5
    assert(approx(F.E, W.v.x * (U.E + W.p)));
    
#ifdef MHD
    W.B = {1.0,2.0,3.0};
    F = W.flux();

    assert(approx(F.rho, 6.0));
    assert(approx(F.p.x, 28.0 + 12/(8*M_PI) ));  // rho v.x^2 + p = 2*9 + 10
    assert(approx(F.p.y, 24.0 - 2/(4*M_PI)));  // rho v.x v.y = 2*3*4
    assert(approx(F.p.z, 30.0- 3/(4*M_PI)));  // rho v.x v.z = 2*3*5
    assert(approx(F.B.x, 0));
    assert(approx(F.B.y, W.v.x * W.B.y - W.v.y * W.B.x));
    assert(approx(F.B.z, W.v.x * W.B.z - W.v.z * W.B.x));

    assert(approx(F.E, W.v.x * (W.energy() + W.p + (W.B*W.B)/(8*M_PI)) - W.v * (W.B.x * W.B)/(4*M_PI) ));

#endif
}
void DRAGON_Test::verify_enthalpy(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 0.0, 10.0);
    
    assert(approx(W.enthalpy(), 25.0));
#ifdef MHD
    W.B = {1.0,2.0,3.0};
    assert(approx(W.enthalpy(), 25.0 + (14.0)/(8*M_PI) ));
#endif
}

void DRAGON_Test::verify_wavespeeds(){
    PrimitiveState W = make_state(2.0, 3.0, 4.0, 0.0, 10.0);
    
    assert(approx(W.cs(), sqrt(_gamma * 5.0)));
#ifdef MHD
    W.B = {1,2,3};
    assert(approx(W.c_alfven(), sqrt(14/(8*M_PI)) ));
    double c2 = _gamma * 5.0, a2 = 14/(8*M_PI), ax2 = 1/(8*M_PI) ;
    assert(approx(W.c_fast(), sqrt((c2+a2)/2 + sqrt(pow(c2+a2,2)/4 - ax2*c2))  ));

#endif
    
}

//MARK: Dimension Swaps
//Primitive
void DRAGON_Test::verify_swaps_P(){
    PrimitiveState W = make_state(2.0, 1.0, 2.0, 3.0, 10.0);
#ifdef MHD
    W.B = {1.0,2.0,3.0};
#endif
    
    PrimitiveState XY = W.swapXY();
    assert(approx(XY.v.x, 2));
    assert(approx(XY.v.y, 1));
    assert(approx(XY.v.z, 3));
#ifdef MHD
    assert(approx(XY.B.x, 2));
    assert(approx(XY.B.y, 1));
    assert(approx(XY.B.z, 3));
#endif

    PrimitiveState XZ = W.swapXZ();
    assert(approx(XZ.v.x, 3));
    assert(approx(XZ.v.y, 2));
    assert(approx(XZ.v.z, 1));
#ifdef MHD
    assert(approx(XZ.B.x, 3));
    assert(approx(XZ.B.y, 2));
    assert(approx(XZ.B.z, 1));
#endif
    
    PrimitiveState YZ = W.swapYZ();
    assert(approx(YZ.v.x, 1));
    assert(approx(YZ.v.y, 3));
    assert(approx(YZ.v.z, 2));
#ifdef MHD
    assert(approx(YZ.B.x, 1));
    assert(approx(YZ.B.y, 3));
    assert(approx(YZ.B.z, 2));
#endif
    
    //Original Unchanged
    assert(approx(W.v.x, 1));
    assert(approx(W.v.y, 2));
    assert(approx(W.v.z, 3));
#ifdef MHD
    assert(approx(W.B.x, 1));
    assert(approx(W.B.y, 2));
    assert(approx(W.B.z, 3));
#endif
}

//Conservative
void DRAGON_Test::verify_swaps_C(){
    ConservativeState U;
    U.p = {1,2,3};
#ifdef MHD
    U.B = {1,2,3};
#endif
    
    ConservativeState _XY = U.swapXY();
    assert(approx(_XY.p.x, 2));
    assert(approx(_XY.p.y, 1));
    assert(approx(_XY.p.z, 3));
#ifdef MHD
    assert(approx(_XY.B.x, 2));
    assert(approx(_XY.B.y, 1));
    assert(approx(_XY.B.z, 3));
#endif
    
    ConservativeState _XZ = U.swapXZ();
    assert(approx(_XZ.p.x, 3));
    assert(approx(_XZ.p.y, 2));
    assert(approx(_XZ.p.z, 1));
#ifdef MHD
    assert(approx(_XZ.B.x, 3));
    assert(approx(_XZ.B.y, 2));
    assert(approx(_XZ.B.z, 1));
#endif

    ConservativeState _YZ = U.swapYZ();
    assert(approx(_YZ.p.x, 1));
    assert(approx(_YZ.p.y, 3));
    assert(approx(_YZ.p.z, 2));
#ifdef MHD
    assert(approx(_YZ.B.x, 1));
    assert(approx(_YZ.B.y, 3));
    assert(approx(_YZ.B.z, 2));
#endif
    
    //Original Unchanged
    assert(approx(U.p.x, 1));
    assert(approx(U.p.y, 2));
    assert(approx(U.p.z, 3));
#ifdef MHD
    assert(approx(U.B.x, 1));
    assert(approx(U.B.y, 2));
    assert(approx(U.B.z, 3));
#endif
}



//MARK: Arithemtic

void DRAGON_Test::verify_add(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
#ifdef MHD
    A.B = {1,2,3};
#endif
    ConservativeState _A = A;

    ConservativeState B;
    B.rho = 10; B.p.x = 20; B.p.y = 30; B.p.z = 40; B.E = 50;
#ifdef MHD
    B.B = {0.1,0.2,0.3};
#endif
    ConservativeState _B = B;

    ConservativeState C = A + B;
    assert(approx(C.rho, 11));
    assert(approx(C.p.x, 22));
    assert(approx(C.p.y, 33));
    assert(approx(C.p.z, 44));
    assert(approx(C.E, 55));
    assert(approx(C.B.x, 1.1));
    assert(approx(C.B.y, 2.2));
    assert(approx(C.B.z, 3.3));

    expect_close(A, _A);   // original unchanged
    expect_close(B, _B);  // original unchanged
    
    A += B;
    assert(approx(A.rho, 11));
    assert(approx(A.p.x, 22));
    assert(approx(A.p.y, 33));
    assert(approx(A.p.z, 44));
    assert(approx(A.E, 55));
    assert(approx(A.B.x, 1.1));
    assert(approx(A.B.y, 2.2));
    assert(approx(A.B.z, 3.3));
    
}
void DRAGON_Test::verify_sub(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
#ifdef MHD
    A.B = {1,2,3};
#endif
    ConservativeState _A = A;

    ConservativeState B;
    B.rho = 10; B.p.x = 20; B.p.y = 30; B.p.z = 40; B.E = 50;
#ifdef MHD
    B.B = {0.1,0.2,0.3};
#endif
    ConservativeState _B = B;
    
    ConservativeState C = B - A;
    assert(approx(C.rho, 9));
    assert(approx(C.p.x, 18));
    assert(approx(C.p.y, 27));
    assert(approx(C.p.z, 36));
    assert(approx(C.E, 45));
    assert(approx(C.B.x, -0.9));
    assert(approx(C.B.y, -1.8));
    assert(approx(C.B.z, -2.7));

    expect_close(A, _A);   // original unchanged
    expect_close(B, _B);  // original unchanged
    
    B -= A;
    assert(approx(B.rho, 9));
    assert(approx(B.p.x, 18));
    assert(approx(B.p.y, 27));
    assert(approx(B.p.z, 36));
    assert(approx(B.E, 45));
    assert(approx(B.B.x, -0.9));
    assert(approx(B.B.y, -1.8));
    assert(approx(B.B.z, -2.7));
}
void DRAGON_Test::verify_mult(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
#ifdef MHD
    A.B = {1,2,3};
#endif
    ConservativeState _A = A;

    double k = 2;
    
    ConservativeState C = k * A;
    assert(approx(C.rho, 2));
    assert(approx(C.p.x, 4));
    assert(approx(C.p.y, 6));
    assert(approx(C.p.z, 8));
    assert(approx(C.E, 10));
    assert(approx(C.B.x, 2));
    assert(approx(C.B.y, 4));
    assert(approx(C.B.z, 6));
    
    C = A * k;
    assert(approx(C.rho, 2));
    assert(approx(C.p.x, 4));
    assert(approx(C.p.y, 6));
    assert(approx(C.p.z, 8));
    assert(approx(C.E, 10));
    assert(approx(C.B.x, 2));
    assert(approx(C.B.y, 4));
    assert(approx(C.B.z, 6));
    
    expect_close(A, _A);   // original unchanged
    
    A *= k;
    assert(approx(A.rho, 2));
    assert(approx(A.p.x, 4));
    assert(approx(A.p.y, 6));
    assert(approx(A.p.z, 8));
    assert(approx(A.E, 10));
    assert(approx(A.B.x, 2));
    assert(approx(A.B.y, 4));
    assert(approx(A.B.z, 6));
    assert(approx(k,2));  // original unchanged
}

void DRAGON_Test::verify_div(){
    ConservativeState A;
    A.rho = 1; A.p.x = 2; A.p.y = 3; A.p.z = 4; A.E = 5;
#ifdef MHD
    A.B = {1,2,3};
#endif
    ConservativeState _A = A;

    double k = 2;
    
    ConservativeState C = A / k;
    assert(approx(C.rho, 0.5));
    assert(approx(C.p.x, 1.0));
    assert(approx(C.p.y, 1.5));
    assert(approx(C.p.z, 2.0));
    assert(approx(C.E, 2.5));
    assert(approx(C.B.x, 0.5));
    assert(approx(C.B.y, 1));
    assert(approx(C.B.z, 1.5));
    
    expect_close(A, _A);   // original unchanged
    
    A /= k;
    assert(approx(A.rho, 0.5));
    assert(approx(A.p.x, 1.0));
    assert(approx(A.p.y, 1.5));
    assert(approx(A.p.z, 2.0));
    assert(approx(A.E, 2.5));
    
    assert(approx(A.B.x, 0.5));
    assert(approx(A.B.y, 1));
    assert(approx(A.B.z, 1.5));
    assert(approx(k,2));  // original unchanged

}
void DRAGON_Test::verify_flux_add() {
    PrimitiveState W= make_state(1.0, 2.0, 3.0, 4.0, 10.0);
#ifdef MHD
    W.B = {0.1,0.2,0.3};
#endif
    ConservativeState U0(W);

    ConservativeState dU;
    dU.rho = 0.1;
    dU.p.x  = 0.2;
    dU.p.y  = 0.3;
    dU.p.z  = 0.4;
    dU.E   = 0.5;
#ifdef MHD
    dU.B.x  = 0.6;
    dU.B.y  = 0.7;
    dU.B.z  = 0.8;
#endif

    W += dU;

    ConservativeState U1(W);
    ConservativeState expected = U0 + dU;

    assert(approx(U1.rho, expected.rho, 1e-12, 1e-12));
    assert(approx(U1.p.x,  expected.p.x,  1e-12, 1e-12));
    assert(approx(U1.p.y,  expected.p.y,  1e-12, 1e-12));
    assert(approx(U1.p.z,  expected.p.z,  1e-12, 1e-12));
    assert(approx(U1.E,   expected.E,   1e-12, 1e-12));
}
