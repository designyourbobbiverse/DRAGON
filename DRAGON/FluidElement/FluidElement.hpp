//
//  FluidElement.hpp
//  DRAGON/FluidElement
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef FluidElement_hpp
#define FluidElement_hpp

#include "Config.h" //#ifdef MHD code only compiles if MHD is defined in Config.h

namespace DRAGON{
struct vec3;
struct PrimitiveState;
struct ConservativeState;

// A 3 dimensional vector, storing the X, Y, and Z components
struct vec3{
    double x, y, z;
    
    vec3(double x=0, double y=0, double z=0): x(x), y(y), z(z){ }
    // Returns a copy with the two specified components swapped
    vec3 swappedXY() const;
    vec3 swappedXZ() const;
    vec3 swappedYZ() const;
    //Swaps the specified components in place.
    void swapXY();
    void swapXZ();
    void swapYZ();

};

// A fluid element in terms of primitive variables: density, velocity, pressure, (magnetic field)
struct PrimitiveState {
    //Density, Pressure
    double rho, p;
    //Velocity
    vec3 v;
    #ifdef MHD
    //Magnetic Field
    //Stored in Gaussian units, be explicit  with 4pi
    vec3 B;
    #endif
    
    PrimitiveState(); //Constructs a zero state
    PrimitiveState(ConservativeState state); //Converts conservative to primitive
    
    
    //Derived Physical Quantities, assuming ideal gas EOS
    double energy() const;//Energy Density
    double enthalpy() const;//Enthalpy per particle
    double cs() const; //Speed of hydrodynamic sound waves
    #ifdef MHD
    double c_alfven() const; //Alfven wave speed
    double c_fast() const; //Wave speed for the fast magnetosonic mode (along x)
    double c_fast_max() const; //Wave speed for the fast magnetosonic mode (fastest of x/y/z)
private:
    double c_fast(double Bk) const; //Wave speed for the fast magnetosonic mode along the k direction
public:
    #endif
    
    //Computes the x-direction flux
    //Flux is stored as a Conservative state equal to [flux] * (1 time unit)/(1 length unit)
    ConservativeState flux() const;
    
    //swapped__() returns a copy with the two specified components swapped
    PrimitiveState swappedXY() const;
    PrimitiveState swappedXZ() const;
    PrimitiveState swappedYZ() const;
    //swap__() swaps the specified components in place.
    void swapXY();
    void swapXZ();
    void swapYZ();
    
    // Checks whether density/pressure are positive, and all vector components are finite
    bool isPhysical() const;
};

// A fluid element in terms of conservated quantities: density, momentum, energy, (magnetic field)
struct ConservativeState {
    //Mass Density, Energy Density
    double rho, E;
    //Momentum density
    vec3 mom;
    #ifdef MHD
    //Magnetic Field
    //Stored in Gaussian units, be explicit  with 4pi
    vec3 B;
    #endif
    
    ConservativeState();//Constructs a zero state
    ConservativeState(const ConservativeState&) = default; //copying
    ConservativeState(PrimitiveState prim);//Converts primitive to conservative
    
    //Computes the pressure, assuming ideal gas EOS
    double pressure() const;
    
    //Computes the x-direction flux
    //Flux is stored as a Conservative state equal to [flux] * (1 time unit)/(1 length unit)
    ConservativeState flux() const;
    ConservativeState flux(vec3 v) const; //More effiecient version if you already know v
    
    //swapped__() returns a copy with the two specified components swapped
    ConservativeState swappedXY() const;
    ConservativeState swappedXZ() const;
    ConservativeState swappedYZ() const;
    //swap__() swaps the specified components in place.
    void swapXY();
    void swapXZ();
    void swapYZ();
    
    // Checks whether density/pressure are positive, and all vector components are finite
    bool isPhysical() const;
    // Checks whether all components are finite
    bool isFinite() const;
};
}

//== verifies that all components are equal, to within an absolute difference of 1e-12 computational units
bool operator==(const DRAGON::PrimitiveState &X, const DRAGON::PrimitiveState &Y);
bool operator==(const DRAGON::ConservativeState &X, const DRAGON::ConservativeState &Y);
bool operator==(const DRAGON::vec3 &v, const DRAGON::vec3 &w);


//Arithmetic (+): Add two conservative states or 3-vectors together
DRAGON::ConservativeState operator+(DRAGON::ConservativeState X, const DRAGON::ConservativeState &Y);
DRAGON::ConservativeState& operator+=(DRAGON::ConservativeState &X, const DRAGON::ConservativeState &Y);
DRAGON::vec3 operator+(DRAGON::vec3 v, const DRAGON::vec3 &w);
DRAGON::vec3& operator+=(DRAGON::vec3 &v, const DRAGON::vec3 &w);
//Arithmetic (-): Subtract one state or 3-vector from another of the same type
DRAGON::ConservativeState operator-(DRAGON::ConservativeState X, const DRAGON::ConservativeState &Y);
DRAGON::ConservativeState& operator-=(DRAGON::ConservativeState &X, const DRAGON::ConservativeState &Y);
DRAGON::vec3 operator-(DRAGON::vec3 v, const DRAGON::vec3 &w);
DRAGON::vec3& operator-=(DRAGON::vec3 &v, const DRAGON::vec3 &w);
//Arithmetic (*): Multiply a conservative state or 3-vector by some scalar.
DRAGON::ConservativeState operator*(DRAGON::ConservativeState X, double a);//state * scalar
DRAGON::ConservativeState operator*(const double &a, DRAGON::ConservativeState X);//scalar * state
DRAGON::ConservativeState& operator*=(DRAGON::ConservativeState &X, double a); //state *= scalar
DRAGON::vec3 operator*(DRAGON::vec3 v, double a);//vector * scalar
DRAGON::vec3 operator*(const double &a, DRAGON::vec3 v); //scalar * vector
DRAGON::vec3& operator*=(DRAGON::vec3 &v, double a); //vector *= scalar
//Arithmetic (*): Dot product
double operator*(const DRAGON::vec3& v, const DRAGON::vec3& w); //Returns v.x*w.x + v.y*w.y + v.z*w.z
//Arithmetic: Cross product
DRAGON::vec3 cross(const DRAGON::vec3& v, const DRAGON::vec3& w);

//Arithmetic (/): Divide a conservative state or 3-vector by some scalar
DRAGON::ConservativeState operator/(DRAGON::ConservativeState X, double a);
DRAGON::ConservativeState& operator/=(DRAGON::ConservativeState &X, double a);
DRAGON::vec3 operator/(DRAGON::vec3 v, double a);
DRAGON::vec3& operator/=(DRAGON::vec3 &v, double a);

//Arithmetic (+): Adjust primative state by an amount equal to flux * (1 time unit) / (1 length unit)
DRAGON::PrimitiveState& operator+=(DRAGON::PrimitiveState &X, const DRAGON::ConservativeState &flux);





#endif /* FluidElement_hpp */
