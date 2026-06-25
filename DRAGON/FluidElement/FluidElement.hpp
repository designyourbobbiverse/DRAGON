//
//  FluidElement.hpp
//  DRAGON/FluidElement
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef FluidElement_hpp
#define FluidElement_hpp

#include "Config.h"

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
    //Stored in Gaussian units, be explicit  with 4*M_PI
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
    vec3 p;
#ifdef MHD
    //Magnetic Field
    //Stored in Gaussian units, be explicit  with 4*M_PI
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


//== verifies that all components are equal, to within an absolute difference of 1e-12 computational units
bool operator==(const PrimitiveState &X, const PrimitiveState &Y);
bool operator==(const ConservativeState &X, const ConservativeState &Y);
bool operator==(const vec3 &v, const vec3 &w);


//Arithmetic (+): Add two conservative states or 3-vectors together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y);
vec3 operator+(vec3 v, const vec3 &w);
vec3& operator+=(vec3 &v, const vec3 &w);
//Arithmetic (-): Subtract one state or 3-vector from another of the same type
ConservativeState operator-(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y);
vec3 operator-(vec3 v, const vec3 &w);
vec3& operator-=(vec3 &v, const vec3 &w);
//Arithmetic (*): Multiply a conservative state or 3-vector by some scalar.
ConservativeState operator*(ConservativeState X, double a);//state * scalar
ConservativeState operator*(const double &a, ConservativeState X);//scalar * state
ConservativeState& operator*=(ConservativeState &X, double a); //state *= scalar
vec3 operator*(vec3 v, double a);//vector * scalar
vec3 operator*(const double &a, vec3 v); //scalar * vector
vec3& operator*=(vec3 &v, double a); //vector *= scalar
//Arithmetic (*): Dot product
double operator*(const vec3& v, const vec3& w); //Returns v.x*w.x + v.y*w.y + v.z*w.z

//Arithmetic (/): Divide a conservative state or 3-vector by some scalar
ConservativeState operator/(ConservativeState X, double a);
ConservativeState& operator/=(ConservativeState &X, double a);
vec3 operator/(vec3 v, double a);
vec3& operator/=(vec3 &v, double a);

//Arithmetic (+): Adjust primative state by an amount equal to flux * (1 time unit) / (1 length unit)
PrimitiveState& operator+=(PrimitiveState &X, const ConservativeState &flux);






#endif /* FluidElement_hpp */
