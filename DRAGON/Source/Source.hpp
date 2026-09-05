//
//  Source.hpp
//  DRAGON/Source
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 04/09/2026.
//

#ifndef Source_hpp
#define Source_hpp


#include "FluidElement/FluidElement.hpp"

namespace DRAGON::Source{

//To define a source term, subclass Source.
class SourceTerm{
public:
    //This function should return the source terms (in conservative form)
    virtual ConservativeState source_density(const PrimitiveState& w, double t) = 0;
    
    //Integrate the source term ODE over dt, starting from initial (w0, t0)
    ConservativeState integrate(double dt, const PrimitiveState& w0, double t0=0); //Chooses RK2 vs RK4 based on Config.h
    ConservativeState rk2(double dt, const PrimitiveState& w0, double t0=0);
    ConservativeState rk4(double dt, const PrimitiveState& w0, double t0=0);
};


//MARK: Source Types
//These are intermediate-parent classes.
//Subclass these to add a source term to one conservative component without needing to manually compute the derived components (e.g. provide force f and the work f*v will be computed for you).

//Mass Sources
//Provide the density source term, and the other components will be computed for you
//Assumes the injected mass is comoving and in thermal equilibrium with the fluid
class MassSource: public SourceTerm{
public:
    virtual double density(const PrimitiveState& w, double t) = 0; //Provide the source term for the rho component
    virtual vec3 velocity(const PrimitiveState& w, double t); //Override this if the source isn't comoving with the fluid
    virtual double thermal_energy(const PrimitiveState& w, double t); //Override this if the source isn't in thermal equilibrium with the fluid
    virtual ConservativeState source_density(const PrimitiveState& w, double t) override;
};

//Momentum Sources for net forces
//Provide the force vector, energy term f*v will be added automatically
class MomentumSource: public SourceTerm{
public:
    virtual vec3 force(const PrimitiveState& w, double t) = 0; //Provide the source term for the rho*v component
    virtual ConservativeState source_density(const PrimitiveState& w, double t) override;
};
//Energy Sources for adding (or removing) energy via isotropic transfer
class EnergySource: public SourceTerm{
public:
    virtual double energy(const PrimitiveState& w, double t) = 0; //Provide the source term fro the E component
    ConservativeState source_density(const PrimitiveState& w, double t) override;
};


}



#endif
