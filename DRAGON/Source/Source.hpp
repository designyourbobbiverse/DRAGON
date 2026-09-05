//
//  Source.hpp
//  DRAGON/Source
//
//  Created by Bobbie Markwick on 04/09/2026.
//

#ifndef Source_hpp
#define Source_hpp


#include "FluidElement/FluidElement.hpp"

namespace DRAGON{

class Source{
public:
    virtual ConservativeState source_density(const PrimitiveState& w, double t) = 0;
    
    //Integrate the source term ODE over dt, starting from initial (w0, t0)
    ConservativeState integrate(double dt, const PrimitiveState& w0, double t0); //Chooses RK2 vs RK4 based on Config.h
    ConservativeState rk2(double dt, const PrimitiveState& w0, double t0);
    ConservativeState rk4(double dt, const PrimitiveState& w0, double t0);
};

//Mass Source is for terms involving changes in mass
//Provide the density source term, and the other components will be computed for you
//Assumes the injected mass is comoving and in thermal equilibrium with the fluid
class MassSource: public Source{
public:
    virtual double thermal_energy(const PrimitiveState& w, double t); //Override this if the source isn't in thermal equilibrium with the fluid
    virtual double density(const PrimitiveState& w, double t) = 0;
    virtual ConservativeState source_density(const PrimitiveState& w, double t) override;
};

//Momentum Source is for net forces
//Provide the force vector, and the energy source f*v will be added automatically
class MomentumSource: public Source{
public:
    virtual vec3 force(const PrimitiveState& w, double t) = 0;
    virtual ConservativeState source_density(const PrimitiveState& w, double t) override;
};
//Source term for adding (or removing) energy via isotropic transfer
class EnergySource: public Source{
public:
    virtual double energy(const PrimitiveState& w, double t) = 0;
    ConservativeState source_density(const PrimitiveState& w, double t) override;
};



}



#endif
