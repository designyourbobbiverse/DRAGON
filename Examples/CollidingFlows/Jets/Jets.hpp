//
//  Jets.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 03/07/2026.
//
#include "Boundary.hpp"

//Jet Boundary condition
//Constructs an inflow of (rho,v,p) on the specified face
//Only overrides the boundary on the specified face, not outside the jet
class Jet : public GhostFill {
public:
    Jet(double rho, double v, double p, double radius, std::string face);
    Jet(double rho, double v, double p, double radius, int face);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
protected:
    double rho, v, p, rj;
    int jetface; //Distinct from Ghostfill.face to allow implicit outflow outside of jet radius
    
};

#ifdef MHD
//MHD Jet Boundary condition
//Similar to Jet, but includes a pressure-balanced toroidal field
class MHDJet : public Jet {
public:
    MHDJet(double rho, double v, double p_amb, double beta, double rm, double rj, std::string face);
    MHDJet(double rho, double v, double p_amb, double beta, double rm, double rj, int face);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
protected:
    double beta, rm;
};
#endif
