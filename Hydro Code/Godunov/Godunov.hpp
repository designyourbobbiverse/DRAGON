//
//  Godunov.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Godunov_hpp
#define Godunov_hpp

#include <stdio.h>
#include "FluidElement.hpp"

struct FluidLine{
    int size;

    FluidLine(int size, int ghosts=0);
    
    void advance(double dt, double dx);
    void advance(double dt, double dx, FluidLine& _L, FluidLine& _R);
    
    PrimitiveState& operator[](int k);
    const PrimitiveState& operator[](int k) const;

    ~FluidLine();
private:
    PrimitiveState* w;
    int ghosts;
    
};

#endif /* Godunov_hpp */
