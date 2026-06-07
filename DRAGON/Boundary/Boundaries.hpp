//
//  Boundary/BoundarySet.tpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#pragma once

namespace Boundary {

//MARK: Setup

    template<BoundaryElement... Bs> Boundaries::Boundaries(Bs&&... bs) : GhostFill(X | Y | Z, true) {
        (append(std::forward<Bs>(bs)), ...);
    }

    inline void Boundaries::resetImplicit(){
        int uncoveredFaces = X | Y | Z;
        for (auto& b : boundaries)  uncoveredFaces &= ~b->get_faces();
        implicits = Outflow(uncoveredFaces);
        stale = false;
    }

//MARK: Apply
    inline void Boundary::Boundaries::apply(Grid1D &grid){
        if(stale) resetImplicit();
        implicits.apply(grid);
        for (auto& b : boundaries)  b->apply(grid);
    }
    inline void Boundary::Boundaries::apply(Grid2D &grid){
        if(stale) resetImplicit();
        implicits.apply(grid);
        for (auto& b : boundaries)  b->apply(grid);
    }
    inline void Boundary::Boundaries::apply(Grid3D &grid){
        if(stale) resetImplicit();
        implicits.apply(grid);
        for (auto& b : boundaries)  b->apply(grid);
    }

//MARK: Append
    template<BoundaryElement B> void Boundaries::append(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, GhostFill>);
        boundaries.push_back(std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
        stale = true;
    }
    inline void Boundaries::append(Boundaries&& bs) {
        for (auto& b : bs.boundaries)  boundaries.push_back(std::move(b));
        bs.boundaries.clear();
        stale = true;
    }
    
//MARK: Prepend
   template<BoundaryElement B> void Boundaries::prepend(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, GhostFill>);
        boundaries.insert(boundaries.begin(), std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
        stale = true;
    }
    inline void Boundaries::prepend(Boundaries&& bs) {
        bs.append(std::move(*this));
        boundaries = std::move(bs.boundaries);
        stale = true;
    }
    
//MARK: + Operators
    template<BoundaryElement A, BoundaryElement B> Boundaries operator+(A&& a, B&& b) {
        Boundaries result;
        result.append(std::forward<A>(a));
        result.append(std::forward<B>(b));
        return result;
    }
    template<BoundaryElement B> Boundaries operator+(Boundaries lhs, B&& rhs) { lhs += rhs; return lhs; }
    template<BoundaryElement B> Boundaries& operator+=(Boundaries& lhs, B&& rhs){
        lhs.append(std::forward<B>(rhs));
        return lhs;
    }
    template<BoundaryElement A> Boundaries operator+(A&& lhs, Boundaries rhs){
        rhs.prepend(std::forward<A>(lhs));
        return rhs;
    }
    inline Boundaries operator+(Boundaries lhs, Boundaries rhs) { lhs += std::move(rhs); return lhs; }
    inline Boundaries& operator+=(Boundaries& lhs, Boundaries rhs) {
        lhs.append(std::move(rhs));
        return lhs;
    }
    

    
    
    
}
