//
//  Boundary/BoundaryList.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#pragma once

#include <vector>
#include <memory>
#include <concepts>
#include <type_traits>
#include <utility>

namespace Boundary {


class BoundaryList : public GhostFill {
public:
    //Can construct something like BoundaryList(Reflective(X), Periodic(Y))
    template<BoundaryElement... Bs> BoundaryList(Bs&&... bs) : GhostFill(X | Y | Z, true) {
        (append(std::forward<Bs>(bs)), ...);
    }
    
    BoundaryList(const BoundaryList&) = delete; //No copying
    BoundaryList(BoundaryList&&) noexcept = default;
    BoundaryList& operator=(BoundaryList&&) noexcept = default;
    template<BoundaryElement B> BoundaryList& operator=(B&& b) {
        clear();
        append(std::forward<B>(b));
        return *this;
    }

    
//MARK: Append
    //Add an element to the end of the list
    template<BoundaryElement B> void append(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, GhostFill>);
        boundaries.push_back(std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
        stale = true;//Require recalculation of implied outflow
    }
    //Add (explicit) elements from another list to the end of the list
    inline void append(BoundaryList&& bs) {
        for (auto& b : bs.boundaries)  boundaries.push_back(std::move(b));
        bs.boundaries.clear();
        stale = true;//Require recalculation of implied outflow
    }
        
//MARK: Prepend
    //Add an element to the beginning of the list (but after any implicit outflow)
   template<BoundaryElement B> void prepend(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, GhostFill>);
        boundaries.insert(boundaries.begin(), std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
        stale = true;//Require recalculation of implied outflow
    }
    //Add (explicit) elements from another list to the beginning of this list (but after any implicit outflow)
    inline void prepend(BoundaryList&& bs) {
        bs.append(std::move(*this));
        boundaries = std::move(bs.boundaries);
        stale = true;//Require recalculation of implied outflow
    }
//MARK: Clear
    //Remove all elements from the list
    void clear() {
        boundaries.clear();
        stale = true;//Require recalculation of implied outflow
    }

//MARK: Apply
    
    void apply(Grid1D &grid) override{
        if(stale) resetImplicit();//Recalculate implicit outflow faces
        implicits.apply(grid);//Apply implicit outflow
        for (auto& b : boundaries)  b->apply(grid);//Apply the explicit boundaries in order
    }
    void apply(Grid2D &grid) override{
        if(stale) resetImplicit();//Recalculate implicit outflow faces
        implicits.apply(grid);//Apply implicit outflow
        for (auto& b : boundaries)  b->apply(grid);//Apply the explicit boundaries in order
    }
    inline void apply(Grid3D &grid) override{
        if(stale) resetImplicit();//Recalculate implicit outflow faces
        implicits.apply(grid);//Apply implicit outflow
        for (auto& b : boundaries)  b->apply(grid);//Apply the explicit boundaries in order
    }
    
private:
    std::vector<std::unique_ptr<GhostFill>>  boundaries;
//MARK: Implied Outflow
    Outflow implicits;
    bool stale = true;//If true, needs to recalculate implicit outflow before applying boundary
    void resetImplicit(){
        int uncoveredFaces = X | Y | Z; //Determine which faces are not handled by any elements in the list
        for (auto& b : boundaries)  uncoveredFaces &= ~b->get_faces();
        implicits = Outflow(uncoveredFaces);
        stale = false;//Reset complete
    }
};


//MARK: + Operators
    template<BoundaryElement A, BoundaryElement B> BoundaryList operator+(A&& a, B&& b) {
        BoundaryList result;
        result.append(std::forward<A>(a));
        result.append(std::forward<B>(b));
        return result;
    }
    template<BoundaryElement B> BoundaryList operator+(BoundaryList lhs, B&& rhs) { lhs += rhs; return lhs; }
    template<BoundaryElement B> BoundaryList& operator+=(BoundaryList& lhs, B&& rhs){
        lhs.append(std::forward<B>(rhs));
        return lhs;
    }
    template<BoundaryElement A> BoundaryList operator+(A&& lhs, BoundaryList rhs){
        rhs.prepend(std::forward<A>(lhs));
        return rhs;
    }
    inline BoundaryList operator+(BoundaryList lhs, BoundaryList rhs) { lhs += std::move(rhs); return lhs; }
    inline BoundaryList& operator+=(BoundaryList& lhs, BoundaryList rhs) {
        lhs.append(std::move(rhs));
        return lhs;
    }
    
    
}
