//
//  Boundary/BoundarySet.tpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#pragma once

namespace Boundary {

    template<typename... Bs> BoundarySet::BoundarySet(Bs&&... bs) : BoundaryType(X | Y | Z, true) {
        (addBoundary(std::forward<Bs>(bs)), ...);

        int uncoveredFaces = X | Y | Z;
        for (auto& b : boundaries)  uncoveredFaces &= ~b->get_faces();

        if (uncoveredFaces) {//TODO: Make this at the beginning
            addBoundary(Outflow(uncoveredFaces));
        }
    }

    template<typename B> void BoundarySet::addBoundary(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, BoundaryType>);
        boundaries.push_back(std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
    }

}
