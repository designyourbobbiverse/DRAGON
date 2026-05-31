//
//  Boundary/BoundarySet.tpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#pragma once

namespace Boundary {

    template<typename... Bs> Boundaries::Boundaries(Bs&&... bs) : GhostFill(X | Y | Z, true) {
        (addBoundary(std::forward<Bs>(bs)), ...);

        int uncoveredFaces = X | Y | Z;
        for (auto& b : boundaries)  uncoveredFaces &= ~b->get_faces();

        if (uncoveredFaces) {
            boundaries.insert(boundaries.begin(), std::make_unique<Outflow>(uncoveredFaces) );
        }
    }

    template<typename B> void Boundaries::addBoundary(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, GhostFill>);
        boundaries.push_back(std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
    }

}
