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
    }

    template<typename B> void Boundaries::addBoundary(B&& b) {
        static_assert(std::derived_from<std::decay_t<B>, GhostFill>);
        boundaries.push_back(std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
        stale = true;
    }

}
