/*
 * Restriction.h
 *
 *  Created on: Mar 5, 2013
 *      Author: kristof
 */

#ifndef PEANOCLAW_INTERSUBGRIDCOMMUNICATION_RESTRICTION_H_
#define PEANOCLAW_INTERSUBGRIDCOMMUNICATION_RESTRICTION_H_

#include "peano/utils/Globals.h"
#include "tarch/logging/Log.h"

#include "tarch/la/Vector.h"

namespace peanoclaw {

  class Patch;

  namespace geometry {
    class Region;
  }

  namespace interSubgridCommunication {
    class Restriction;
  }
}

namespace peanoclaw {
  namespace interSubgridCommunication {
    /**
     * Returns the region of the region where the two given patches overlap.
     */
    inline double calculateOverlappingRegion(
        const tarch::la::Vector<DIMENSIONS, double>& position1,
        const tarch::la::Vector<DIMENSIONS, double>& size1,
        const tarch::la::Vector<DIMENSIONS, double>& position2,
        const tarch::la::Vector<DIMENSIONS, double>& size2
    ) {
      double region = 1.0;

      for(int d = 0; d < DIMENSIONS; d++) {
        double overlappingInterval =
            std::min(position1(d)+size1(d), position2(d)+size2(d))
        - std::max(position1(d), position2(d));
        region *= overlappingInterval;

        region = std::max(region, 0.0);
      }

      return region;
    }
  }
 }

class peanoclaw::interSubgridCommunication::Restriction {

  public:

    virtual ~Restriction(){}

    /**
     * Restricts data from a fine patch to a coarse patch.
     */
    virtual void restrictSolution (
      peanoclaw::Patch& source,
      peanoclaw::Patch& destination,
      bool              restrictOnlyOverlappedRegions
    ) = 0;
};

#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_RESTRICTION_H_ */
