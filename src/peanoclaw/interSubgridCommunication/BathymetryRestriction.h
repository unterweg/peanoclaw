/*
 * BathymetryRestriction.h
 *
 *  Created on: Jan 26, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYRESTRICTION_H_
#define PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYRESTRICTION_H_

#include "peanoclaw/interSubgridCommunication/Restriction.h"
#include "peanoclaw/grid/SubgridAccessor.h"

namespace peanoclaw {
  namespace interSubgridCommunication {
    class BathymetryRestriction;
  }
}

class peanoclaw::interSubgridCommunication::BathymetryRestriction : public Restriction {

  private:
    void restrictRegion (
        peanoclaw::Patch&                  source,
        const peanoclaw::geometry::Region& sourceRegion,
        peanoclaw::grid::SubgridAccessor&  sourceIterator,
        peanoclaw::Patch&                  destination,
        peanoclaw::grid::SubgridAccessor&  destinationIterator
    );

  public:
    virtual ~BathymetryRestriction(){};

    /**
     * Restricts data from a fine patch to a coarse patch.
     */
    virtual void restrictSolution (
      peanoclaw::Patch& source,
      peanoclaw::Patch& destination,
      bool              restrictOnlyOverlappedRegions
    );
};


#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYRESTRICTION_H_ */
