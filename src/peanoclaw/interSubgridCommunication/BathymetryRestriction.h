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

/**
 * This restriction type is used for the Storegga scenario.
 */
class peanoclaw::interSubgridCommunication::BathymetryRestriction : public Restriction {

  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    double _seaLevel;

    void restrictRegion (
        peanoclaw::Patch&                     source,
        const peanoclaw::geometry::Region&    sourceRegion,
        peanoclaw::grid::SubgridIterator<3>&  sourceIterator,
        peanoclaw::Patch&                     destination,
        peanoclaw::grid::SubgridIterator<3>&  destinationIterator
    );

  public:
    BathymetryRestriction(double seaLevel);

    virtual ~BathymetryRestriction(){};

    /**
     * Restricts data from a fine patch to a coarse patch.
     */
    int restrictSolution (
      peanoclaw::Patch& source,
      peanoclaw::Patch& destination,
      bool              restrictOnlyOverlappedRegions
    );

    void postProcessRestriction(
      peanoclaw::Patch& destination,
      bool              restrictOnlyOverlappedRegions
    ) const;
};


#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYRESTRICTION_H_ */
