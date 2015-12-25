/*
 * DefaultRestriction.h
 *
 *  Created on: Mar 7, 2013
 *      Author: kristof
 */

#ifndef PEANOCLAW_INTERSUBGRIDCOMMUNICATION_DEFAULTRESTRICTION_H_
#define PEANOCLAW_INTERSUBGRIDCOMMUNICATION_DEFAULTRESTRICTION_H_

#include "peanoclaw/grid/SubgridAccessor.h"
#include "peanoclaw/interSubgridCommunication/Restriction.h"
#include "peanoclaw/tests/GhostLayerCompositorTest.h"
#include "peanoclaw/tests/GridLevelTransferTest.h"

#include "tarch/la/Vector.h"

#include <algorithm>

namespace peanoclaw {
  namespace interSubgridCommunication {
    class DefaultRestriction;

    template<int NumberOfUnknowns>
    class DefaultRestrictionTemplate;
  }
}

template<int NumberOfUnknowns>
class peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate {

  private:
    /**
     * Logging device for the trace macros.
     */
    static tarch::logging::Log  _log;

    friend class peanoclaw::tests::GhostLayerCompositorTest;
    friend class peanoclaw::tests::GridLevelTransferTest;

    /**
     * Restricts the data from the given source to the destination averaging
     * over the source cells that overlap with a destination cell.
     *
     * @param source The patch to get the data from
     * @param destination The patch to restrict the data to
     * @param destinationTime The time that should be used as the destination
     * time when interpolating in time
     * @param restrictToUOld Decides wether to restrict to UOld or to UNew
     */
    void restrictRegion (
      peanoclaw::Patch&                                   source,
      peanoclaw::Patch&                                   destination,
      peanoclaw::grid::SubgridIterator<NumberOfUnknowns>& sourceIterator,
      peanoclaw::grid::SubgridIterator<NumberOfUnknowns>& destinationIterator,
      const peanoclaw::geometry::Region&                  region
    );

  public:
    int restrictSolution (
      peanoclaw::Patch& source,
      peanoclaw::Patch&       destination,
      bool restrictOnlyOverlappedRegions
    );
};

/**
 * Performs a normal averaging restriction between grid levels.
 */
class peanoclaw::interSubgridCommunication::DefaultRestriction
    : public peanoclaw::interSubgridCommunication::Restriction {
  public:
    int restrictSolution (
      peanoclaw::Patch& source,
      peanoclaw::Patch&       destination,
      bool restrictOnlyOverlappedRegions
    ) {
      switch(source.getUnknownsPerSubcell()) {
        case 1:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<1> transfer1;
          return transfer1.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 2:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<2> transfer2;
          return transfer2.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 3:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<3> transfer3;
          return transfer3.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 4:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<4> transfer4;
          return transfer4.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 5:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<5> transfer5;
          return transfer5.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 6:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<6> transfer6;
          return transfer6.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 7:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<7> transfer7;
          return transfer7.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 8:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<8> transfer8;
          return transfer8.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 9:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<9> transfer9;
          return transfer9.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        case 10:
          peanoclaw::interSubgridCommunication::DefaultRestrictionTemplate<10> transfer10;
          return transfer10.restrictSolution(source, destination, restrictOnlyOverlappedRegions);
        default:
          assertionFail("Number of unknowns " << source.getUnknownsPerSubcell() << " not supported!");
      }
    }

    virtual void postProcessRestriction(
      peanoclaw::Patch& destination,
      bool              restrictOnlyOverlappedRegions
    ) const {}
};

#include "peanoclaw/interSubgridCommunication/DefaultRestriction.cpph"

#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_DEFAULTRESTRICTION_H_ */
