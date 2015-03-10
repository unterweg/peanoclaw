/*
 * BathymetryRestriction.cpp
 *
 *  Created on: Jan 26, 2015
 *      Author: kristof
 */
#include "peanoclaw/interSubgridCommunication/BathymetryRestriction.h"

#include "peanoclaw/Patch.h"

#include "tarch/Assertions.h"

void peanoclaw::interSubgridCommunication::BathymetryRestriction::restrictRegion (
  peanoclaw::Patch&                  source,
  const peanoclaw::geometry::Region& sourceRegion,
  peanoclaw::grid::SubgridAccessor&  sourceIterator,
  peanoclaw::Patch&                  destination,
  peanoclaw::grid::SubgridAccessor&  destinationIterator
) {
//  dfor(localSourceSubcellIndex, sourceRegion._size) {
//    tarch::la::Vector<DIMENSIONS, int> sourceSubcellIndex = localSourceSubcellIndex + sourceRegion._offset;
//  }
}


void peanoclaw::interSubgridCommunication::BathymetryRestriction::restrictSolution (
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool              restrictOnlyOverlappedRegions
) {
  peanoclaw::grid::SubgridAccessor sourceAccessor = source.getAccessor();
//  peanoclaw::grid::SubgridIterator<NumberOfUnknowns> sourceIterator
//    = sourceAccessor.getSubgridIterator<NumberOfUnknowns>(tarch::la::Vector<DIMENSIONS, int>(0), source.getSubdivisionFactor());
  peanoclaw::grid::SubgridAccessor destinationAccessor = destination.getAccessor();
//  peanoclaw::grid::SubgridIterator<NumberOfUnknowns> destinationIterator
//    = destinationAccessor.getSubgridIterator<NumberOfUnknowns>(tarch::la::Vector<DIMENSIONS, int>(0), destination.getSubdivisionFactor());

  if(restrictOnlyOverlappedRegions) {
    peanoclaw::geometry::Region regions[DIMENSIONS_TIMES_TWO];
    int numberOfRegionsToProcess = peanoclaw::geometry::Region::getRegionsOverlappedByNeighboringGhostlayers(
        destination.getLowerNeighboringGhostlayerBounds(),
        destination.getUpperNeighboringGhostlayerBounds(),
        source.getPosition(),
        source.getSize(),
        source.getSubcellSize(),
        source.getSubdivisionFactor(),
        regions
    );

    logDebug("restrictSolution(...)", "Restriction from patch " << source << std::endl << " to patch " << destination);
    for( int i = 0; i < numberOfRegionsToProcess; i++ ) {

      assertion1(tarch::la::allGreaterEquals(regions[i]._size, 0), regions[i]);

      logDebug("restrictSolution(...)", "Restricting region [" << regions[i]._offset << "], [" << regions[i]._size << "]");
      if(tarch::la::allGreater(regions[i]._size, tarch::la::Vector<DIMENSIONS, int>(0))) {
        restrictRegion(source, regions[i], sourceAccessor, destination, destinationAccessor);
      }
    }
  } else {
    //Restrict complete Patch
    peanoclaw::geometry::Region region;
    region._offset = tarch::la::Vector<DIMENSIONS, int>(0);
    region._size = source.getSubdivisionFactor();
    restrictRegion(source, region, sourceAccessor, destination, destinationAccessor);
  }
}



