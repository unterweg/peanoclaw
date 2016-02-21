/*
 * BathymetryRestriction.cpp
 *
 *  Created on: Jan 26, 2015
 *      Author: kristof
 */
#include "peanoclaw/interSubgridCommunication/BathymetryRestriction.h"

#include "peanoclaw/interSubgridCommunication/BathymetryTransformation.h"
#include "peanoclaw/Patch.h"

#include "tarch/Assertions.h"

tarch::logging::Log peanoclaw::interSubgridCommunication::BathymetryRestriction::_log("peanoclaw::interSubgridCommunication::BathymetryRestriction");

void peanoclaw::interSubgridCommunication::BathymetryRestriction::restrictRegion (
  peanoclaw::Patch&                     source,
  const peanoclaw::geometry::Region&    sourceRegion,
  peanoclaw::grid::SubgridIterator<3>&  sourceIterator,
  peanoclaw::Patch&                     destination,
  peanoclaw::grid::SubgridIterator<3>&  destinationIterator
) {
  logTraceInWith2Arguments("restrictRegion", source.toString(), destination.toString());
  assertionEquals(source.getUnknownsPerSubcell(), destination.getUnknownsPerSubcell());
  assertion(destination.isVirtual());

  double epsilon = 1e-12;
  const int NUMBER_OF_UNKNOWNS = 3;

  //TODO unterweg restricting to interval [0, 1]
  #ifdef PEANOCLAW_USE_ASCEND_FOR_RESTRICTION
  double destinationTimeUOld = destination.getTimeIntervals().getMinimalFineGridTime();
  double destinationTimeUNew = destinationTimeUOld + destination.getTimeIntervals().getMaximalFineGridTimestep();
  #else
  double destinationTimeUOld = 0.0;// destination.getTimeUOld();
  double destinationTimeUNew = 1.0;// destination.getTimeUNew();
  #endif

  const tarch::la::Vector<DIMENSIONS, double> sourceSubcellSize = tarch::la::multiplyComponents(source.getSize(), tarch::la::invertEntries(source.getSubdivisionFactor().convertScalar<double>()));
  const tarch::la::Vector<DIMENSIONS, double> destinationSubcellSize = tarch::la::multiplyComponents(destination.getSize(), tarch::la::invertEntries(destination.getSubdivisionFactor().convertScalar<double>()));
  const double destinationSubcellRegion = tarch::la::volume(destinationSubcellSize);
  const tarch::la::Vector<DIMENSIONS, double> sourcePosition = source.getPosition();

  //Time factor
  double timeFactorUOld = 1.0;
  double timeFactorUNew = 1.0;
  if(tarch::la::greater(source.getTimeIntervals().getTimestepSize(), 0.0)) {
    timeFactorUOld = (destinationTimeUOld - source.getTimeIntervals().getTimeUOld()) / (source.getTimeIntervals().getTimeUNew() - source.getTimeIntervals().getTimeUOld());
    timeFactorUNew = (destinationTimeUNew - source.getTimeIntervals().getTimeUOld()) / (source.getTimeIntervals().getTimeUNew() - source.getTimeIntervals().getTimeUOld());
  }

  assertion4(timeFactorUOld == timeFactorUOld, destinationTimeUOld, destinationTimeUNew, source.getTimeIntervals().getCurrentTime(), source.getTimeIntervals().getTimestepSize());
  assertion4(timeFactorUNew == timeFactorUNew, destinationTimeUOld, destinationTimeUNew, source.getTimeIntervals().getCurrentTime(), source.getTimeIntervals().getTimestepSize());

  //Destination region
  peanoclaw::geometry::Region destinationRegion = sourceRegion.mapToPatch(source, destination, epsilon);
  destinationIterator.restart(destinationRegion._offset, destinationRegion._size);

  //Loop through region in destination
  while(destinationIterator.moveToNextCell()) {
    tarch::la::Vector<DIMENSIONS, int> destinationSubcellIndex = destinationIterator.getCellIndex();
    tarch::la::Vector<DIMENSIONS, double> destinationSubcellPosition = destinationIterator.getCellPosition();

    assertion4(tarch::la::allGreaterEquals(destinationSubcellIndex, tarch::la::Vector<DIMENSIONS, int>(0))
              && tarch::la::allGreater(destination.getSubdivisionFactor(), destinationSubcellIndex), sourceRegion, destinationRegion, destination, source);

    //Get region for single destination cell
    peanoclaw::geometry::Region subcellRegion = sourceRegion.mapCellToPatch(sourcePosition, sourceSubcellSize, destinationSubcellSize, destinationSubcellIndex, destinationSubcellPosition, epsilon);
    sourceIterator.restart(subcellRegion._offset, subcellRegion._size);
    assertion7(tarch::la::allGreaterEquals(subcellRegion._size, 0), subcellRegion, sourcePosition, sourceSubcellSize, destinationSubcellSize, destinationSubcellPosition, source, destination);

    //Loop through region in source
    while(sourceIterator.moveToNextCell()) {
      tarch::la::Vector<DIMENSIONS, double> sourceSubcellPosition = sourceIterator.getCellPosition();

      assertion6(tarch::la::allGreaterEquals(sourceIterator.getCellIndex(), tarch::la::Vector<DIMENSIONS, int>(0))
                && tarch::la::allGreater(source.getSubdivisionFactor(), sourceIterator.getCellIndex()),
                destinationSubcellIndex, subcellRegion, sourceRegion, destinationRegion, destination, source);

      double overlapRegion = calculateOverlappingRegion(
        destinationSubcellPosition,
        destinationSubcellSize,
        sourceSubcellPosition,
        sourceSubcellSize
      );

      tarch::la::Vector<NUMBER_OF_UNKNOWNS, double> sourceValuesUOld = sourceIterator.getUnknownsUOld();
      tarch::la::Vector<NUMBER_OF_UNKNOWNS, double> sourceValuesUNew = sourceIterator.getUnknownsUNew();

      tarch::la::Vector<NUMBER_OF_UNKNOWNS, double> deltaUOld = (sourceValuesUOld * (1.0 - timeFactorUOld) + sourceValuesUNew * timeFactorUOld)
                                                                * overlapRegion / destinationSubcellRegion;
      tarch::la::Vector<NUMBER_OF_UNKNOWNS, double> deltaUNew = (sourceValuesUOld * (1.0 - timeFactorUNew) + sourceValuesUNew * timeFactorUNew)
                                                                * overlapRegion / destinationSubcellRegion;

      destinationIterator.setUnknownsUOld(destinationIterator.getUnknownsUOld() + deltaUOld);
      destinationIterator.setUnknownsUNew(destinationIterator.getUnknownsUNew() + deltaUNew);
    }
  }

  assertion8(!destination.containsNaN(), timeFactorUNew, timeFactorUOld, destination, destination.toStringUNew(), destination.toStringUOldWithGhostLayer(), source, source.toStringUNew(), source.toStringUOldWithGhostLayer());

  logTraceOut("restrictRegion");
}

peanoclaw::interSubgridCommunication::BathymetryRestriction::BathymetryRestriction(
  double seaLevel
) : _seaLevel(seaLevel) {
}


int peanoclaw::interSubgridCommunication::BathymetryRestriction::restrictSolution (
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool              restrictOnlyOverlappedRegions
) {
  const int NUMBER_OF_UNKNOWNS = 3;
  assertionEquals(source.getUnknownsPerSubcell(), NUMBER_OF_UNKNOWNS);
  peanoclaw::grid::SubgridAccessor sourceAccessor = source.getAccessor();
  peanoclaw::grid::SubgridIterator<NUMBER_OF_UNKNOWNS> sourceIterator
    = sourceAccessor.getSubgridIterator<NUMBER_OF_UNKNOWNS>(tarch::la::Vector<DIMENSIONS, int>(0), source.getSubdivisionFactor());
  peanoclaw::grid::SubgridAccessor destinationAccessor = destination.getAccessor();
  peanoclaw::grid::SubgridIterator<NUMBER_OF_UNKNOWNS> destinationIterator
    = destinationAccessor.getSubgridIterator<NUMBER_OF_UNKNOWNS>(tarch::la::Vector<DIMENSIONS, int>(0), destination.getSubdivisionFactor());

  BathymetryTransformation transformation(_seaLevel);
  peanoclaw::geometry::Region sourceRegion(tarch::la::Vector<DIMENSIONS,int>(0), source.getSubdivisionFactor());
  transformation.transformWaterHeight(source, sourceRegion, true, false); //UOld
  transformation.transformWaterHeight(source, sourceRegion, false, false); //UNew

  int numberOfRestrictedCells = 0;
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
        restrictRegion(source, regions[i], sourceIterator, destination, destinationIterator);
        numberOfRestrictedCells += tarch::la::volume(regions[i]._size);
      }
    }
  } else {
    //Restrict complete Patch
    peanoclaw::geometry::Region region;
    region._offset = tarch::la::Vector<DIMENSIONS, int>(0);
    region._size = source.getSubdivisionFactor();
    restrictRegion(source, region, sourceIterator, destination, destinationIterator);
    numberOfRestrictedCells += tarch::la::volume(region._size);
  }

  transformation.transformWaterHeight(source, sourceRegion, true, true); //UOld
  transformation.transformWaterHeight(source, sourceRegion, false, true); //UNew

  return numberOfRestrictedCells;
}

void peanoclaw::interSubgridCommunication::BathymetryRestriction::postProcessRestriction(
  peanoclaw::Patch& destination,
  bool              restrictOnlyOverlappedRegions
) const {
  peanoclaw::geometry::Region destinationRegion(tarch::la::Vector<DIMENSIONS,int>(0), destination.getSubdivisionFactor());
  BathymetryTransformation transformation(_seaLevel);
  transformation.transformWaterHeight(destination, destinationRegion, true, true); //UOld
  transformation.transformWaterHeight(destination, destinationRegion, false, true); //UNew
}



