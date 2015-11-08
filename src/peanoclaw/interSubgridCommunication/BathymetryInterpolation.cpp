/*
 * BathymetryInterpolation.cpp
 *
 *  Created on: Oct 4, 2015
 *      Author: kristof
 */
#include "peanoclaw/interSubgridCommunication/BathymetryInterpolation.h"

#include "peanoclaw/Patch.h"
#include "peanoclaw/interSubgridCommunication/BathymetryTransformation.h"
#include "peanoclaw/interSubgridCommunication/DefaultInterpolation.h"
#include "peanoclaw/grid/SubgridAccessor.h"
#include "peano/utils/Loop.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif

void peanoclaw::interSubgridCommunication::BathymetryInterpolation::performDefaultInterpolation(
  const tarch::la::Vector<DIMENSIONS, int>&    destinationSize,
  const tarch::la::Vector<DIMENSIONS, int>&    destinationOffset,
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool interpolateToUOld,
  bool interpolateToCurrentTime,
  bool useTimeUNewOrTimeUOld
) {
  switch(source.getUnknownsPerSubcell()) {
    case 1:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<1> transfer1;
        transfer1.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 2:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<2> transfer2;
        transfer2.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 3:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<3> transfer3;
        transfer3.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 4:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<4> transfer4;
        transfer4.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 5:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<5> transfer5;
        transfer5.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 6:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<6> transfer6;
        transfer6.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 7:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<7> transfer7;
        transfer7.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 8:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<8> transfer8;
        transfer8.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 9:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<9> transfer9;
        transfer9.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    case 10:
      {
        peanoclaw::interSubgridCommunication::DefaultInterpolationTemplate<10> transfer10;
        transfer10.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
      }
      break;
    default:
      assertionFail("Number of unknowns " << source.getUnknownsPerSubcell() << " not supported!");
  }
}

void peanoclaw::interSubgridCommunication::BathymetryInterpolation::interpolateSolution (
  const tarch::la::Vector<DIMENSIONS, int>&    destinationSize,
  const tarch::la::Vector<DIMENSIONS, int>&    destinationOffset,
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool interpolateToUOld,
  bool interpolateToCurrentTime,
  bool useTimeUNewOrTimeUOld
) {
  peanoclaw::grid::SubgridAccessor sourceAccessor = source.getAccessor();
  peanoclaw::grid::SubgridAccessor destinationAccessor = destination.getAccessor();
  tarch::la::Vector<DIMENSIONS,int> sourceSubdivisionFactor = source.getSubdivisionFactor();

  peanoclaw::geometry::Region destinationRegion(destinationOffset, destinationSize);
  peanoclaw::geometry::Region sourceRegion = destinationRegion.mapToPatch(destination, source);

  //Increase sourceRegion by one cell in each direction.
  for(int d = 0; d < DIMENSIONS; d++) {
    if(sourceRegion._offset[d] > 0) {
      sourceRegion._offset[d] = sourceRegion._offset[d]-1;
      sourceRegion._size[d] = std::min(sourceSubdivisionFactor[d], sourceRegion._offset[d] + sourceRegion._size[d] + 2) - sourceRegion._offset[d];
    } else {
      sourceRegion._size[d] = std::min(sourceSubdivisionFactor[d], sourceRegion._offset[d] + sourceRegion._size[d] + 1) - sourceRegion._offset[d];;
    }
  }

  //Source: Water Height above Sea Floor -> Absolute Water Height
  BathymetryTransformation transformation(_seaLevel);
  transformation.transformWaterHeight(source, sourceRegion, true, false); //UOld
  transformation.transformWaterHeight(source, sourceRegion, false, false); // UNew

  //Interpolate
  performDefaultInterpolation(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);

  //Source: Absolute Water Height -> Water Height above Sea Floor
  transformation.transformWaterHeight(source, sourceRegion, true, true); //UOld
  transformation.transformWaterHeight(source, sourceRegion, false, true); // UNew

  //Destination: Absolute Water Height -> Water Height above Sea Floor
  transformation.transformWaterHeight(destination, destinationRegion, interpolateToUOld, true);
}
