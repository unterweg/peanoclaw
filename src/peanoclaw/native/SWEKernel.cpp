/*
 * PyClaw.cpp
 *
 *  Created on: Feb 18, 2012
 *      Author: kristof
 */

#include "peanoclaw/native/SWEKernel.h"

#include "peanoclaw/Patch.h"
#include "peanoclaw/geometry/Region.h"
#include "peanoclaw/grid/aspects/BoundaryIterator.h"
#include "peanoclaw/grid/boundaryConditions/ReflectingBoundaryCondition.h"
#include "peanoclaw/interSubgridCommunication/DefaultTransfer.h"
#include "peanoclaw/native/SWE_WavePropagationBlock_patch.hh"

#include "tarch/multicore/MulticoreDefinitions.h"
#include "tarch/timing/Watch.h"
#include "tarch/parallel/Node.h"

using namespace peanoclaw::grid::boundaryConditions;

tarch::logging::Log peanoclaw::native::SWEKernel::_log("peanoclaw::native::SWEKernel");

void peanoclaw::native::SWEKernel::transformWaterHeight(
  peanoclaw::Patch& subgrid,
  const peanoclaw::geometry::Region&       region,
  bool              modifyUOld,
  bool              absoluteToAboveSeaFloor
) const {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();
  double sign = absoluteToAboveSeaFloor ? -1 : +1;
  if(modifyUOld) {
    dfor(internalSubcellIndex, region._size) {
      tarch::la::Vector<DIMENSIONS,int> subcellIndex = internalSubcellIndex + region._offset;
      accessor.setValueUOld(
        subcellIndex,
        0,
        accessor.getValueUOld(subcellIndex, 0) + sign * accessor.getParameterWithGhostlayer(subcellIndex, 0)
      );
    }
  } else {
    dfor(internalSubcellIndex, region._size) {
      tarch::la::Vector<DIMENSIONS,int> subcellIndex = internalSubcellIndex + region._offset;
      accessor.setValueUNew(
        subcellIndex,
        0,
        accessor.getValueUNew(subcellIndex, 0) + sign * accessor.getParameterWithGhostlayer(subcellIndex, 0)
      );
    }
  }
}

void peanoclaw::native::SWEKernel::advanceBlockInTime(
  SWE_WavePropagationBlock_patch& block,
  peanoclaw::Patch& subgrid,
  double maximumTimestepSize
) {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();

  block.setArrays(
        subgrid,
        reinterpret_cast<float*>(accessor.getUOldWithGhostLayerArray(0)),
        reinterpret_cast<float*>(accessor.getUOldWithGhostLayerArray(1)),
        reinterpret_cast<float*>(accessor.getUOldWithGhostLayerArray(2)),
        reinterpret_cast<float*>(accessor.getParameterWithGhostLayerArray(0))
      );

  block.computeNumericalFluxes();

  double dt = fmin(block.getMaxTimestep(), maximumTimestepSize);
  double estimatedNextTimestepSize = block.getMaxTimestep();

  block.updateUnknowns(dt);

  peanoclaw::interSubgridCommunication::DefaultTransfer transfer;
  transfer.swapUNewAndUOld(subgrid);

  assertion4(
      tarch::la::greater(subgrid.getTimeIntervals().getTimestepSize(), 0.0)
      || tarch::la::greater(estimatedNextTimestepSize, 0.0)
      || tarch::la::equals(maximumTimestepSize, 0.0)
      || tarch::la::equals(subgrid.getTimeIntervals().getEstimatedNextTimestepSize(), 0.0),
      subgrid, maximumTimestepSize, estimatedNextTimestepSize, subgrid.toStringUNew());
  assertion(subgrid.getTimeIntervals().getTimestepSize() < std::numeric_limits<double>::infinity());

  if (tarch::la::greater(dt, 0.0)) {
    subgrid.getTimeIntervals().advanceInTime();
    subgrid.getTimeIntervals().setTimestepSize(dt);
  }
  subgrid.getTimeIntervals().setEstimatedNextTimestepSize(estimatedNextTimestepSize);
}

peanoclaw::native::SWEKernel::SWEKernel(
  peanoclaw::native::scenarios::SWEScenario& scenario,
  peanoclaw::interSubgridCommunication::DefaultTransfer* transfer,
  peanoclaw::interSubgridCommunication::Interpolation*   interpolation,
  peanoclaw::interSubgridCommunication::Restriction*     restriction,
  peanoclaw::interSubgridCommunication::FluxCorrection*  fluxCorrection
) : Numerics(
    transfer,
    scenario.hasCustomInterpolation() ? scenario.getCustomInterpolation() : interpolation,
    scenario.hasCustomRestriction() ? scenario.getCustomRestriction() : restriction,
    fluxCorrection
),
_totalSolverCallbackTime(0.0),
_scenario(scenario),
_cachedSubdivisionFactor(-1),
_cachedGhostlayerWidth(-1),
_cachedBlock(0)
{

}

peanoclaw::native::SWEKernel::~SWEKernel()
{
}


void peanoclaw::native::SWEKernel::initializePatch(
  Patch& patch
) {
  logTraceIn( "initializePatch(...)" );

  _scenario.initializePatch(patch);

  logTraceOut( "initializePatch(...)" );
}

void peanoclaw::native::SWEKernel::update(Patch& subgrid) {
  _scenario.update(subgrid);
}

void peanoclaw::native::SWEKernel::solveTimestep(
  Patch& subgrid,
  double maximumTimestepSize,
  bool useDimensionalSplitting,
  tarch::la::Vector<DIMENSIONS_TIMES_TWO, bool> domainBoundaryFlags
) {
  logTraceInWith2Arguments( "solveTimestep(...)", maximumTimestepSize, useDimensionalSplitting);

  assertion2(tarch::la::greater(maximumTimestepSize, 0.0), "Timestepsize == 0 should be checked outside.", subgrid.getTimeIntervals().getMinimalNeighborTimeConstraint());

  tarch::timing::Watch solverWatch("", "", false);
  solverWatch.startTimer();

  #ifdef SharedMemoryParallelisation
  SWE_WavePropagationBlock_patch block(subgrid);
  advanceBlockInTime(
    block,
    subgrid,
    maximumTimestepSize
  );
  #else
  if(
    _cachedSubdivisionFactor != subgrid.getSubdivisionFactor()
    || _cachedGhostlayerWidth != subgrid.getGhostlayerWidth()
  ) {
    //SWE_WavePropagationBlock_patch swe(patch);
    _cachedBlock = std::auto_ptr<SWE_WavePropagationBlock_patch>(new SWE_WavePropagationBlock_patch(subgrid));
    _cachedSubdivisionFactor = subgrid.getSubdivisionFactor();
    _cachedGhostlayerWidth = subgrid.getGhostlayerWidth();
  }

  advanceBlockInTime(
    *_cachedBlock,
    subgrid,
    maximumTimestepSize
  );
  #endif

  solverWatch.stopTimer();
  _totalSolverCallbackTime += solverWatch.getCalendarTime();

  logTraceOut( "solveTimestep(...)");
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::native::SWEKernel::getDemandedMeshWidth(Patch& patch, bool isInitializing) {
  return _scenario.computeDemandedMeshWidth(patch, isInitializing);
}

void peanoclaw::native::SWEKernel::addPatchToSolution(Patch& patch) {
}

void peanoclaw::native::SWEKernel::fillBoundaryLayer(Patch& subgrid, int dimension, bool setUpper) {
  logTraceInWith3Arguments("fillBoundaryLayerInPyClaw", patch, dimension, setUpper);

  logDebug("fillBoundaryLayerInPyClaw", "Setting left boundary for " << patch.getPosition() << ", dim=" << dimension << ", setUpper=" << setUpper);

  //std::cout << "------ setUpper " << setUpper << " dimension " << dimension << std::endl;
  //std::cout << patch << std::endl;
  //std::cout << "++++++" << std::endl;
  //std::cout << patch.toStringUOldWithGhostLayer() << std::endl;
  //std::cout << "||||||" << std::endl;

  // implement a wall boundary
  tarch::la::Vector<DIMENSIONS, int> src_subcellIndex;
  tarch::la::Vector<DIMENSIONS, int> dest_subcellIndex;

  peanoclaw::grid::SubgridAccessor& accessor = subgrid.getAccessor();

  //Fill default boundary condition
  ReflectingBoundaryCondition reflectingBoundaryCondition;
  peanoclaw::grid::aspects::BoundaryIterator<ReflectingBoundaryCondition> defaultBoundaryIterator(reflectingBoundaryCondition);
  defaultBoundaryIterator.iterate(subgrid, accessor, dimension, setUpper);

  //Fill scenario boundary condition
  peanoclaw::grid::aspects::BoundaryIterator<peanoclaw::native::scenarios::SWEScenario> scenarioBoundaryIterator(_scenario);
  scenarioBoundaryIterator.iterate(subgrid, accessor, dimension, setUpper);

  logTraceOut("fillBoundaryLayerInPyClaw");
}

void peanoclaw::native::SWEKernel::interpolateSolution (
  const tarch::la::Vector<DIMENSIONS, int>& destinationSize,
  const tarch::la::Vector<DIMENSIONS, int>& destinationOffset,
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool interpolateToUOld,
  bool interpolateToCurrentTime,
  bool useTimeUNewOrTimeUOld
) const {
  if(_scenario.hasCustomInterpolation()) {
    //Interpolate
    Numerics::interpolateSolution (
      destinationSize,
      destinationOffset,
      source,
      destination,
      interpolateToUOld,
      interpolateToCurrentTime,
      useTimeUNewOrTimeUOld
    );
  } else {
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
    transformWaterHeight(source, sourceRegion, true, false); //UOld
    transformWaterHeight(source, sourceRegion, false, false); // UNew

    //Interpolate
    Numerics::interpolateSolution (
      destinationSize,
      destinationOffset,
      source,
      destination,
      interpolateToUOld,
      interpolateToCurrentTime,
      useTimeUNewOrTimeUOld
    );

    //Source: Absolute Water Height -> Water Height above Sea Floor
    transformWaterHeight(source, sourceRegion, true, true); //UOld
    transformWaterHeight(source, sourceRegion, false, true); // UNew

    //Destination: Absolute Water Height -> Water Height above Sea Floor
    transformWaterHeight(destination, destinationRegion, interpolateToUOld, true);
  }
}

void peanoclaw::native::SWEKernel::restrictSolution (
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool              restrictOnlyOverlappedRegions
) const {
  if(_scenario.hasCustomRestriction()) {
    Numerics::restrictSolution(source, destination, restrictOnlyOverlappedRegions);
  } else {
    peanoclaw::geometry::Region sourceRegion(tarch::la::Vector<DIMENSIONS,int>(0), source.getSubdivisionFactor());

    transformWaterHeight(source, sourceRegion, true, false); //UOld
    transformWaterHeight(source, sourceRegion, false, false); //UNew

    Numerics::restrictSolution(
      source,
      destination,
      restrictOnlyOverlappedRegions
    );

    transformWaterHeight(source, sourceRegion, true, true); //UOld
    transformWaterHeight(source, sourceRegion, false, true); //UNew
  }
}

void peanoclaw::native::SWEKernel::postProcessRestriction(
  peanoclaw::Patch& destination,
  bool              restrictOnlyOverlappedRegions
) const {
  if(_scenario.hasCustomRestriction()) {
    Numerics::postProcessRestriction(destination, restrictOnlyOverlappedRegions);
  } else {
    peanoclaw::geometry::Region destinationRegion(tarch::la::Vector<DIMENSIONS,int>(0), destination.getSubdivisionFactor());
    transformWaterHeight(destination, destinationRegion, true, true); //UOld
    transformWaterHeight(destination, destinationRegion, false, true); //UNew
  }
}

