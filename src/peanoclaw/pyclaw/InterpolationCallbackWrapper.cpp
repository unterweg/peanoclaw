/*
 * InterpolationCallback.cpp
 *
 *  Created on: Mar 7, 2013
 *      Author: kristof
 */
#include "peanoclaw/pyclaw/InterpolationCallbackWrapper.h"
#include "peanoclaw/pyclaw/PyClawState.h"
#include "peanoclaw/Patch.h"

peanoclaw::pyclaw::InterpolationCallbackWrapper::InterpolationCallbackWrapper(
  InterPatchCommunicationCallback interpolationCallback
) : _interpolationCallback(interpolationCallback) {

}

peanoclaw::pyclaw::InterpolationCallbackWrapper::~InterpolationCallbackWrapper() {

}

void peanoclaw::pyclaw::InterpolationCallbackWrapper::interpolateSolution (
  const tarch::la::Vector<DIMENSIONS, int>&    destinationSize,
  const tarch::la::Vector<DIMENSIONS, int>&    destinationOffset,
  peanoclaw::Patch& source,
  peanoclaw::Patch& destination,
  bool interpolateToUOld,
  bool interpolateToCurrentTime,
  bool useTimeUNewOrTimeUOld
) {

  PyClawState sourceState(source);
  PyClawState destinationState(destination);

  _interpolationCallback(
    sourceState._q,
    sourceState._qbc,
    sourceState._aux,
    source.getSubdivisionFactor()(0),
    source.getSubdivisionFactor()(1),
    #ifdef Dim3
    source.getSubdivisionFactor()(2),
    #else
    0,
    #endif
    source.getSize()(0),
    source.getSize()(1),
    #ifdef Dim3
    source.getSize()(2),
    #else
    0.0,
    #endif
    source.getPosition()(0),
    source.getPosition()(1),
    #ifdef Dim3
    source.getPosition()(2),
    #else
    0.0,
    #endif
    source.getTimeIntervals().getCurrentTime(),
    source.getTimeIntervals().getTimestepSize(),

    destinationState._q,
    destinationState._qbc,
    destinationState._aux,
    destination.getSubdivisionFactor()(0),
    destination.getSubdivisionFactor()(1),
    #ifdef Dim3
    destination.getSubdivisionFactor()(2),
    #else
    0,
    #endif
    destination.getSize()(0),
    destination.getSize()(1),
    #ifdef Dim3
    destination.getSize()(2),
    #else
    0.0,
    #endif
    destination.getPosition()(0),
    destination.getPosition()(1),
    #ifdef Dim3
    destination.getPosition()(2),
    #else
    0.0,
    #endif
    destination.getTimeIntervals().getCurrentTime(),
    destination.getTimeIntervals().getTimestepSize(),
    source.getUnknownsPerSubcell(),
    source.getNumberOfParametersWithoutGhostlayerPerSubcell()
  );
}



