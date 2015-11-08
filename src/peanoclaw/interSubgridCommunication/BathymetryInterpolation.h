/*
 * BathymetryInterpolation.h
 *
 *  Created on: Oct 4, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYINTERPOLATION_H_
#define PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYINTERPOLATION_H_

#include "peanoclaw/geometry/Region.h"
#include "peanoclaw/interSubgridCommunication/Interpolation.h"

#include "peano/utils/Dimensions.h"

namespace peanoclaw {
  namespace interSubgridCommunication {
    class BathymetryInterpolation;

    template<int NumberOfUnknowns>
    class BathymetryInterpolationTemplate;
  }
}

/**
 * Implements a d-linear interpolation of grid values.
 */
class peanoclaw::interSubgridCommunication::BathymetryInterpolation : public peanoclaw::interSubgridCommunication::Interpolation {
  private:
    double _seaLevel;

    void performDefaultInterpolation(
      const tarch::la::Vector<DIMENSIONS, int>&    destinationSize,
      const tarch::la::Vector<DIMENSIONS, int>&    destinationOffset,
      peanoclaw::Patch& source,
      peanoclaw::Patch& destination,
      bool interpolateToUOld,
      bool interpolateToCurrentTime,
      bool useTimeUNewOrTimeUOld
    );

  public:
    BathymetryInterpolation(double seaLevel) : _seaLevel(seaLevel) {}

    void interpolateSolution (
      const tarch::la::Vector<DIMENSIONS, int>&    destinationSize,
      const tarch::la::Vector<DIMENSIONS, int>&    destinationOffset,
      peanoclaw::Patch& source,
      peanoclaw::Patch& destination,
      bool interpolateToUOld,
      bool interpolateToCurrentTime,
      bool useTimeUNewOrTimeUOld
    );
};


#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYINTERPOLATION_H_ */
