/*
 * NonOvershootingInterpolation.h
 *
 *  Created on: Oct 12, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_INTERSUBGRIDCOMMUNICATION_NONOVERSHOOTINGINTERPOLATION_H_
#define PEANOCLAW_INTERSUBGRIDCOMMUNICATION_NONOVERSHOOTINGINTERPOLATION_H_

#include "peanoclaw/interSubgridCommunication/Interpolation.h"
#include "peanoclaw/Patch.h"

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

namespace peanoclaw {
  namespace interSubgridCommunication {
    class NonOvershootingInterpolation;

    template<int NumberOfUnknowns>
    class NonOvershootingInterpolationTemplate;
  }
}

/**
 * Default implementation for interpolation grid values from coarse
 * to fine subgrids.
 */
template<int NumberOfUnknowns>
class peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate {
  private:
    /**
     * Logging device for the trace macros.
     */
    static tarch::logging::Log  _log;

    int _signLookupTable[TWO_POWER_D_TIMES_D];

    void interpolateUnknowns(
      peanoclaw::grid::SubgridAccessor& sourceAccessor,
      peanoclaw::grid::SubgridAccessor& destinationAccessor,
      const tarch::la::Vector<DIMENSIONS,int>& neighborIndexInSourcePatch,
      const tarch::la::Vector<DIMENSIONS,int>& subcellIndexInDestinationPatch,
      double temporalFactor,
      double spatialFactor,
      bool interpolateToUOld
    ) const;

    /**
     * Sign is used to account for extrapolation if destination cell center is outside the source cell.
     */
    template<int Dimensionality>
    double getSpatialFactor(
      const tarch::la::Vector<Dimensionality, double>& destinationSubcellCenter,
      const tarch::la::Vector<Dimensionality, double>& neighborPositionInSourcePatch,
      const tarch::la::Vector<Dimensionality, double>& sourceSubcellSize,
      const tarch::la::Vector<Dimensionality, double>& destinationSubcellSize,
      const tarch::la::Vector<Dimensionality, int>& offset,
      const bool extrapolating
    ) const;

    template<int Dimensionality>
    tarch::la::Vector<Dimensionality, int> getSubcellInSourceSubgrid(
      const tarch::la::Vector<Dimensionality, double>& subcellPositionInDestinationSubgrid,
      const tarch::la::Vector<Dimensionality, double>& sourcePosition,
      const tarch::la::Vector<Dimensionality, double>& sourceSubcellSize,
      const tarch::la::Vector<Dimensionality, int>& sourceSubdivisionFactor,
      const tarch::la::Vector<Dimensionality, double>& inverseSourceSubcellSize,
      bool& extrapolating
    ) const;

  public:
    NonOvershootingInterpolationTemplate();

    /**
     * @see peanoclaw::interSubgridCommunication::Interpolation
     *
     * TODO unterweg: useTimeUNewOrTimeUOld is only used because the
     * restriction now is done to the minimalNeighborTimeInterval. However,
     * during grid refinement the interpolation has to take the actual currentTime
     * and timestepSize of the subgrid into account and not timeUNew and timeUOld
     * from the TimeIntervals class. This could be circumvented if the restriction
     * is carried out in the ascend event and can be done to the minimalFineGrid
     * interval which becomes the new currentTime and timestepSize of the virtual
     * subgrid anyway.
     */
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

/**
 * Implements a d-linear interpolation of grid values.
 */
class peanoclaw::interSubgridCommunication::NonOvershootingInterpolation : public peanoclaw::interSubgridCommunication::Interpolation {
  public:
    void interpolateSolution (
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
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<1> transfer1;
            transfer1.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 2:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<2> transfer2;
            transfer2.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 3:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<3> transfer3;
            transfer3.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 4:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<4> transfer4;
            transfer4.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 5:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<5> transfer5;
            transfer5.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 6:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<6> transfer6;
            transfer6.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 7:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<7> transfer7;
            transfer7.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 8:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<8> transfer8;
            transfer8.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 9:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<9> transfer9;
            transfer9.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        case 10:
          {
            peanoclaw::interSubgridCommunication::NonOvershootingInterpolationTemplate<10> transfer10;
            transfer10.interpolateSolution(destinationSize, destinationOffset, source, destination, interpolateToUOld, interpolateToCurrentTime, useTimeUNewOrTimeUOld);
          }
          break;
        default:
          assertionFail("Number of unknowns " << source.getUnknownsPerSubcell() << " not supported!");
      }
    }
};

#include "peanoclaw/interSubgridCommunication/NonOvershootingInterpolation.cpph"

#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_NONOVERSHOOTINGINTERPOLATION_H_ */
