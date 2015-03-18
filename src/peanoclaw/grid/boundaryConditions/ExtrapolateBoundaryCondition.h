/*
 * ExtrapolateBoundaryCondition.h
 *
 *  Created on: Feb 24, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_GRID_BOUNDARYCONDITIONS_EXTRAPOLATEBOUNDARYCONDITION_H_
#define PEANOCLAW_GRID_BOUNDARYCONDITIONS_EXTRAPOLATEBOUNDARYCONDITION_H_

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

namespace peanoclaw {
  class Patch;
  namespace grid {
    class SubgridAccessor;
    namespace boundaryConditions {
      class ExtrapolateBoundaryCondition;
    }
  }
}

class peanoclaw::grid::boundaryConditions::ExtrapolateBoundaryCondition {
  public:
    void setBoundaryCondition(
      peanoclaw::Patch& subgrid,
      peanoclaw::grid::SubgridAccessor& accessor,
      int dimension,
      bool setUpper,
      const tarch::la::Vector<DIMENSIONS,int>& sourceSubcellIndex,
      const tarch::la::Vector<DIMENSIONS,int>& destinationSubcellIndex
    );
};

#endif /* PEANOCLAW_GRID_BOUNDARYCONDITIONS_EXTRAPOLATEBOUNDARYCONDITION_H_ */
