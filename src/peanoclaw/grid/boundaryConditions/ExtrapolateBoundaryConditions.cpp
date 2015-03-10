/*
 * ExtrapolateBoundaryConditions.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: kristof
 */
#include "peanoclaw/grid/boundaryConditions/ExtrapolateBoundaryCondition.h"

#include "peanoclaw/Patch.h"
#include "peanoclaw/grid/SubgridAccessor.h"

#include "tarch/la/Vector.h"

void peanoclaw::grid::boundaryConditions::ExtrapolateBoundaryCondition::setBoundaryCondition(
    peanoclaw::Patch& subgrid,
    peanoclaw::grid::SubgridAccessor& accessor,
    int dimension,
    bool setUpper,
    const tarch::la::Vector<DIMENSIONS,int>& sourceSubcellIndex,
    const tarch::la::Vector<DIMENSIONS,int>& destinationSubcellIndex
) {
   //Copy
   for(int unknown = 0; unknown < subgrid.getUnknownsPerSubcell(); unknown++) {
     accessor.setValueUOld(destinationSubcellIndex, unknown, accessor.getValueUOld(sourceSubcellIndex, unknown));
   }
}


