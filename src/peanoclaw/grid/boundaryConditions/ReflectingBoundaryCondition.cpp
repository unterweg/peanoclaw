/*
 * ReflectingBoundaryConditions.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: kristof
 */
#include "peanoclaw/grid/boundaryConditions/ReflectingBoundaryCondition.h"

#include "peanoclaw/Patch.h"
#include "peanoclaw/grid/SubgridAccessor.h"

void peanoclaw::grid::boundaryConditions::ReflectingBoundaryCondition::setBoundaryCondition(
    peanoclaw::Patch& subgrid,
    peanoclaw::grid::SubgridAccessor& accessor,
    int dimension,
    bool setUpper,
    const tarch::la::Vector<DIMENSIONS,int>& sourceSubcellIndex,
    const tarch::la::Vector<DIMENSIONS,int>& destinationSubcellIndex
) {
   for(int unknown = 0; unknown < subgrid.getUnknownsPerSubcell(); unknown++) {
     if(unknown == 1 + dimension) {
       // Reflect
       accessor.setValueUOld(destinationSubcellIndex, unknown, -accessor.getValueUOld(sourceSubcellIndex, unknown));
     } else {
       //Copy
       accessor.setValueUOld(destinationSubcellIndex, unknown, accessor.getValueUOld(sourceSubcellIndex, unknown));
     }
   }
}



