/*
 * BathymetryTransformation.cpp
 *
 *  Created on: Oct 4, 2015
 *      Author: kristof
 */
#include "peanoclaw/interSubgridCommunication/BathymetryTransformation.h"

#include "peanoclaw/Patch.h"

#include "peano/utils/Loop.h"

peanoclaw::interSubgridCommunication::BathymetryTransformation::BathymetryTransformation(
  double seaLevel
) : _seaLevel(seaLevel) {

}

void peanoclaw::interSubgridCommunication::BathymetryTransformation::transformWaterHeight(
  peanoclaw::Patch& subgrid,
  const peanoclaw::geometry::Region&       region,
  bool              modifyUOld,
  bool              absoluteToAboveSeaFloor
) const {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();
  if(modifyUOld) {
    dfor(internalSubcellIndex, region._size) {
      tarch::la::Vector<DIMENSIONS,int> subcellIndex = internalSubcellIndex + region._offset;
      double value = transformValue(accessor.getValueUOld(subcellIndex, 0), accessor.getParameterWithGhostlayer(subcellIndex, 0), absoluteToAboveSeaFloor);
      accessor.setValueUOld(
        subcellIndex,
        0,
        value
      );
    }
  } else {
    dfor(internalSubcellIndex, region._size) {
      tarch::la::Vector<DIMENSIONS,int> subcellIndex = internalSubcellIndex + region._offset;
      double value = transformValue(accessor.getValueUNew(subcellIndex, 0), accessor.getParameterWithGhostlayer(subcellIndex, 0), absoluteToAboveSeaFloor);
      accessor.setValueUNew(
        subcellIndex,
        0,
        value
      );
    }
  }
}




