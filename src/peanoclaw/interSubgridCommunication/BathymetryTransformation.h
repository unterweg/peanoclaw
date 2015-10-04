/*
 * BathymetryTransformation.h
 *
 *  Created on: Oct 4, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYTRANSFORMATION_H_
#define PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYTRANSFORMATION_H_

#include "peanoclaw/geometry/Region.h"

namespace peanoclaw {
  class Patch;

  namespace interSubgridCommunication {
    class BathymetryTransformation;
  }
}

class peanoclaw::interSubgridCommunication::BathymetryTransformation {
  private:
    double _seaLevel;

    double transformValue(
      double waterHeight,
      double bathymetry,
      bool absoluteToAboveSeaFloor
    ) const {
      if(absoluteToAboveSeaFloor) {
        return std::max(0.0, waterHeight - bathymetry);
      } else {
        if(tarch::la::greater(waterHeight, 0.0)) {
          return bathymetry + waterHeight;
        } else {
          return _seaLevel;
        }
      }
    }

  public:
    BathymetryTransformation(double seaLevel);

    void transformWaterHeight(
      peanoclaw::Patch& subgrid,
      const peanoclaw::geometry::Region&       region,
      bool              modifyUOld,
      bool              absoluteToAboveSeaFloor
    ) const;
};


#endif /* PEANOCLAW_INTERSUBGRIDCOMMUNICATION_BATHYMETRYTRANSFORMATION_H_ */
