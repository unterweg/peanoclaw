/*
 * Storegga.h
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_NATIVE_SCENARIOS_STOREGGA_H_
#define PEANOCLAW_NATIVE_SCENARIOS_STOREGGA_H_

#include "peanoclaw/Patch.h"
#include "peanoclaw/native/scenarios/SWEScenario.h"
#include "peanoclaw/native/scenarios/GEBCO.h"

#include <vector>

namespace peanoclaw {
  namespace native {
    namespace scenarios {
      class Storegga;
    }
  }
}

class peanoclaw::native::scenarios::Storegga : public peanoclaw::native::scenarios::SWEScenario {
  private:
//    tarch::la::Vector<DIMENSIONS,double> _domainSize;
    peanoclaw::native::scenarios::GEBCO   _gebco;
    tarch::la::Vector<DIMENSIONS, double> _minimalMeshWidth;
    tarch::la::Vector<DIMENSIONS, double> _maximalMeshWidth;
    tarch::la::Vector<DIMENSIONS,int>    _subdivisionFactor;
    double                               _endTime;
    double                               _globalTimestepSize;

    tarch::la::Vector<DIMENSIONS, double> _depressionCenter;
    tarch::la::Vector<DIMENSIONS, double> _plumeCenter;
    double                                _slideRadius;
    double                                _slideVelocity;
    double                                _slideDepth;

  public:
    Storegga(
      const std::vector<std::string>& arguments
    );

    virtual void initializePatch(peanoclaw::Patch& patch);
    virtual tarch::la::Vector<DIMENSIONS,double> computeDemandedMeshWidth(peanoclaw::Patch& patch, bool isInitializing);

    virtual void update(peanoclaw::Patch& patch);

    virtual void setBoundaryCondition(
      peanoclaw::Patch& subgrid,
      peanoclaw::grid::SubgridAccessor& accessor,
      int dimension,
      bool setUpper,
      const tarch::la::Vector<DIMENSIONS,int>& sourceSubcellIndex,
      const tarch::la::Vector<DIMENSIONS,int>& destinationSubcellIndex
    );

    //PeanoClaw-Scenario
    tarch::la::Vector<DIMENSIONS,double> getDomainOffset() const;
    tarch::la::Vector<DIMENSIONS,double> getDomainSize() const;
    tarch::la::Vector<DIMENSIONS,double> getInitialMinimalMeshWidth() const;
    tarch::la::Vector<DIMENSIONS,int>    getSubdivisionFactor() const;
    double                               getGlobalTimestepSize() const;
    double                               getEndTime() const;

    //pure SWE-Scenario
    virtual float getWaterHeight(float x, float y);
    virtual float waterHeightAtRest();
    virtual float endSimulation();
};


#endif /* PEANOCLAW_NATIVE_SCENARIOS_STOREGGA_H_ */
