/*
 * Storegga.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */
#include "peanoclaw/native/scenarios/Storegga.h"

#include "peanoclaw/grid/aspects/BoundaryIterator.h"
#include "peanoclaw/grid/boundaryConditions/ExtrapolateBoundaryCondition.h"

#include "peano/utils/Loop.h"

using namespace peanoclaw::grid::boundaryConditions;

#ifdef PEANOCLAW_USE_NETCDF
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#endif

peanoclaw::native::scenarios::Storegga::Storegga(
  const std::vector<std::string>& arguments
) : _gebco() {
  #ifdef Dim2
  if(arguments.size() != 5) {
    std::cerr << "Expected arguments for Scenario 'storegga': finestSubgridTopology coarsestSubgridTopology subdivisionFactor endTime globalTimestepSize" << std::endl
        << "\tGot " << arguments.size() << " arguments." << std::endl;
    throw "";
  }

  double finestSubgridTopologyPerDimension = atof(arguments[0].c_str());
  _minimalMeshWidth = _gebco.getSize()/ finestSubgridTopologyPerDimension;

  double coarsestSubgridTopologyPerDimension = atof(arguments[1].c_str());
  _maximalMeshWidth = _gebco.getSize() / coarsestSubgridTopologyPerDimension;

  //TODO unterweg debug
  std::cout << "max: " << _maximalMeshWidth << " min: " << _minimalMeshWidth << " size: " << _gebco.getSize() << " coarsest: "
      << coarsestSubgridTopologyPerDimension << std::endl;

  _subdivisionFactor = tarch::la::Vector<DIMENSIONS,int>(atoi(arguments[2].c_str()));

  _endTime = atof(arguments[3].c_str());

  _globalTimestepSize = atof(arguments[4].c_str());
  #endif
}

void peanoclaw::native::scenarios::Storegga::initializePatch(peanoclaw::Patch& subgrid) {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();

  update(subgrid);

  dfor(subcellIndex, subgrid.getSubdivisionFactor()) {
    double bathymetry = accessor.getParameterWithGhostlayer(subcellIndex, 0);
    accessor.setValueUNew(subcellIndex, 0, bathymetry > 0.0 ? 0.0 : -bathymetry);
    accessor.setValueUNew(subcellIndex, 0, 1);
    accessor.setValueUNew(subcellIndex, 1, 0.0);
    accessor.setValueUNew(subcellIndex, 2, 0.0);
  }
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::Storegga::computeDemandedMeshWidth(
  peanoclaw::Patch& patch,
  bool isInitializing
) {
  return _maximalMeshWidth;
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::Storegga::getDomainOffset() const {
  #ifdef Dim2
  return _gebco.getOffset();
  #else
  return 0;
  #endif
}

void peanoclaw::native::scenarios::Storegga::update(peanoclaw::Patch& subgrid) {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();
  dfor(subcellIndex, subgrid.getSubdivisionFactor() + 2 * subgrid.getGhostlayerWidth()) {
    tarch::la::Vector<DIMENSIONS,double> cellPosition = subgrid.getSubcellCenter(subcellIndex - subgrid.getGhostlayerWidth());
    //accessor.setParameterWithGhostlayer(subcellIndex - subgrid.getGhostlayerWidth(), 0, _gebco.getBathymetry(cellPosition));
    accessor.setParameterWithGhostlayer(subcellIndex - subgrid.getGhostlayerWidth(), 0, 0);
  }
}


void peanoclaw::native::scenarios::Storegga::setBoundaryCondition(
  peanoclaw::Patch& subgrid,
  peanoclaw::grid::SubgridAccessor& accessor,
  int dimension,
  bool setUpper,
  const tarch::la::Vector<DIMENSIONS,int>& sourceSubcellIndex,
  const tarch::la::Vector<DIMENSIONS,int>& destinationSubcellIndex
) {
  for(int unknown = 0; unknown < subgrid.getUnknownsPerSubcell(); unknown++) {
    accessor.setValueUOld(destinationSubcellIndex, unknown, accessor.getValueUOld(sourceSubcellIndex, unknown));
  }
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::Storegga::getDomainSize() const {
  #ifdef Dim2
  return _gebco.getSize();
  #else
  return 0;
  #endif
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::Storegga::getInitialMinimalMeshWidth() const {
  return _maximalMeshWidth;
}

tarch::la::Vector<DIMENSIONS,int> peanoclaw::native::scenarios::Storegga::getSubdivisionFactor() const {
  return _subdivisionFactor;
}

double peanoclaw::native::scenarios::Storegga::getGlobalTimestepSize() const {
  return _globalTimestepSize;
}

double peanoclaw::native::scenarios::Storegga::getEndTime() const {
  return _endTime;
}

float peanoclaw::native::scenarios::Storegga::getWaterHeight(float x, float y) {
  return 0.0;
}

float peanoclaw::native::scenarios::Storegga::waterHeightAtRest() {
  return 0.0;
}

float peanoclaw::native::scenarios::Storegga::endSimulation() {
  return _endTime;
}

