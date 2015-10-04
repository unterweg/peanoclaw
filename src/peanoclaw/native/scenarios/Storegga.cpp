/*
 * Storegga.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */
#include "peanoclaw/native/scenarios/Storegga.h"

#include "peanoclaw/geometry/Region.h"
#include "peanoclaw/grid/aspects/BoundaryIterator.h"
#include "peanoclaw/grid/boundaryConditions/ExtrapolateBoundaryCondition.h"
#include "peanoclaw/Numerics.h"

#include "peanoclaw/interSubgridCommunication/BathymetryInterpolation.h"
#include "peanoclaw/interSubgridCommunication/BathymetryRestriction.h"

#include "peano/utils/Loop.h"

#include <limits>

using namespace peanoclaw::grid::boundaryConditions;

#ifdef PEANOCLAW_USE_NETCDF
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#endif

peanoclaw::native::scenarios::Storegga::Storegga(
  const std::vector<std::string>& arguments
) : _gebco(),
    _slideRadius(150),
    _slideVelocity(15/5),
    _slideDepth(25) {
  #ifdef Dim2
  assignList(_depressionCenter) = -3, 64.5;
  _depressionCenter = _gebco.mapLatitudeLongitudeToMeters(_depressionCenter);
  assignList(_plumeCenter) = -4, 65.5;
  _plumeCenter = _gebco.mapLatitudeLongitudeToMeters(_plumeCenter);

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
      << coarsestSubgridTopologyPerDimension << " finest: " << finestSubgridTopologyPerDimension << " depressionCenter=" << _depressionCenter << std::endl;

  _subdivisionFactor = tarch::la::Vector<DIMENSIONS,int>(atoi(arguments[2].c_str()));

  _endTime = atof(arguments[3].c_str());

  _globalTimestepSize = atof(arguments[4].c_str());
  #endif
}

void peanoclaw::native::scenarios::Storegga::initializePatch(peanoclaw::Patch& subgrid) {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();

  update(subgrid);

  dfor(subcellIndex, subgrid.getSubdivisionFactor()) {
    tarch::la::Vector<DIMENSIONS,double> cellPosition = subgrid.getSubcellCenter(subcellIndex- subgrid.getGhostlayerWidth());

    double depressionRadius = tarch::la::norm2(cellPosition - _depressionCenter);
    double plumeRadius = tarch::la::norm2(cellPosition - _plumeCenter);
    double depressionFactor = sqrt(std::max(0.0, (_slideRadius - depressionRadius) / _slideRadius));
    double plumeFactor = sqrt(std::max(0.0, (_slideRadius - plumeRadius) / _slideRadius));
    double waterSurface = -depressionFactor * _slideDepth + plumeFactor * _slideDepth;

    double bathymetry = accessor.getParameterWithGhostlayer(subcellIndex, 0);
    double waterHeight = bathymetry > waterSurface ? 0.0 : waterSurface - bathymetry;
    accessor.setValueUNew(subcellIndex, 0, waterHeight);
    accessor.setValueUNew(subcellIndex, 1, -(depressionFactor + plumeFactor) * _slideVelocity * waterHeight);
    accessor.setValueUNew(subcellIndex, 2, (depressionFactor + plumeFactor) * _slideVelocity * waterHeight);
  }
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::Storegga::computeDemandedMeshWidth(
  peanoclaw::Patch& subgrid,
  bool isInitializing
) {
  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();

  double minBathymetry;
  double maxBathymetry;
  _gebco.getMinAndMaxBathymetry(
    subgrid.getPosition(),
    subgrid.getSize(),
    minBathymetry,
    maxBathymetry
  );

  double minWaterHeight = std::numeric_limits<double>::max();
  double maxWaterHeight = -std::numeric_limits<double>::max();
  dfor(subcellIndex, subgrid.getSubdivisionFactor()) {
    double absoluteWaterHeight = std::min(0.0, accessor.getParameterWithGhostlayer(subcellIndex, 0)) + accessor.getValueUNew(subcellIndex, 0);
    minWaterHeight = std::min(minWaterHeight, absoluteWaterHeight);
    maxWaterHeight = std::max(maxWaterHeight, absoluteWaterHeight);
  }

  if(/*minBathymetry * maxBathymetry < 0.0 ||*/ (maxWaterHeight - minWaterHeight) > 0.5) {
    //Refine along coastline
    return _minimalMeshWidth;
  } else if (maxWaterHeight - minWaterHeight > 1.0) {
    return subgrid.getSubcellSize();
  } else {
    return _maximalMeshWidth;
  }
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
  dfor(localSubcellIndex, subgrid.getSubdivisionFactor() + 2 * subgrid.getGhostlayerWidth()) {
    tarch::la::Vector<DIMENSIONS,int> subcellIndex = localSubcellIndex - tarch::la::Vector<DIMENSIONS,int>(subgrid.getGhostlayerWidth());
    tarch::la::Vector<DIMENSIONS,double> cellPosition = subgrid.getSubcellCenter(subcellIndex);
//    accessor.setParameterWithGhostlayer(subcellIndex, 0, _gebco.getBathymetry(cellPosition));

    tarch::la::Vector<DIMENSIONS,double> cellSize = subgrid.getSubcellSize();
    tarch::la::Vector<DIMENSIONS,double> cellOffset = cellPosition - cellSize / 2.0;
    accessor.setParameterWithGhostlayer(subcellIndex, 0, _gebco.getAveragedBathymetry(cellOffset, cellSize));
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
  double bathymetry = accessor.getParameterWithGhostlayer(destinationSubcellIndex, 0);
  accessor.setParameterWithGhostlayer(destinationSubcellIndex, 0, bathymetry);

  accessor.setValueUOld(destinationSubcellIndex, 0, bathymetry < 0.0 ? -bathymetry : 0.0);

  for(int unknown = 1; unknown < subgrid.getUnknownsPerSubcell(); unknown++) {
//    double sourceHeight = accessor.getValueUOld(sourceSubcellIndex, 0);
//    double destinationHeight = accessor.getValueUOld(destinationSubcellIndex, 0);

//    if(sourceHeight < 1) {
//      accessor.setValueUOld(destinationSubcellIndex, unknown, 0.0);
//    } else {
//      accessor.setValueUOld(destinationSubcellIndex, unknown, accessor.getValueUOld(sourceSubcellIndex, unknown) * (destinationHeight/sourceHeight));
//    }
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

peanoclaw::interSubgridCommunication::Interpolation* peanoclaw::native::scenarios::Storegga::getCustomInterpolation() const {
  return new peanoclaw::interSubgridCommunication::BathymetryInterpolation(0.0);
}

peanoclaw::interSubgridCommunication::Restriction* peanoclaw::native::scenarios::Storegga::getCustomRestriction() const {
  return new peanoclaw::interSubgridCommunication::BathymetryRestriction(0.0);
}

