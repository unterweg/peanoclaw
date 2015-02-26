/*
 * GEBCO.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include "peanoclaw/native/scenarios/GEBCO.h"

#include "peano/utils/Loop.h"
#include "tarch/Assertions.h"

#include <limits>

tarch::logging::Log peanoclaw::native::scenarios::GEBCO::_log("peanoclaw::native::scenarios::GEBCO");

const double peanoclaw::native::scenarios::GEBCO::EARTH_RADIUS = 6378.137 * 1000.0;

#ifdef PEANOCLAW_USE_NETCDF
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

void peanoclaw::native::scenarios::GEBCO::readRange(NcFile& ncFile, const std::string& rangeName, double range[2]) {
  NcVar rangeVar=ncFile.getVar(rangeName);
  assertion1(!rangeVar.isNull, "No range of name " << rangeName << ".");
  rangeVar.getVar(range);
}
#endif

tarch::la::Vector<2, int> peanoclaw::native::scenarios::GEBCO::normalizeIndex(
  const tarch::la::Vector<2, int>& index
) const {
  tarch::la::Vector<2, int> normalizedIndex;
  for(int d = 0;  d < 2; d++) {
    normalizedIndex[d] = std::max(0, std::min(_numberOfCells[d]-1, index[d]));
  }
  return normalizedIndex;
}

tarch::la::Vector<2, int> peanoclaw::native::scenarios::GEBCO::mapToGEBCO(
  const tarch::la::Vector<2,double>& position
) const {
  tarch::la::Vector<2, int> cellIndex = tarch::la::multiplyComponents(
                                          position,
                                          tarch::la::invertEntries(_cellSize)
                                        ).convertScalar<int>();
  return normalizeIndex(cellIndex);
}

double peanoclaw::native::scenarios::GEBCO::getBathymetry(const tarch::la::Vector<2,int>& cellIndex) const {
  return _bathymetry[cellIndex[0] + (_numberOfCells[1] - 1 - cellIndex[1]) * _numberOfCells[0]];
}

peanoclaw::native::scenarios::GEBCO::GEBCO()
  : _minBathymetry(-1),
    _maxBathymetry(-1)
{
  #ifdef PEANOCLAW_USE_NETCDF
  try {
    std::string bathymetryFileName("bathymetry.nc");
    NcFile ncFile(bathymetryFileName, NcFile::read);

    double xRange[2];
    readRange(ncFile, "x_range", xRange);
    double yRange[2];
    readRange(ncFile, "y_range", yRange);
    double zRange[2];
    readRange(ncFile, "z_range", zRange);

    _offset[0] = 0;
    _offset[1] = 0;
    _coordinateOrigin[0] = xRange[0];
    _coordinateOrigin[1] = xRange[1];
    _size[0] = xRange[1] - xRange[0];
    _size[1] = yRange[1] - yRange[0];
    _size = mapLatitudeLongitudeToMeters(_size);
    _minBathymetry = zRange[0];
    _maxBathymetry = zRange[1];

    int dimensions[2];
    NcVar dimensionsVar = ncFile.getVar("dimension");
    assertion1(!dimensionsVar.isNull(), "No dimensions in file.");
    dimensionsVar.getVar(dimensions);
    _numberOfCells[0] = dimensions[0];
    _numberOfCells[1] = dimensions[1];
    _cellSize = tarch::la::multiplyComponents(_size, tarch::la::invertEntries(_numberOfCells.convertScalar<double>()));
//    _cellSize[0] = (xRange[1] - xRange[0]) / _numberOfCells[0];
//    _cellSize[1] = (yRange[1] - yRange[0]) / _numberOfCells[1];

    _bathymetry.resize(_numberOfCells[0] * _numberOfCells[1]);

    NcVar bathymetryVar = ncFile.getVar("z");
    assertion1(!bathymetryVar.isNull(), "No z data given in file.");
    bathymetryVar.getVar(_bathymetry.data());
  } catch(NcException& e) {
    logError("GEBCO", "Error while reading bathymetry: " << e.what());
  }
  #endif
}

tarch::la::Vector<2, double> peanoclaw::native::scenarios::GEBCO::getOffset() const {
  return _offset;
}

tarch::la::Vector<2, double> peanoclaw::native::scenarios::GEBCO::getSize() const {
  return _size;
}

double peanoclaw::native::scenarios::GEBCO::getBathymetry(
  const tarch::la::Vector<DIMENSIONS, double>& position
) const {
  #ifdef Dim2
  tarch::la::Vector<2, int> cellIndex = mapToGEBCO(position);
  return getBathymetry(cellIndex);
  #else
  return 0.0;
  #endif
}

double peanoclaw::native::scenarios::GEBCO::getAveragedBathymetry(
  const tarch::la::Vector<DIMENSIONS,double>& offset,
  const tarch::la::Vector<DIMENSIONS,double>& size
) const {
  #ifdef Dim2
  tarch::la::Vector<2,int> lowerBound = mapToGEBCO(offset);
  tarch::la::Vector<2,int> upperBound = mapToGEBCO(offset + size - tarch::la::Vector<2,double>(1e-8))
                                        + tarch::la::Vector<2,int>(1);

  double bathymetry = 0.0;
  int numberOfCells = 0;
  dfor(localCellIndex, upperBound - lowerBound) {
    tarch::la::Vector<2,int> cellIndex = localCellIndex + lowerBound;
    bathymetry += getBathymetry(cellIndex);
    numberOfCells++;
  }
  return bathymetry / numberOfCells;
  #else
  return 0.0;
  #endif
}

double peanoclaw::native::scenarios::GEBCO::getMaximalBathymetryDifference(
  const tarch::la::Vector<DIMENSIONS,double>& offset,
  const tarch::la::Vector<DIMENSIONS,double>& size
) const {
  double minBathymetry;
  double maxBathymetry;
  getMinAndMaxBathymetry(offset, size, minBathymetry, maxBathymetry);
  return maxBathymetry - minBathymetry;
}

void peanoclaw::native::scenarios::GEBCO::getMinAndMaxBathymetry(
  const tarch::la::Vector<DIMENSIONS,double>& offset,
  const tarch::la::Vector<DIMENSIONS,double>& size,
  double& minBathymetry,
  double& maxBathymetry
) const {
  tarch::la::Vector<2,int> lowerBound = mapToGEBCO(offset);
  tarch::la::Vector<2,int> upperBound = mapToGEBCO(offset + size - tarch::la::Vector<2,double>(1e-8))
                                        + tarch::la::Vector<2,int>(1);

  minBathymetry = std::numeric_limits<double>::max();
  maxBathymetry = -std::numeric_limits<double>::max();
  dfor(localCellIndex, upperBound - lowerBound) {
    tarch::la::Vector<2,int> cellIndex = localCellIndex + lowerBound;
    double bathymetry = getBathymetry(cellIndex);
    minBathymetry = std::min(minBathymetry, bathymetry);
    maxBathymetry = std::max(maxBathymetry, bathymetry);
  }
}

tarch::la::Vector<2, double> peanoclaw::native::scenarios::GEBCO::mapLatitudeLongitudeToMeters(
  const tarch::la::Vector<2,double>& coordinates
) const {
  tarch::la::Vector<2,double> radians = coordinates;
  for(int d = 0; d < 2; d++) {
    radians[d] *= M_PI / 180.0;
  }
  tarch::la::Vector<2,double> meters;
  meters[0] = 2 * M_PI * EARTH_RADIUS * (coordinates[0] / 360.0) * cos(coordinates[1] * M_PI / 180.0);
  meters[1] = 2 * M_PI * EARTH_RADIUS * (coordinates[1] / 360.0);
  return meters / 1000.0;
}

tarch::la::Vector<2, double> peanoclaw::native::scenarios::GEBCO::mapMetersToLatitudeLongitude(
  const tarch::la::Vector<2,double>& meters
) const {
  tarch::la::Vector<2,double> coordinates;
  coordinates[1] = meters[1] * (360.0 / (2 * M_PI * EARTH_RADIUS));
  coordinates[0] = meters[0] * (360.0 / (2 * M_PI * EARTH_RADIUS * cos(coordinates[1] * M_PI / 180.0)));
  return coordinates;
}
