/*
 * GEBCO.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */
#include "peanoclaw/native/scenarios/GEBCO.h"

#include "tarch/Assertions.h"

tarch::logging::Log peanoclaw::native::scenarios::GEBCO::_log("peanoclaw::native::scenarios::GEBCO");

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

peanoclaw::native::scenarios::GEBCO::GEBCO() {
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

    _offset[0] = xRange[0];
    _offset[1] = yRange[1];
    _size[0] = xRange[1] - xRange[0];
    _size[1] = yRange[1] - yRange[0];
    _minBathymetry = zRange[0];
    _maxBathymetry = zRange[1];

    int dimensions[2];
    NcVar dimensionsVar = ncFile.getVar("dimension");
    assertion1(!dimensionsVar.isNull(), "No dimensions in file.");
    dimensionsVar.getVar(dimensions);
    _numberOfCells[0] = dimensions[0];
    _numberOfCells[1] = dimensions[1];
    _cellSize = tarch::la::multiplyComponents(_size, tarch::la::invertEntries(_numberOfCells.convertScalar<double>()));

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
  tarch::la::Vector<2, int> cellIndex = tarch::la::multiplyComponents(
                                           (position - _offset),
                                           tarch::la::invertEntries(_cellSize)
                                        ).convertScalar<int>();
  for(int d = 0; d < DIMENSIONS; d++){
    cellIndex[d] = std::min(_numberOfCells[d]-1, std::max(cellIndex[d], 0));
  }
  cellIndex[1] = _numberOfCells[1] - 1 - cellIndex[1];
  return _bathymetry[cellIndex[0] + cellIndex[1] * _numberOfCells[0]];
  #else
  return 0.0;
  #endif
}
