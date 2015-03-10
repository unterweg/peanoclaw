/*
 * GEBCO.h
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_NATIVE_SCENARIOS_GEBCO_H_
#define PEANOCLAW_NATIVE_SCENARIOS_GEBCO_H_

#include "peanoclaw/geometry/Region.h"

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include <vector>

#ifdef PEANOCLAW_USE_NETCDF
#include <netcdf>
#endif

namespace peanoclaw {
  namespace native {
    namespace scenarios {
      class GEBCO;
    }
  }
}

class peanoclaw::native::scenarios::GEBCO {
  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    static const double EARTH_RADIUS;
    static const double SCALING;

    std::vector<double>          _bathymetry;
    tarch::la::Vector<2, double> _coordinateOrigin; //Lower-left corner for latitude-longitude
    tarch::la::Vector<2, double> _offset;           //Offset in meters (always zero)
    tarch::la::Vector<2, double> _size;             //Size of domain in meters
    double                       _minBathymetry;
    double                       _maxBathymetry;
    tarch::la::Vector<2, int>    _numberOfCells;
    tarch::la::Vector<2, double> _cellSize;         //GEBCO-cell size in meters

    #ifdef PEANOCLAW_USE_NETCDF
    void readRange(netCDF::NcFile& ncFile, const std::string& rangeName, double range[2]);
    #endif

    /**
     * Normalizes the index such that it points to a cell within the GEBCO file.
     */
    tarch::la::Vector<2, int> normalizeIndex(const tarch::la::Vector<2, int>& index) const;

    /**
     * Maps the given position to the cell in the GEBCO file.
     */
    tarch::la::Vector<2, int> mapToGEBCO(const tarch::la::Vector<2,double>& position) const;

    /**
     * Maps the given region to a selection of GEBCO cells.
     */
    void mapToGEBCO(
      const tarch::la::Vector<2,double>& offset,
      const tarch::la::Vector<2,double>& size,
      tarch::la::Vector<2,int>& lowerBound,
      tarch::la::Vector<2,int>& upperBound
    ) const;

    /**
     * Returns the bathymetry for the given cell.
     */
    double getBathymetry(const tarch::la::Vector<2,int>& cellIndex) const;

  public:
    GEBCO();

    tarch::la::Vector<2, double> getOffset() const;

    tarch::la::Vector<2, double> getSize() const;

    double getBathymetry(const tarch::la::Vector<DIMENSIONS, double>& position) const;

    /**
     * Averages over all cells of the GEBCO file that are covered by the given region.
     */
    double getAveragedBathymetry(
      const tarch::la::Vector<DIMENSIONS,double>& offset,
      const tarch::la::Vector<DIMENSIONS,double>& size
    ) const;

    /**
     * Estimates the maximal difference in bathymetry within the given region.
     */
    double getMaximalBathymetryDifference(
      const tarch::la::Vector<DIMENSIONS,double>& offset,
      const tarch::la::Vector<DIMENSIONS,double>& size
    ) const;

    /**
     * Returns the minimal and maximal bathymetry within the given region.
     */
    void getMinAndMaxBathymetry(
      const tarch::la::Vector<DIMENSIONS,double>& offset,
      const tarch::la::Vector<DIMENSIONS,double>& size,
      double& minBathymetry,
      double& maxBathymetry
    ) const;

    /**
     * Maps world coordinates to meters
     */
    tarch::la::Vector<2, double> mapLatitudeLongitudeToMeters(const tarch::la::Vector<2,double>& coordinates) const;

    /**
     * Maps meters to world coordinates
     */
    tarch::la::Vector<2, double> mapMetersToLatitudeLongitude(const tarch::la::Vector<2,double>& meters) const;
};


#endif /* PEANOCLAW_NATIVE_SCENARIOS_GEBCO_H_ */
