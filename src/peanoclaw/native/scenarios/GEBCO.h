/*
 * GEBCO.h
 *
 *  Created on: Feb 23, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_NATIVE_SCENARIOS_GEBCO_H_
#define PEANOCLAW_NATIVE_SCENARIOS_GEBCO_H_

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

    std::vector<double>          _bathymetry;
    tarch::la::Vector<2, double> _offset;
    tarch::la::Vector<2, double> _size;
    double                       _minBathymetry;
    double                       _maxBathymetry;
    tarch::la::Vector<2, int>    _numberOfCells;
    tarch::la::Vector<2, double> _cellSize;

    #ifdef PEANOCLAW_USE_NETCDF
    void readRange(netCDF::NcFile& ncFile, const std::string& rangeName, double range[2]);
    #endif

  public:
    GEBCO();

    tarch::la::Vector<2, double> getOffset() const;

    tarch::la::Vector<2, double> getSize() const;

    double getBathymetry(const tarch::la::Vector<DIMENSIONS, double>& position) const;
};


#endif /* PEANOCLAW_NATIVE_SCENARIOS_GEBCO_H_ */
