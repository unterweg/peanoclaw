/*
 * SWASHESParameters.h
 *
 *  Created on: Oct 17, 2014
 *      Author: kristof
 */

#ifndef PEANOCLAW_NATIVE_SCENARIOS_SWASHES_SWASHESPARAMETERS_H_
#define PEANOCLAW_NATIVE_SCENARIOS_SWASHES_SWASHESPARAMETERS_H_

#ifdef PEANOCLAW_SWASHES
#include "SWASHES/parameters.hpp"
#endif

#include <string>
#include <vector>

namespace peanoclaw {
  namespace native {
    namespace scenarios {
      namespace swashes {
        class SWASHESParameters;
      }
    }
  }
}

class peanoclaw::native::scenarios::swashes::SWASHESParameters
#ifdef PEANOCLAW_SWASHES
: public SWASHES::Parameters
#endif
{
  private:
    char**       _parameterCStrings;
    std::string* _parameterStrings;

    char** getFilledParameterStrings(int numberOfCellsX, int choice);

  public:
    SWASHESParameters(int numberOfCellsX, int choice);

    virtual ~SWASHESParameters();
};


#endif /* PEANOCLAW_NATIVE_SCENARIOS_SWASHES_SWASHESPARAMETERS_H_ */
