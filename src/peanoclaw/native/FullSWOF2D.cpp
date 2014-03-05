/*
 * PyClaw.cpp
 *
 *  Created on: Feb 18, 2012
 *      Author: kristof
 */

#include "peanoclaw/Patch.h"
#include "tarch/timing/Watch.h"
#include "tarch/parallel/Node.h"

#include "peanoclaw/native/FullSWOF2D.h"

#include "peanoclaw/native/FullSWOF_2D/Headers/libschemes/choice_scheme.hpp"

tarch::logging::Log peanoclaw::native::FullSWOF2D::_log("peanoclaw::native::FullSWOF2D");

peanoclaw::native::FullSWOF2D::FullSWOF2D(
  SWEKernelScenario& scenario,
  peanoclaw::interSubgridCommunication::Interpolation*  interpolation,
  peanoclaw::interSubgridCommunication::Restriction*    restriction,
  peanoclaw::interSubgridCommunication::FluxCorrection* fluxCorrection
) : Numerics(interpolation, restriction, fluxCorrection),
_scenario(scenario),
_totalSolverCallbackTime(0.0)
{
  //import_array();

    // TODO: manually set number of cells and so on ...
}

peanoclaw::native::FullSWOF2D::~FullSWOF2D()
{
}


double peanoclaw::native::FullSWOF2D::initializePatch(
  Patch& patch
) {
  logTraceIn( "initializePatch(...)");

  _scenario.initializePatch(patch);
  double demandedMeshWidth = _scenario.computeDemandedMeshWidth(patch);

  logTraceOutWith1Argument( "initializePatch(...)", demandedMeshWidth);
  return demandedMeshWidth;
}

double peanoclaw::native::FullSWOF2D::solveTimestep(Patch& patch, double maximumTimestepSize, bool useDimensionalSplitting) {
  logTraceInWith2Arguments( "solveTimestep(...)", maximumTimestepSize, useDimensionalSplitting);

  assertion2(tarch::la::greater(maximumTimestepSize, 0.0), "Timestepsize == 0 should be checked outside.", patch.getTimeIntervals().getMinimalNeighborTimeConstraint());

  tarch::timing::Watch pyclawWatch("", "", false);
  pyclawWatch.startTimer();
  double dtAndEstimatedNextDt[2];

  // kick off the computation here -----

  FullSWOF2D_Parameters par(patch);
  std::cout << "parameters read" << std::endl;

  Choice_scheme *wrapper_scheme = new Choice_scheme(par);
  Scheme *scheme = wrapper_scheme->getInternalScheme();
  std::cout << "scheme chosen" << std::endl;

  // overwrite internal values
  tarch::la::Vector<DIMENSIONS,int> subdivisionFactor = patch.getSubdivisionFactor();
  tarch::la::Vector<DIMENSIONS,int> subcellIndex;

  /** Water height.*/
  TAB& h = scheme->getH();
  for (int x = -1; x < subdivisionFactor(0)+1; x++) {
        for (int y = -1; y < subdivisionFactor(1)+1; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            h[x+1][y+1] = patch.getValueUOld(subcellIndex, 0);
        }
  }
 
  /** X Velocity.*/
  TAB& u = scheme->getU();
  for (int x = -1; x < subdivisionFactor(0)+1; x++) {
        for (int y = -1; y < subdivisionFactor(1)+1; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            u[x+1][y+1] = patch.getValueUOld(subcellIndex, 1);
        }
  }

  /** Y Velocity.*/
  TAB& v = scheme->getV();
  for (int x = -1; x < subdivisionFactor(0)+1; x++) {
        for (int y = -1; y < subdivisionFactor(1)+1; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            v[x+1][y+1] = patch.getValueUOld(subcellIndex, 2);
        }
  }

  /** Discharge.*/
  TAB& q1 = scheme->getQ1();
  for (int x = 0; x < subdivisionFactor(0)+1; x++) {
        for (int y = -1; y < subdivisionFactor(1)+1; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            // h * u
            //q1[x+1][y+1] = patch.getValueUOld(subcellIndex, 0) * patch.getValueUOld(subcellIndex, 1);
            q1[x+1][y+1] = patch.getValueUOld(subcellIndex, 3);
        }
  }


  /** Discharge.*/
  TAB& q2 = scheme->getQ2();
  for (int x = 0; x < subdivisionFactor(0)+1; x++) {
        for (int y = -1; y < subdivisionFactor(1)+1; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            // h * v
            //q2[x+1][y+1] = patch.getValueUOld(subcellIndex, 0) * patch.getValueUOld(subcellIndex, 2);
            q2[x+1][y+1] = patch.getValueUOld(subcellIndex, 4);
        }
  }

  /** Topography.*/
  TAB& z = scheme->getZ();
  for (int x = -1; x < subdivisionFactor(0)+1; x++) {
        for (int y = -1; y < subdivisionFactor(1)+1; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            z[x+1][y+1] = patch.getValueUOld(subcellIndex, 5);
        }
  }

  // kick off computation!
  scheme->setTimestep(maximumTimestepSize);
  scheme->calcul();

  // copy back internal values
  
  /** Water height after one step of the scheme.*/
  TAB& hs = scheme->getHs();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 0, hs[x+1][y+1]);
        }
  }
 
  /** X Velocity after one step of the scheme.*/
  TAB& us = scheme->getUs();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 1, us[x+1][y+1]);
        }
  }

  /** Y Velocity after one step of the scheme.*/
  TAB& vs = scheme->getVs();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 2, vs[x+1][y+1]);
        }
  }

  /** Discharge after one step of the scheme.*/
  TAB& qs1 = scheme->getQs1();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 3, qs1[x+1][y+1]);
        }
  }


  /** Discharge after one step of the scheme.*/
  TAB& qs2 = scheme->getQs2();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 4, qs2[x+1][y+1]);
        }
  }

  double dt = fmin(scheme->getTimestep(), maximumTimestepSize);
  double estimatedNextTimestepSize = scheme->getTimestep();

  std::cout << "\nComputation finished!" << endl;
  delete wrapper_scheme;
  // computation is done -> back to peanoclaw 

  pyclawWatch.stopTimer();
  _totalSolverCallbackTime += pyclawWatch.getCalendarTime();

  assertion4(
      tarch::la::greater(patch.getTimeIntervals().getTimestepSize(), 0.0)
      || tarch::la::greater(estimatedNextTimestepSize, 0.0)
      || tarch::la::equals(maximumTimestepSize, 0.0)
      || tarch::la::equals(patch.getTimeIntervals().getEstimatedNextTimestepSize(), 0.0),
      patch, maximumTimestepSize, estimatedNextTimestepSize, patch.toStringUNew());
  assertion(patch.getTimeIntervals().getTimestepSize() < std::numeric_limits<double>::infinity());

  if (tarch::la::greater(dt, 0.0)) {
    patch.getTimeIntervals().advanceInTime();
    patch.getTimeIntervals().setTimestepSize(dt);
  }
  patch.getTimeIntervals().setEstimatedNextTimestepSize(estimatedNextTimestepSize);

  logTraceOut( "solveTimestep(...)");


  return _scenario.computeDemandedMeshWidth(patch);
}

void peanoclaw::native::FullSWOF2D::addPatchToSolution(Patch& patch) {
}

void peanoclaw::native::FullSWOF2D::fillBoundaryLayer(Patch& patch, int dimension, bool setUpper) {
  logTraceInWith3Arguments("fillBoundaryLayerInPyClaw", patch, dimension, setUpper);

  logDebug("fillBoundaryLayerInPyClaw", "Setting left boundary for " << patch.getPosition() << ", dim=" << dimension << ", setUpper=" << setUpper);

   //std::cout << "------ setUpper " << setUpper << " dimension " << dimension << std::endl;
   //std::cout << patch << std::endl;
   //std::cout << "++++++" << std::endl;
   //std::cout << patch.toStringUOldWithGhostLayer() << std::endl;
   //std::cout << "||||||" << std::endl;

   // implement a wall boundary
    tarch::la::Vector<DIMENSIONS, int> src_subcellIndex;
    tarch::la::Vector<DIMENSIONS, int> dest_subcellIndex;

    if (dimension == 0) {
        for (int yi = -1; yi < patch.getSubdivisionFactor()(1)+1; yi++) {
            int xi = setUpper ? patch.getSubdivisionFactor()(0) : -1;
            src_subcellIndex(0) = xi;
            src_subcellIndex(1) = yi;
            src_subcellIndex(dimension) += setUpper ? -1 : +1; 

            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
     
            for (int unknown=0; unknown < patch.getUnknownsPerSubcell(); unknown++) {
                double q = patch.getValueUOld(src_subcellIndex, unknown);

                if (unknown == dimension + 1) {
                    patch.setValueUOld(dest_subcellIndex, unknown, -q);
                } else {
                    patch.setValueUOld(dest_subcellIndex, unknown, q);
                }
            }
        }

    } else {
        for (int xi = -1; xi < patch.getSubdivisionFactor()(0)+1; xi++) {
            int yi = setUpper ? patch.getSubdivisionFactor()(1) : -1;
            src_subcellIndex(0) = xi;
            src_subcellIndex(1) = yi;
            src_subcellIndex(dimension) += setUpper ? -1 : +1; 

            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
     
            for (int unknown=0; unknown < patch.getUnknownsPerSubcell(); unknown++) {
                double q = patch.getValueUOld(src_subcellIndex, unknown);

                if (unknown == dimension + 1) {
                    patch.setValueUOld(dest_subcellIndex, unknown, -q);
                } else {
                    patch.setValueUOld(dest_subcellIndex, unknown, q);
                }
            }
        }
    }

   //std::cout << "++++++" << std::endl;
   //std::cout << patch.toStringUOldWithGhostLayer() << std::endl;
   //std::cout << "||||||" << std::endl;

      logTraceOut("fillBoundaryLayerInPyClaw");
}

peanoclaw::native::FullSWOF2D_Parameters::FullSWOF2D_Parameters(Patch& patch) {
    // seed parameters based on Input file
    setparameters("./fullswof2d_parameters.txt");

    // now override the peanoclaw specific ones
    Nxcell = patch.getSubdivisionFactor()(0)+2; // FullSWOF2D is not aware of a ghostlayer, so we cheat now!
    Nycell = patch.getSubdivisionFactor()(1)+2; // FullSWOF2D is not aware of a ghostlayer, so we cheat now!

    //scheme_type = 1; // 1= fixed cfl => we get timestamp and maximum timestep
    //                 // 2=fixed dt

    // TODO: we probably have to provide dx / dy values directly
    // FullSOWF2D uses L and l just to compute dx and dy
    L = patch.getSubcellSize()(0)*Nxcell; // length of domain in x direction (TODO)
    l = patch.getSubcellSize()(1)*Nycell; // length of domain in y direction (TODO)


    //order = 1; // TODO: there is a second order one, but we have to port it to peanoclaw first
}

peanoclaw::native::FullSWOF2D_Parameters::~FullSWOF2D_Parameters() {

}
