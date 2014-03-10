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

  FullSWOF2D_Parameters par(patch, maximumTimestepSize);
  //std::cout << "parameters read" << std::endl;

  Choice_scheme *wrapper_scheme = new Choice_scheme(par);
  Scheme *scheme = wrapper_scheme->getInternalScheme();

  // overwrite internal values
  tarch::la::Vector<DIMENSIONS,int> subdivisionFactor = patch.getSubdivisionFactor();
  tarch::la::Vector<DIMENSIONS,int> subcellIndex;

  // FullSWOF2D has a mixture of 0->nxcell+1 and 1->nxcell

  int ghostlayerWidth = patch.getGhostlayerWidth();
  int fullswofGhostlayerWidth = ghostlayerWidth - 1;

  /** Water height.*/
  TAB& h = scheme->getH();
  for (int x = -ghostlayerWidth; x < subdivisionFactor(0)+ghostlayerWidth; x++) {
        for (int y = -ghostlayerWidth; y < subdivisionFactor(1)+ghostlayerWidth; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            h[x+ghostlayerWidth][y+ghostlayerWidth] = patch.getValueUOld(subcellIndex, 0);
        }
  }
 
  /** X Velocity.*/
  TAB& u = scheme->getU();
  for (int x = -ghostlayerWidth; x < subdivisionFactor(0)+ghostlayerWidth; x++) {
        for (int y = -ghostlayerWidth; y < subdivisionFactor(1)+ghostlayerWidth; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            u[x+ghostlayerWidth][y+ghostlayerWidth] = patch.getValueUOld(subcellIndex, 1);
        }
  }

  /** Y Velocity.*/
  TAB& v = scheme->getV();
  for (int x = -ghostlayerWidth; x < subdivisionFactor(0)+ghostlayerWidth; x++) {
        for (int y = -ghostlayerWidth; y < subdivisionFactor(1)+ghostlayerWidth; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            v[x+ghostlayerWidth][y+ghostlayerWidth] = patch.getValueUOld(subcellIndex, 2);
        }
  }
 
  /** Topography.*/
  TAB& z = scheme->getZ();
  for (int x = -ghostlayerWidth; x < subdivisionFactor(0)+ghostlayerWidth; x++) {
        for (int y = -ghostlayerWidth; y < subdivisionFactor(1)+ghostlayerWidth; y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            z[x+ghostlayerWidth][y+ghostlayerWidth] = patch.getValueUOld(subcellIndex, 3);
        }
  }


  /** compute Discharge. (1->nxcell) */
  TAB& q1 = scheme->getQ1();
  for (int x = -(ghostlayerWidth-1); x < subdivisionFactor(0)+(ghostlayerWidth-1); x++) {
    for (int y = -(ghostlayerWidth-1); y < subdivisionFactor(1)+(ghostlayerWidth-1); y++) {
        subcellIndex(0) = x;
        subcellIndex(1) = y;
        q1[x+ghostlayerWidth][y+ghostlayerWidth] = patch.getValueUOld(subcellIndex, 4);
    }
  }

  /** compute Discharge. (1->nycell)*/
  TAB& q2 = scheme->getQ2();
  for (int x = -(ghostlayerWidth-1); x < subdivisionFactor(0)+(ghostlayerWidth-1); x++) {
    for (int y = -(ghostlayerWidth-1); y < subdivisionFactor(1)+(ghostlayerWidth-1); y++) {
        subcellIndex(0) = x;
        subcellIndex(1) = y;
        q2[x+ghostlayerWidth][y+ghostlayerWidth] = patch.getValueUOld(subcellIndex, 5);
    }
  }

  // kick off computation!
  scheme->setTimestep(maximumTimestepSize);
  do {
    scheme->resetN();
    scheme->calcul();
    if (scheme->getVerif() == 0) {
        std::cout << "scheme retry activated!" << std::endl;
    }
  } while (scheme->getVerif() == 0); // internal error detection of FullSWOF2D

  // copy back internal values but skip ghostlayer
  
  /** Water height after one step of the scheme.*/
  TAB& hs = scheme->getHs();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 0, hs[x+ghostlayerWidth][y+ghostlayerWidth]);
        }
  }
 
  /** X Velocity after one step of the scheme.*/
  TAB& us = scheme->getUs();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 1, us[x+ghostlayerWidth][y+ghostlayerWidth]);
        }
  }

  /** Y Velocity after one step of the scheme.*/
  TAB& vs = scheme->getVs();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 2, vs[x+ghostlayerWidth][y+ghostlayerWidth]);
        }
  }

  /** Topography.*/
  TAB& znew = scheme->getZ();
  for (int x = 0; x < subdivisionFactor(0); x++) {
        for (int y = 0; y < subdivisionFactor(1); y++) {
            subcellIndex(0) = x;
            subcellIndex(1) = y;
            patch.setValueUNew(subcellIndex, 3, znew[x+ghostlayerWidth][y+ghostlayerWidth]);
        }
  }
 
  /** compute Discharge. (1->nxcell) */
  TAB& qs1 = scheme->getQs1();
  for (int x = 0; x < subdivisionFactor(0); x++) {
    for (int y = 0; y < subdivisionFactor(1); y++) {
        subcellIndex(0) = x;
        subcellIndex(1) = y;
        patch.setValueUNew(subcellIndex, 4, qs1[x+ghostlayerWidth][y+ghostlayerWidth]);
    }
  }

  /** compute Discharge. (1->nycell)*/
  TAB& qs2 = scheme->getQs2();
  for (int x = 0; x < subdivisionFactor(0); x++) {
    for (int y = 0; y < subdivisionFactor(1); y++) {
        subcellIndex(0) = x;
        subcellIndex(1) = y;
        patch.setValueUNew(subcellIndex, 5, qs2[x+ghostlayerWidth][y+ghostlayerWidth]);
    }
  }

  double dt = min(scheme->getTimestep(), maximumTimestepSize);
  double estimatedNextTimestepSize = dt; //scheme->getTimestep();

  // BENCHMARK
  //dt = 0.00001;
  //estimatedNextTimestepSize = 0.00001;

  //std::cout << "\nComputation finished!" << endl;
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

#if 1
   // implement a wall boundary
    tarch::la::Vector<DIMENSIONS, int> src_subcellIndex;
    tarch::la::Vector<DIMENSIONS, int> dest_subcellIndex;

    if (dimension == 0) { // left and right boundary
        for (int yi = -2; yi < patch.getSubdivisionFactor()(1)+2; yi++) {
            int xi = setUpper ? patch.getSubdivisionFactor()(0)-1 : 0;
            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
            patch.setValueUOld(dest_subcellIndex, 3, 1000.0);

            // only mirror orthogonal velocities
            patch.setValueUOld(dest_subcellIndex, 1, -patch.getValueUOld(dest_subcellIndex, 2));
            patch.setValueUOld(dest_subcellIndex, 4, -patch.getValueUOld(dest_subcellIndex, 5));

            xi = setUpper ? patch.getSubdivisionFactor()(0)-2 : 1;
            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
            patch.setValueUOld(dest_subcellIndex, 3, 1000.0);
 
            // only mirror orthogonal velocities
            patch.setValueUOld(dest_subcellIndex, 1, -patch.getValueUOld(dest_subcellIndex, 2));
            patch.setValueUOld(dest_subcellIndex, 4, -patch.getValueUOld(dest_subcellIndex, 5));
        }
    } else { // top and bottom boundary
        for (int xi = -2; xi < patch.getSubdivisionFactor()(0)+2; xi++) {
            int yi = setUpper ? patch.getSubdivisionFactor()(1)-1 : 0;
            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
            patch.setValueUOld(dest_subcellIndex, 3, 1000.0);
            
            // only mirror orthogonal velocities
            patch.setValueUOld(dest_subcellIndex, 2, -patch.getValueUOld(dest_subcellIndex, 2));
            patch.setValueUOld(dest_subcellIndex, 5, -patch.getValueUOld(dest_subcellIndex, 5));

            yi = setUpper ? patch.getSubdivisionFactor()(1)-2 : 1;
            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
            patch.setValueUOld(dest_subcellIndex, 3, 1000.0);
 
            // only mirror orthogonal velocities
            patch.setValueUOld(dest_subcellIndex, 2, -patch.getValueUOld(dest_subcellIndex, 2));
            patch.setValueUOld(dest_subcellIndex, 5, -patch.getValueUOld(dest_subcellIndex, 5));
        }
    }
#else

    // implement an open boundary
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
                    patch.setValueUOld(dest_subcellIndex, unknown, 0.0);
                } else {
                    patch.setValueUOld(dest_subcellIndex, unknown, 0.0);
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
                    patch.setValueUOld(dest_subcellIndex, unknown, 0.0);
                } else {
                    patch.setValueUOld(dest_subcellIndex, unknown, 0.0);
                }
            }
        }
    }

#endif

   //std::cout << "++++++" << std::endl;
   //std::cout << patch.toStringUOldWithGhostLayer() << std::endl;
   //std::cout << "||||||" << std::endl;

      logTraceOut("fillBoundaryLayerInPyClaw");
}

peanoclaw::native::FullSWOF2D_Parameters::FullSWOF2D_Parameters(Patch& patch, double maximumTimestepSize) {
    // seed parameters based on Input file
    //setparameters("./fullswof2d_parameters.txt");
 
    int ghostlayerWidth = patch.getGhostlayerWidth();

    // now override the peanoclaw specific ones
    Nxcell = patch.getSubdivisionFactor()(0)+2*(ghostlayerWidth-1); // FullSWOF2D is not aware of a ghostlayer, so we cheat now!
    Nycell = patch.getSubdivisionFactor()(1)+2*(ghostlayerWidth-1); // FullSWOF2D is not aware of a ghostlayer, so we cheat now!

    scheme_type = 1; // 1= fixed cfl => we get timestamp and maximum timestep
                     // 2=fixed dt

    // TODO: we probably have to provide dx / dy values directly
    // FullSOWF2D uses L and l just to compute dx and dy
    L = patch.getSubcellSize()(0)*Nxcell; // length of domain in x direction (TODO: i multiply be Nxcell because the code will divide later on...)
    l = patch.getSubcellSize()(1)*Nycell; // length of domain in y direction (TODO: i multiply be Nycell because the code will divide later on...)

    //T = maximumTimestepSize; this is a bad idea: the algorithm internally works similar to peanoclaw: either dt caused by cfl or the remaining part in the time interval
    //(see line 198 of order2.cpp)
    T = 1000; // that should be enough as we only do one timestep anyway
 
    order = 2; // order 1 or order 2

   cfl_fix = 0.5; // TODO: what is this actually good for? (working setting 0.5 with order 1, 0.8 got stuck with order 1 and 2)
   dt_fix = 0.05; // TODO: what is this actually good for
   nbtimes = 1;
   dx = patch.getSubcellSize()(0);
   dy = patch.getSubcellSize()(1);

   // TODO: actually deacticated by macro
   Lbound = 3;
   Rbound = 3;
   Bbound = 3;
   Tbound = 3;
	
   flux = 2; // 1 = rusanov 2 = HLL 3 = HLL2
   rec = 1;
   fric = 2; // Friction law (0=NoFriction 1=Manning 2=Darcy-Weisbach)  <fric>:: 0
   lim = 1;
   inf = 1; // Infiltration model (0=No Infiltration 1=Green-Ampt)
   topo = 2; // flat topogrophy, just prevent it from loading data as we will fill the data in later on
   huv_init = 2; // initialize h,v and u to 0
   rain = 2; // 2: rain is generated, basically its just an auxillary array which is filled with time dependent data (we can couple this evolveToTime)
   amortENO = 0.25;
   modifENO = 0.9;
   //frotcoef is not actually used though it is existing in parameters
   friccoef = 0;

  //for infiltration model
  Kc_init = 2;
  Kc_coef = 1.8e-5;

  Ks_init = 2;
  Ks_coef = 0.0;

  dtheta_init = 2;
  dtheta_coef = 0.254;

  Psi_init = 2;
  Psi_coef = 0.167;

  zcrust_init = 2;
  zcrust_coef = 0.01;

  imax_init = 2;
  imax_coef = 5.7e-4;

  //  SCALAR  Kssoil;

  left_imp_discharge = 0.0;
  left_imp_h = 0.1;
  right_imp_discharge = 0.0;
  right_imp_h = 0.1;
  bottom_imp_discharge = 0.0;
  bottom_imp_h = 0.1;
  top_imp_discharge = 0.0;
  top_imp_h = 0.1;

  output_format = 0; // disable all output
}

peanoclaw::native::FullSWOF2D_Parameters::~FullSWOF2D_Parameters() {

}
