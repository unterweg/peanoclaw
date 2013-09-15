/*
 * PeanoClawLibraryRunner.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: Kristof Unterweger
 */

#include "peanoclaw/runners/PeanoClawLibraryRunner.h"

#include "peanoclaw/runners/Runner.h"
#include "peanoclaw/repositories/RepositoryFactory.h"
#include "peanoclaw/State.h"
#include "peanoclaw/configurations/PeanoClawConfigurationForSpacetreeGrid.h"
#include "peanoclaw/records/CellDescription.h"
#include "peanoclaw/records/VertexDescription.h"
#include "peanoclaw/records/Data.h"
#include "peanoclaw/records/Cell.h"
#include "peanoclaw/records/Vertex.h"
#include "peanoclaw/statistics/LevelStatistics.h"

#include "peano/heap/Heap.h"
#include "peano/utils/UserInterface.h"
#include "peano/peano.h"

#if defined(Parallel)
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/messages/ForkMessage.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"
#include "mpibalancing/OracleForOnePhaseControlLoopWrapper.h"
#endif

#ifdef SharedTBB
#include "tarch/multicore/tbb/Core.h"
#endif

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

tarch::logging::Log peanoclaw::runners::PeanoClawLibraryRunner::_log("peanoclaw::runners::PeanoClawLibraryRunner");

peanoclaw::runners::PeanoClawLibraryRunner::PeanoClawLibraryRunner(
  peanoclaw::configurations::PeanoClawConfigurationForSpacetreeGrid& configuration,
  peanoclaw::Numerics& numerics,
  const tarch::la::Vector<DIMENSIONS, double>& domainOffset,
  const tarch::la::Vector<DIMENSIONS, double>& domainSize,
  const tarch::la::Vector<DIMENSIONS, double>& initialMinimalMeshWidth,
  const tarch::la::Vector<DIMENSIONS, int>& subdivisionFactor,
  int defaultGhostLayerWidth,
  int unknownsPerSubcell,
  int auxiliarFieldsPerSubcell,
  double initialTimestepSize,
  bool useDimensionalSplittingOptimization
) :
  _plotNumber(1),
  _configuration(configuration),
  _iterationTimer("peanoclaw::runners::PeanoClawLibraryRunner", "iteration", false),
  _totalRuntime(0.0),
  _numerics(numerics),
  _validateGrid(false)
{
  #ifndef Asserts
  _validateGrid = false;
  #endif

  //Parallel configuration
  #ifdef Parallel
  //tarch::parallel::Node::getInstance().setTimeOutWarning(10);
  //tarch::parallel::Node::getInstance().setDeadlockTimeOut(20);
  #endif

  //Multicore configuration
  #ifdef SharedMemoryParallelisation
  tarch::multicore::tbb::Core::getInstance().configure(1);
  #endif

  //User interface
  peano::utils::UserInterface userInterface;
  userInterface.writeHeader();

  //Initialize heap data
  peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().setName("CellDescription");
  peano::heap::Heap<peanoclaw::records::Data>::getInstance().setName("Data");
  peano::heap::Heap<peanoclaw::records::VertexDescription>::getInstance().setName("VertexDescription");
  peano::heap::Heap<peanoclaw::statistics::LevelStatistics>::getInstance().setName("LevelStatistics");
 
  initializeParallelEnvironment();

#ifdef SharedMemoryParallelisation
  peano::datatraversal::autotuning::Oracle::getInstance().setOracle( new peano::datatraversal::autotuning::OracleForOnePhaseDummy(
    true, // multithreading
    false,
    1, // splitTheThree
    false, // pipelineDescendProcessing
    false, // pipelineAscendProcessing
    tarch::la::aPowI(DIMENSIONS,3*3*3*3/2), // smallestGrainSizeForAscendDescend
    3, // grainSizeForAsendDescend
    tarch::la::aPowI(DIMENSIONS,9/2), // smallestGrainSizeForEnterLeaveCell
    2, // grainSizeForEnterLevelCell
    tarch::la::aPowI(DIMENSIONS,3*3*3*3+1), // smallestGrainSizeForTouchFirstLast
    64, // grainSizeForTouchFirstLast
    tarch::la::aPowI(DIMENSIONS,3*3*3), // smallestGrainSizeForSplitLoadStore
    8, // grainSizeForSplitLoadStore
    -1, // adapterNumber
   peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling // methodTrace
    )
  );
  /*peano::datatraversal::autotuning::Oracle::getInstance().setOracle( 
          new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize()
          //new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize( peano::datatraversal::autotuning::CallEnterCellOnRegularStationaryGrid  )
  );*/
#endif

  //tarch::la::Vector<DIMENSIONS,double> boundingBoxOffset(-1.0/7.0);
  tarch::la::Vector<DIMENSIONS,double> boundingBoxOffset = domainOffset;
  tarch::la::Vector<DIMENSIONS,double> boundingBoxSize = domainSize;
  boundingBoxSize *= 1.0; // (9.0/7.0);

  //Initialize pseudo geometry
  _geometry =
    new  peano::geometry::Hexahedron (
      domainSize,  // width
      domainOffset // offset
    );
  _repository =
    peanoclaw::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      *_geometry,
      domainSize,   // domainSize,
      domainOffset  // computationalDomainOffset
    );

  logInfo("PeanoClawLibraryRunner", "Initial values: "
      << "Domain size = [" << domainSize << "]"
      << ", default subdivision factor = " << subdivisionFactor
      << ", default ghostlayer width = " << defaultGhostLayerWidth
      << ", unknowns per cell = " << unknownsPerSubcell);

  assertion(_repository != 0);
  State& state = _repository->getState();
  state.setDefaultSubdivisionFactor(subdivisionFactor);
  state.setDefaultGhostLayerWidth(defaultGhostLayerWidth);
  state.setUnknownsPerSubcell(unknownsPerSubcell);
  state.setAuxiliarFieldsPerSubcell(auxiliarFieldsPerSubcell);
  //tarch::la::Vector<DIMENSIONS, double> initialMinimalSubcellSize = tarch::la::multiplyComponents(initialMinimalMeshWidth, subdivisionFactor.convertScalar<double>());
  tarch::la::Vector<DIMENSIONS, double> initialMinimalSubcellSize = initialMinimalMeshWidth;

  state.setInitialMinimalMeshWidth(initialMinimalSubcellSize);
  state.setNumerics(numerics);
  state.resetTotalNumberOfCellUpdates();
  state.setInitialTimestepSize(initialTimestepSize);
  state.setDomain(domainOffset, domainSize);
  state.setUseDimensionalSplittingOptimization(useDimensionalSplittingOptimization && !_configuration.disableDimensionalSplittingOptimization());

#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
     tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();  
  }
#endif

  tarch::la::Vector<DIMENSIONS, double> current_initialMinimalSubcellSize(0.1);
  tarch::la::Vector<DIMENSIONS, double> next_initialMinimalSubcellSize(0.1);
  
  if (tarch::la::oneGreater(next_initialMinimalSubcellSize,initialMinimalSubcellSize)) {
        current_initialMinimalSubcellSize = next_initialMinimalSubcellSize;
    } else {
        current_initialMinimalSubcellSize = initialMinimalSubcellSize;
    }

#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
 
#if 1
    state.setIsInitializing(true);
    while (tarch::la::oneGreater(current_initialMinimalSubcellSize,initialMinimalSubcellSize) || initialMinimalSubcellSize == current_initialMinimalSubcellSize) {
        state.setInitialMinimalMeshWidth(current_initialMinimalSubcellSize);

        std::cout << "current initial minimal subcell size: " << current_initialMinimalSubcellSize << std::endl;
  
        bool hasChangedState = false;

       hasChangedState = false;
       if(_validateGrid) {
            _repository->switchToInitialiseAndValidateGrid();
        } else {
            _repository->switchToInitialiseGrid();
        }
        _repository->iterate();
        hasChangedState |= (!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced());
        //std::cout << "first iteration block (1)" << _repository->getState() << std::endl;

        _repository->switchToRemesh();

        do {
            // TODO: not sure if we need all 3 iterations of them or if it is sufficient to check for a stationary grid

            hasChangedState = false;

            _repository->iterate();
            //hasChangedState |= !_repository->getState().isGridStationary();
            hasChangedState |= (!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced());
            std::cout << "first iteration block (1)" << _repository->getState() << " | working_nodes " << tarch::parallel::NodePool::getInstance().getNumberOfWorkingNodes() << std::endl;

           /* _repository->iterate();
            //hasChangedState |= !_repository->getState().isGridStationary();
            hasChangedState |= (!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced());
            std::cout << "first iteration block (2)" << _repository->getState() << " | working_nodes " << tarch::parallel::NodePool::getInstance().getNumberOfWorkingNodes() << std::endl;

            _repository->iterate();
            //hasChangedState |= !_repository->getState().isGridStationary();
             hasChangedState |= (!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced());
            std::cout << "first iteration block (3)" << _repository->getState() << " | working_nodes " << tarch::parallel::NodePool::getInstance().getNumberOfWorkingNodes() << std::endl;*/
        } while (hasChangedState);
 
        next_initialMinimalSubcellSize = current_initialMinimalSubcellSize * (1.0/subdivisionFactor[0]);
        std::cout << "next initial minimal subcell size " << next_initialMinimalSubcellSize 
                    << " " << current_initialMinimalSubcellSize 
                    << " " << initialMinimalSubcellSize 
                    << std::endl;

        if (tarch::la::oneGreater(next_initialMinimalSubcellSize,initialMinimalSubcellSize)) {
            current_initialMinimalSubcellSize = next_initialMinimalSubcellSize;
        } else if (initialMinimalSubcellSize == current_initialMinimalSubcellSize) {
            current_initialMinimalSubcellSize *= (1.0/3.0); // abort loop
        } else {
            current_initialMinimalSubcellSize = initialMinimalSubcellSize;
        }
    }
#endif

#if 0 // works
      state.setIsInitializing(true);
      state.setInitialMinimalMeshWidth(initialMinimalSubcellSize);
  
        bool hasChangedState = false;

        do {
           hasChangedState = false;
           if(_validateGrid) {
                _repository->switchToInitialiseAndValidateGrid();
            } else {
                _repository->switchToInitialiseGrid();
            }
            _repository->iterate();
            hasChangedState |= !_repository->getState().isGridStationary();
            std::cout << "first iteration block (1)" << _repository->getState() << std::endl;

            _repository->switchToRemesh();

            _repository->iterate();
            hasChangedState |= !_repository->getState().isGridStationary();
            std::cout << "first iteration block (1)" << _repository->getState() << std::endl;
 
            _repository->iterate();
            hasChangedState |= !_repository->getState().isGridStationary();
            std::cout << "first iteration block (1)" << _repository->getState() << std::endl;
 
            _repository->iterate();
            hasChangedState |= !_repository->getState().isGridStationary();
            std::cout << "first iteration block (1)" << _repository->getState() << std::endl;
        } while (hasChangedState);
#endif

      do {
        state.setInitialRefinementTriggered(false);
        if(_validateGrid) {
          _repository->switchToInitialiseAndValidateGrid();
        } else {
          _repository->switchToInitialiseGrid();
        }
        std::cout << "second iteration block (1)" << _repository->getState() << std::endl;
        _repository->iterate();
      } while(state.getInitialRefinementTriggered());
      state.setIsInitializing(false);

      _repository->getState().setPlotNumber(0);
      if(_configuration.plotAtOutputTimes() || _configuration.plotSubsteps()) {
        _repository->switchToPlot(); _repository->iterate();
      }
  }
#endif

      //std::cout << "done with step by step grid construction: final initial Meshwidth " << state.getInitialMinimalMeshWidth() << std::endl;

#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      _repository->getState().setPlotNumber(0);
      if(_configuration.plotAtOutputTimes() || _configuration.plotSubsteps()) {
        _repository->switchToPlot(); _repository->iterate();
      }
  }
#endif
}

peanoclaw::runners::PeanoClawLibraryRunner::~PeanoClawLibraryRunner()
{
  logTraceIn("~PeanoClawLibraryRunner");
//  logInfo("~PeanoClawLibraryRunner()", "Total number of cell updates: " << _repository->getState().getTotalNumberOfCellUpdates());

  //Delete remaining heap data
  peano::heap::Heap<peanoclaw::records::CellDescription>& cellDescriptionHeap = peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance();
//  if(cellDescriptionHeap.getNumberOfAllocatedEntries() > 0) {
//    logWarning("~PeanoClawLibraryRunner()", "The heap for CellDescriptions still contains " << cellDescriptionHeap.getNumberOfAllocatedEntries() << " undeleted entries.");
//  }
  peano::heap::Heap<peanoclaw::records::Data>& dataHeap = peano::heap::Heap<peanoclaw::records::Data>::getInstance();
//  if(dataHeap.getNumberOfAllocatedEntries() > 0) {
//    logWarning("~PeanoClawLibraryRunner()", "The heap for patch data still contains " << dataHeap.getNumberOfAllocatedEntries() << " undeleted entries.");
//  }
  peano::heap::Heap<peanoclaw::records::VertexDescription>& vertexDescriptionHeap = peano::heap::Heap<peanoclaw::records::VertexDescription>::getInstance();
//  if(vertexDescriptionHeap.getNumberOfAllocatedEntries() > 0) {
//    logWarning("~PeanoClawLibraryRunner()", "The heap for VertexDescriptions still contains " << vertexDescriptionHeap.getNumberOfAllocatedEntries() << " undeleted entries.");
//  }
  cellDescriptionHeap.deleteAllData();
  dataHeap.deleteAllData();
  vertexDescriptionHeap.deleteAllData();

  #ifdef Parallel
  if(tarch::parallel::Node::getInstance().isGlobalMaster()) {
  #endif
    _repository->getState().plotTotalStatistics();

    _repository->logIterationStatistics();
    _repository->terminate();
  #ifdef Parallel
  }
  #endif
  delete _repository;
  delete _geometry;

  cellDescriptionHeap.restart();
  dataHeap.restart();
  vertexDescriptionHeap.restart();

  #ifdef Parallel
  //tarch::parallel::NodePool::getInstance().terminate();
  #endif
  //peano::shutdownParallelEnvironment();
  //peano::shutdownSharedMemoryEnvironment();

  logTraceOut("~PeanoClawLibraryRunner");
}

void peanoclaw::runners::PeanoClawLibraryRunner::evolveToTime(
  double time
) {
  logTraceIn("evolveToTime");
  
  bool plotSubsteps = _configuration.plotSubsteps()
      || (_configuration.plotSubstepsAfterOutputTime() != -1 && _configuration.plotSubstepsAfterOutputTime() <= _plotNumber);
 
  configureGlobalTimestep(time);
  do {
      runNextPossibleTimestep();
  } while(!_repository->getState().getAllPatchesEvolvedToGlobalTimestep());

  if(_configuration.plotAtOutputTimes() && !plotSubsteps) {
    _repository->getState().setPlotNumber(_plotNumber);
    if(_validateGrid) {
      _repository->switchToPlotAndValidateGrid(); _repository->iterate();
    } else {
      _repository->switchToPlot(); _repository->iterate();
    }
    _plotNumber++;
  } else if (!_configuration.plotAtOutputTimes() && !plotSubsteps) {
    _plotNumber++;
  }
  logTraceOut("evolveToTime");
}

void peanoclaw::runners::PeanoClawLibraryRunner::initializeParallelEnvironment() {
 /* #if defined(Parallel)
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().setStrategy( new tarch::parallel::FCFSNodePoolStrategy() );
  }
  tarch::parallel::NodePool::getInstance().restart();

  peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
    new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(true)
  );

  // have to be the same for all ranks
  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(2048);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(2048);
  #endif*/
 
  //tarch::parallel::NodePool::getInstance().restart();
  
  peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
        //new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(true)
        new mpibalancing::OracleForOnePhaseControlLoopWrapper(true, _controlLoopLoadBalancer)
  );
}

void peanoclaw::runners::PeanoClawLibraryRunner::gatherCurrentSolution() {
  logTraceIn("gatherCurrentSolution");
  assertion(_repository != 0);

  if(_validateGrid) {
    _repository->switchToGatherCurrentSolutionAndValidateGrid(); _repository->iterate();
  } else {
    _repository->switchToGatherCurrentSolution(); _repository->iterate();
  }
  logTraceOut("gatherCurrentSolution");
}

int peanoclaw::runners::PeanoClawLibraryRunner::runWorker() {
  #ifdef Parallel
  int newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob();
  while ( newMasterNode != tarch::parallel::NodePool::JobRequestMessageAnswerValues::Terminate ) {
    if ( newMasterNode >= tarch::parallel::NodePool::JobRequestMessageAnswerValues::NewMaster ) {
      peano::parallel::messages::ForkMessage forkMessage;
      forkMessage.receive(tarch::parallel::NodePool::getInstance().getMasterRank(),tarch::parallel::NodePool::getInstance().getTagForForkMessages(), true);
      _repository->restart(
        forkMessage.getH(),
        forkMessage.getDomainOffset(),
        forkMessage.getLevel(),
        forkMessage.getPositionOfFineGridCellRelativeToCoarseGridCell()
      );

      bool continueToIterate = true;
      while (continueToIterate) {
        switch (_repository->continueToIterate()) {
          case peanoclaw::repositories::Repository::Continue:
            _repository->iterate();
            break;
          case peanoclaw::repositories::Repository::Terminate:
            continueToIterate = false;
            break;
          case peanoclaw::repositories::Repository::RunGlobalStep:
            //runGlobalStep();
            break;
        }
      }

      // insert your postprocessing here
      // -------------------------------

      // -------------------------------

      _repository->terminate();
    }
    else if ( newMasterNode == tarch::parallel::NodePool::JobRequestMessageAnswerValues::RunAllNodes ) {
      //runGlobalStep();
    }
    newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob();
  }

  #endif
  return 0;
}

const peanoclaw::State& peanoclaw::runners::PeanoClawLibraryRunner::getState() {
    return _repository->getState();
}

void peanoclaw::runners::PeanoClawLibraryRunner::configureGlobalTimestep(double time) {
  _repository->getState().setGlobalTimestepEndTime(time);
  _repository->getState().setNumerics(_numerics);
  _repository->getState().setPlotNumber(_plotNumber);
}

void peanoclaw::runners::PeanoClawLibraryRunner::runNextPossibleTimestep() {
    bool plotSubsteps = _configuration.plotSubsteps()
      || (_configuration.plotSubstepsAfterOutputTime() != -1 && _configuration.plotSubstepsAfterOutputTime() <= _plotNumber);

    logInfo("evolveToTime", "Solving timestep " << (_plotNumber-1) << " with maximum global time interval ("
        << _repository->getState().getStartMaximumGlobalTimeInterval() << ", " << _repository->getState().getEndMaximumGlobalTimeInterval() << ")"
        << " and minimum global time interval (" << _repository->getState().getStartMinimumGlobalTimeInterval() << ", " << _repository->getState().getEndMinimumGlobalTimeInterval() << ")");
    _iterationTimer.startTimer();

    _repository->getState().resetGlobalTimeIntervals();
    _repository->getState().resetMinimalTimestep();
    _repository->getState().setAllPatchesEvolvedToGlobalTimestep(true);


    if(plotSubsteps) {
      _repository->getState().setPlotNumber(_plotNumber);
      if(_validateGrid) {
        _repository->switchToSolveTimestepAndPlotAndValidateGrid();
      } else {
        _repository->switchToSolveTimestepAndPlot();
      }
     
      std::cout << "before solve timestep: " << _repository->getState() << std::endl;
      _repository->iterate();
      std::cout << "after solve timestep: " << _repository->getState() << std::endl;

      /*do {
        _repository->iterate();
      } while (!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced());*/

      _plotNumber++;
    } else {
      if(_validateGrid) {
        _repository->switchToSolveTimestepAndValidateGrid();
      } else {
        _repository->switchToSolveTimestep();
      }
 
      std::cout << "before solve timestep: " << _repository->getState() << std::endl;
      _repository->iterate();
      std::cout << "after solve timestep: " << _repository->getState() << std::endl;

      /*_repository->switchToRemesh();
      do {
        _repository->iterate();
      } while (!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced());*/

    }


#if 0
      bool hasChangedState = false;
      do {
            hasChangedState = false;
            _repository->switchToRemesh();
            _repository->iterate();
            std::cout << "remesh iteration block (1)" << _repository->getState() << std::endl;
            hasChangedState |= !_repository->getState().isGridStationary();

            /*_repository->iterate(); 
            std::cout << "remesh iteration block (2)" << _repository->getState() << std::endl;
            hasChangedState |= !_repository->getState().isGridStationary();

            _repository->iterate(); 
            std::cout << "remesh iteration block (3)" << _repository->getState() << std::endl;
            hasChangedState |= !_repository->getState().isGridStationary();

            _repository->iterate();
            std::cout << "remesh iteration block (4)" << _repository->getState() << std::endl;
            hasChangedState |= !_repository->getState().isGridStationary();*/
      } while(!_repository->getState().isGridStationary() || !_repository->getState().isGridBalanced() || hasChangedState);
#endif


    _repository->getState().plotStatisticsForLastGridIteration();

    _iterationTimer.stopTimer();
    _totalRuntime += (double)_iterationTimer.getCPUTicks() / (double)CLOCKS_PER_SEC;
    logInfo("evolveToTime", "Wallclock time for this grid iteration/Total runtime: " << _iterationTimer.getCalendarTime() << "s/" << _totalRuntime << "s");
    logInfo("evolveToTime", "Minimal timestep for this grid iteration: " << _repository->getState().getMinimalTimestep());

    assertion(_repository->getState().getMinimalTimestep() < std::numeric_limits<double>::infinity());
}
