#include "peano/heap/Heap.h"

#include "mappings/SolveTimestep.h"
#include "Patch.h"


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::SolveTimestep::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::SolveTimestep::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::SolveTimestep::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::SolveTimestep::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::SolveTimestep::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::SolveTimestep::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


tarch::logging::Log                peanoclaw::mappings::SolveTimestep::_log( "peanoclaw::mappings::SolveTimestep" ); 


peanoclaw::mappings::SolveTimestep::SolveTimestep() {
  logTraceIn( "SolveTimestep()" );
  // @todo Insert your code here
  logTraceOut( "SolveTimestep()" );
}


peanoclaw::mappings::SolveTimestep::~SolveTimestep() {
  logTraceIn( "~SolveTimestep()" );
  // @todo Insert your code here
  logTraceOut( "~SolveTimestep()" );
}


#if defined(SharedMemoryParallelisation)
peanoclaw::mappings::SolveTimestep::SolveTimestep(const SolveTimestep&  masterThread) {
  logTraceIn( "SolveTimestep(SolveTimestep)" );
  // @todo Insert your code here
  logTraceOut( "SolveTimestep(SolveTimestep)" );
}


void peanoclaw::mappings::SolveTimestep::mergeWithWorkerThread(const SolveTimestep& workerThread) {
  logTraceIn( "mergeWithWorkerThread(SolveTimestep)" );
  // @todo Insert your code here
  logTraceOut( "mergeWithWorkerThread(SolveTimestep)" );
}
#endif


void peanoclaw::mappings::SolveTimestep::createHangingVertex(
      peanoclaw::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      peanoclaw::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      peanoclaw::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createHangingVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::destroyHangingVertex(
      const peanoclaw::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "destroyHangingVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::createInnerVertex(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createInnerVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createInnerVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::createBoundaryVertex(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createBoundaryVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createBoundaryVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::destroyVertex(
      const peanoclaw::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "destroyVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::createCell(
      peanoclaw::Cell&                 fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "createCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
 
  if(static_cast<int>(_levelStatistics.size()) < fineGridVerticesEnumerator.getLevel()) {
    _levelStatistics.push_back(peanoclaw::statistics::LevelInformation());
  }
  peanoclaw::statistics::LevelInformation& levelInformation = _levelStatistics.at(fineGridVerticesEnumerator.getLevel()-1);
  levelInformation._createdPatches++;

  logTraceOutWith1Argument( "createCell(...)", fineGridCell );
}


void peanoclaw::mappings::SolveTimestep::destroyCell(
      const peanoclaw::Cell&           fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "destroyCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
 
  if(static_cast<int>(_levelStatistics.size()) < fineGridVerticesEnumerator.getLevel()) {
    _levelStatistics.push_back(peanoclaw::statistics::LevelInformation());
  }
  peanoclaw::statistics::LevelInformation& levelInformation = _levelStatistics.at(fineGridVerticesEnumerator.getLevel()-1);
  levelInformation._destroyedPatches++;

  logTraceOutWith1Argument( "destroyCell(...)", fineGridCell );
}

#ifdef Parallel
void peanoclaw::mappings::SolveTimestep::mergeWithNeighbour(
  peanoclaw::Vertex&  vertex,
  const peanoclaw::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );
  // @todo Insert your code here
  logTraceOut( "mergeWithNeighbour(...)" );
}

void peanoclaw::mappings::SolveTimestep::prepareSendToNeighbour(
  peanoclaw::Vertex&  vertex,
  int  toRank,
  int  level
) {
  logTraceInWith3Arguments( "prepareSendToNeighbour(...)", vertex, toRank, level );
  // @todo Insert your code here
  logTraceOut( "prepareSendToNeighbour(...)" );
}

void peanoclaw::mappings::SolveTimestep::prepareCopyToRemoteNode(
  peanoclaw::Vertex&  localVertex,
  int  toRank
) {
  logTraceInWith2Arguments( "prepareCopyToRemoteNode(...)", localVertex, toRank );
  // @todo Insert your code here
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void peanoclaw::mappings::SolveTimestep::prepareCopyToRemoteNode(
  peanoclaw::Cell&  localCell,
  int  toRank
) {
  logTraceInWith2Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank );
  // @todo Insert your code here
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void peanoclaw::mappings::SolveTimestep::mergeWithRemoteDataDueToForkOrJoin(
  peanoclaw::Vertex&  localVertex,
  const peanoclaw::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
  logTraceInWith6Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localVertex, masterOrWorkerVertex, fromRank, x, h, level );
  // @todo Insert your code here
  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void peanoclaw::mappings::SolveTimestep::mergeWithRemoteDataDueToForkOrJoin(
  peanoclaw::Cell&  localCell,
  const peanoclaw::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
  logTraceInWith3Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localCell, masterOrWorkerCell, fromRank );
  // @todo Insert your code here
  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void peanoclaw::mappings::SolveTimestep::prepareSendToWorker(
  peanoclaw::Cell&                 fineGridCell,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  logTraceIn( "prepareSendToWorker(...)" );
  // @todo Insert your code here
  logTraceOut( "prepareSendToWorker(...)" );
}

void peanoclaw::mappings::SolveTimestep::prepareSendToMaster(
  peanoclaw::Cell&     localCell,
  peanoclaw::Vertex *  vertices,
  const peano::grid::VertexEnumerator&  verticesEnumerator
) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "prepareSendToMaster(...)" );
}


void peanoclaw::mappings::SolveTimestep::mergeWithMaster(
  const peanoclaw::Cell&           workerGridCell,
  peanoclaw::Vertex * const        workerGridVertices,
 const peano::grid::VertexEnumerator& workerEnumerator,
  peanoclaw::Cell&                 fineGridCell,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
  const peanoclaw::State&          workerState,
  peanoclaw::State&                masterState
) {
  logTraceIn( "mergeWithMaster(...)" );
  // @todo Insert your code here
  logTraceOut( "mergeWithMaster(...)" );
}


void peanoclaw::mappings::SolveTimestep::receiveDataFromMaster(
  peanoclaw::Cell&                    receivedCell, 
  peanoclaw::Vertex *                 receivedVertices,
  const peano::grid::VertexEnumerator&    verticesEnumerator
) {
  logTraceInWith2Arguments( "receiveDataFromMaster(...)", receivedCell.toString(), verticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "receiveDataFromMaster(...)" );
}


void peanoclaw::mappings::SolveTimestep::mergeWithWorker(
  peanoclaw::Cell&           localCell, 
  const peanoclaw::Cell&     receivedMasterCell
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );
  // @todo Insert your code here
  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}


void peanoclaw::mappings::SolveTimestep::mergeWithWorker(
  peanoclaw::Vertex&        localVertex,
  const peanoclaw::Vertex&  receivedMasterVertex
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localVertex.toString(), receivedMasterVertex.toString() );
  // @todo Insert your code here
  logTraceOutWith1Argument( "mergeWithWorker(...)", localVertex.toString() );
}
#endif

void peanoclaw::mappings::SolveTimestep::touchVertexFirstTime(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexFirstTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "touchVertexFirstTime(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::touchVertexLastTime(
      peanoclaw::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexLastTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
 
  // Application driven refinement control
  if(
      fineGridVertex.shouldRefine()
      //&& fineGridVertex.getRefinementControl() != peanoclaw::Vertex::Records::Refined // MARK roland
      && fineGridVertex.getRefinementControl() == peanoclaw::Vertex::Records::Unrefined // MARK roland
    ) {
    fineGridVertex.refine();
  } else if (
      fineGridVertex.shouldErase()
      && fineGridVertex.getCurrentAdjacentCellsHeight() == 1
    ) {
    fineGridVertex.erase();
  }
  fineGridVertex.setShouldRefine(false);
  fineGridVertex.resetSubcellsEraseVeto();

  //Predefined refinement
//    double radius = 0.2 + _averageGlobalTimeInterval / 0.2 * 0.3;
    double radius = 0.15 + _averageGlobalTimeInterval / 0.5 * 0.3;

    bool shouldBeRefined
     = (tarch::la::oneGreater(fineGridH, _initialMinimalMeshWidth)) ||
      ((std::abs(tarch::la::norm2(fineGridX-0.5)-radius) < fineGridH(0) / 2.0)
          && (tarch::la::oneGreater(fineGridH, _initialMinimalMeshWidth / (std::pow(3.0, (double)_additionalLevelsForPredefinedRefinement)))));

    if(shouldBeRefined) {
      if(fineGridVertex.getRefinementControl() == peanoclaw::records::Vertex::Unrefined) {
        //TODO unterweg debug
//        std::cout << "Refining vertex at " << fineGridX << " on level " << (coarseGridVerticesEnumerator.getLevel() + 1) << std::endl;
//        fineGridVertex.refine();
      }
    } else {
      if(fineGridVertex.getRefinementControl() == peanoclaw::records::Vertex::Refined
          && fineGridVertex.getAdjacentCellsHeightOfPreviousIteration() == 1
          && fineGridVertex.getCurrentAdjacentCellsHeight() == 1
        ) {
        //TODO unterweg debug
//        std::cout << "Coarsening vertex at " << fineGridX << " on level " << (coarseGridVerticesEnumerator.getLevel() + 1) << std::endl;
//        fineGridVertex.erase();
      }
    }

  logTraceOutWith1Argument( "touchVertexLastTime(...)", fineGridVertex );
}


void peanoclaw::mappings::SolveTimestep::enterCell(
      peanoclaw::Cell&                 fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
 
  //Create patch
  Patch patch = Patch(
    fineGridVerticesEnumerator.getVertexPosition(0),
    fineGridVerticesEnumerator.getCellSize(),
    fineGridCell
  );

  //Histogram
  assertion(static_cast<int>(_levelStatistics.size()) >= coarseGridVerticesEnumerator.getLevel());
  if(static_cast<int>(_levelStatistics.size()) < fineGridVerticesEnumerator.getLevel()) {
    _levelStatistics.push_back(peanoclaw::statistics::LevelInformation());
  }
  peanoclaw::statistics::LevelInformation& levelInformation = _levelStatistics.at(fineGridVerticesEnumerator.getLevel()-1);
  if(patch.isLeaf()) {
    levelInformation._numberOfPatches++;
    levelInformation._numberOfCells += (patch.getSubdivisionFactor() * patch.getSubdivisionFactor());
    levelInformation._area += tarch::la::volume(patch.getSize());
  }

  //Solve timestep for this patch
  if(fineGridCell.isLeaf()) {
    #ifdef Asserts
    CellDescription& cellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(fineGridCell.getCellDescriptionIndex()).at(0);
    double startTime = cellDescription.getTime();
    double endTime = cellDescription.getTime() + cellDescription.getTimestepSize();
    assertionEquals2(patch.getCurrentTime(), startTime, startTime, patch.toString());
    assertionEquals2(patch.getCurrentTime() + patch.getTimestepSize(), endTime, endTime, patch.toString());
    #endif

    assertion(patch.isLeaf() || patch.isVirtual());

    //Perform timestep
    double maximumTimestepSize = _globalTimestepEndTime - (patch.getCurrentTime() + patch.getTimestepSize());
    if (
         tarch::la::greater(maximumTimestepSize, 0.0)
         && patch.isAllowedToAdvanceInTime()
         && !patch.shouldSkipNextGridIteration()
       ) {

      // Copy uNew to uOld
//      patch.copyUNewToUOld();

      // Set boundary conditions (For hanging nodes the isBoundary() is not valid. Therefore, we simulate it by checking for the domainoffset and -size.)
      if((fineGridVertices[fineGridVerticesEnumerator(0)].isBoundary()
              || !tarch::la::allGreater(fineGridVerticesEnumerator.getVertexPosition(0), _domainOffset) || !tarch::la::allGreater(_domainOffset+_domainSize, fineGridVerticesEnumerator.getVertexPosition(0)))
              && fineGridVertices[fineGridVerticesEnumerator(0)].getAdjacentCellDescriptionIndex(1) == -1) {
        _pyClaw->fillLeftBoundaryLayer(patch);
      }
      if((fineGridVertices[fineGridVerticesEnumerator(2)].isBoundary()
          || !tarch::la::allGreater(fineGridVerticesEnumerator.getVertexPosition(2), _domainOffset) || !tarch::la::allGreater(_domainOffset+_domainSize, fineGridVerticesEnumerator.getVertexPosition(2)))
          && fineGridVertices[fineGridVerticesEnumerator(2)].getAdjacentCellDescriptionIndex(0) == -1) {
        _pyClaw->fillUpperBoundaryLayer(patch);
      }
      if((fineGridVertices[fineGridVerticesEnumerator(1)].isBoundary()
          || !tarch::la::allGreater(fineGridVerticesEnumerator.getVertexPosition(1), _domainOffset) || !tarch::la::allGreater(_domainOffset+_domainSize, fineGridVerticesEnumerator.getVertexPosition(1)))
          && fineGridVertices[fineGridVerticesEnumerator(1)].getAdjacentCellDescriptionIndex(0) == -1) {
        _pyClaw->fillRightBoundaryLayer(patch);
      }
      if((fineGridVertices[fineGridVerticesEnumerator(0)].isBoundary()
          || !tarch::la::allGreater(fineGridVerticesEnumerator.getVertexPosition(0), _domainOffset) || !tarch::la::allGreater(_domainOffset+_domainSize, fineGridVerticesEnumerator.getVertexPosition(0)))
          && fineGridVertices[fineGridVerticesEnumerator(0)].getAdjacentCellDescriptionIndex(2) == -1) {
        _pyClaw->fillLowerBoundaryLayer(patch);
      }

      //Do one timestep...
      double requiredMeshWidth = _pyClaw->solveTimestep(patch, maximumTimestepSize, _useDimensionalSplitting);
      patch.setDemandedMeshWidth(requiredMeshWidth);

      //Coarse grid correction
      for(int i = 0; i < TWO_POWER_D; i++) {
        if(fineGridVertices[fineGridVerticesEnumerator(i)].isHangingNode()) {
          fineGridVertices[fineGridVerticesEnumerator(i)].applyCoarseGridCorrection(*_pyClaw);
        }
      }

      //Statistics
      if(tarch::la::greater(patch.getTimestepSize(), 0.0)) {
        levelInformation._numberOfCellUpdates += (patch.getSubdivisionFactor() * patch.getSubdivisionFactor());
        _minimalTimestep = std::min(_minimalTimestep, patch.getTimestepSize());
      }

      //Probes
      for(std::vector<peanoclaw::statistics::Probe>::iterator i = _probeList.begin();
          i != _probeList.end();
          i++) {
        i->plotDataIfContainedInPatch(patch);
      }

      logDebug("enterCell(...)", "New time interval of patch " << fineGridVerticesEnumerator.getCellCenter() << " on level " << fineGridVerticesEnumerator.getLevel() << " is [" << patch.getCurrentTime() << ", " << (patch.getCurrentTime() + patch.getTimestepSize()) << "]");
    } else {
      //TODO unterweg debug
//      std::cout << "Not advancing in time due to _globalTimestepEndTime=" << _globalTimestepEndTime
//          << " currenttime=" << patch.getCurrentTime() << " timestepsize=" << patch.getTimestepSize()
//          << " maximumTimestep=" << maximumTimestepSize << " isAllowedToAdvance=" << patch.isAllowedToAdvanceInTime()
//          << " minimalNeighborTime=" << patch.getMinimalNeighborTime() << std::endl;
      logDebug("enterCell(...)", "Unchanged time interval of patch " << fineGridVerticesEnumerator.getCellCenter() << " on level " << fineGridVerticesEnumerator.getLevel() << " is [" << patch.getCurrentTime() << ", " << (patch.getCurrentTime() + patch.getTimestepSize()) << "]");
      patch.setSkipNextGridIteration(false);
    }

    //Refinement criterion
    if(tarch::la::oneGreater(patch.getSubcellSize(), tarch::la::Vector<DIMENSIONS, double>(patch.getDemandedMeshWidth()))) {
      // Refine
      for(int i = 0; i < TWO_POWER_D; i++) {
        if(!fineGridVertices[fineGridVerticesEnumerator(i)].isHangingNode()) {
          fineGridVertices[fineGridVerticesEnumerator(i)].setShouldRefine(true);
          coarseGridVertices[coarseGridVerticesEnumerator(i)].setSubcellEraseVeto(i);
        }
      }
    } else if (!tarch::la::oneGreater(patch.getSubcellSize() * 3.0, tarch::la::Vector<DIMENSIONS, double>(patch.getDemandedMeshWidth()))) {
      // Coarsen
    } else {
      for(int i = 0; i < TWO_POWER_D; i++) {
        coarseGridVertices[coarseGridVerticesEnumerator(i)].setSubcellEraseVeto(i);
        if(fineGridVertices[fineGridVerticesEnumerator(i)].isHangingNode()
            && !coarseGridVertices[coarseGridVerticesEnumerator(i)].isHangingNode()) {
          coarseGridVertices[coarseGridVerticesEnumerator(i)].setShouldRefine(true);
        }
      }
    }

    //Stopping criterion for global timestep
    if(tarch::la::smaller(patch.getCurrentTime() + patch.getTimestepSize(), _globalTimestepEndTime)) {
      _allPatchesEvolvedToGlobalTimestep = false;
    }

    assertion2(!tarch::la::smaller(patch.getCurrentTime(), startTime), patch, startTime);
    assertion2(!tarch::la::smaller(patch.getCurrentTime() + patch.getTimestepSize(), endTime), patch.getCurrentTime() + patch.getTimestepSize(), endTime);

    #ifdef Asserts

    // roland MARK
    if(patch.containsNaN()) {
      logError("", "Invalid solution in patch " << patch.toString()
          << std::endl << patch.toStringUNew() << std::endl << patch.toStringUOldWithGhostLayer());
      //TODO unterweg
      Patch coarsePatch(peano::heap::Heap<CellDescription>::getInstance().getData(coarseGridCell.getCellDescriptionIndex()).at(0));
      std::cout << "Coarse Patch:" << std::endl << coarsePatch.toString() << std::endl << coarsePatch.toStringUNew() << std::endl;

      for(int i = 0; i < TWO_POWER_D; i++) {
        std::cout << "Vertex" << i <<": " << fineGridVertices[fineGridVerticesEnumerator(i)].toString() << std::endl;
        for(int j = 0; j < TWO_POWER_D; j++) {
          if(fineGridVertices[fineGridVerticesEnumerator(i)].getAdjacentCellDescriptionIndex(j) != -1) {
            CellDescription& cellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(fineGridVertices[fineGridVerticesEnumerator(i)].getAdjacentCellDescriptionIndex(j)).at(0);
            Patch neighborPatch(cellDescription);
            std::cout << i << " " << j << std::endl
                << neighborPatch.toString() << std::endl
                << neighborPatch.toStringUNew() << std::endl
                << neighborPatch.toStringUOldWithGhostLayer() << std::endl;
          } else {
            std::cout << i << ": Invalid patch";
          }
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      assertion(false);
      throw "";
    }
    #endif

    logTraceOutWith2Arguments( "enterCell(...)", cellDescription.getTimestepSize(), cellDescription.getTime() + cellDescription.getTimestepSize() );
  } else {
    logTraceOut( "enterCell(...)" );
  }

  //Statistics
  _startMaximumLocalTimeInterval = std::min(patch.getCurrentTime(), _startMaximumLocalTimeInterval);
  _endMaximumLocalTimeInterval = std::max(patch.getCurrentTime() + patch.getTimestepSize(), _endMaximumLocalTimeInterval);
  _startMinimumLocalTimeInterval = std::max(patch.getCurrentTime(), _startMinimumLocalTimeInterval);
  _endMinimumLocalTimeInterval = std::min(patch.getCurrentTime() + patch.getTimestepSize(), _endMinimumLocalTimeInterval);

  patch.increaseAgeByOneGridIteration();
}


void peanoclaw::mappings::SolveTimestep::leaveCell(
      peanoclaw::Cell&           fineGridCell,
      peanoclaw::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "leaveCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // @todo Insert your code here
  logTraceOutWith1Argument( "leaveCell(...)", fineGridCell );
}


void peanoclaw::mappings::SolveTimestep::beginIteration(
  peanoclaw::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );
 
  _pyClaw = &solverState.getPyClaw();
  _globalTimestepEndTime = solverState.getGlobalTimestepEndTime();
  _allPatchesEvolvedToGlobalTimestep = solverState.getAllPatchesEvolvedToGlobalTimestep();
  _startMaximumLocalTimeInterval = std::numeric_limits<double>::max();
  _endMaximumLocalTimeInterval = -std::numeric_limits<double>::max();
  _startMinimumLocalTimeInterval = -std::numeric_limits<double>::max();
  _endMinimumLocalTimeInterval = std::numeric_limits<double>::max();
//  _numberOfCellUpdates = 0;
  _minimalTimestep = std::numeric_limits<double>::max();
  _domainOffset = solverState.getDomainOffset();
  _domainSize = solverState.getDomainSize();
  _initialMinimalMeshWidth = solverState.getInitialMinimalMeshWidth();
  _additionalLevelsForPredefinedRefinement = solverState.getAdditionalLevelsForPredefinedRefinement();
  _probeList = solverState.getProbeList();
  _useDimensionalSplitting = solverState.useDimensionalSplitting();

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void peanoclaw::mappings::SolveTimestep::endIteration(
  peanoclaw::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
 
  solverState.setAllPatchesEvolvedToGlobalTimestep(
    solverState.getAllPatchesEvolvedToGlobalTimestep()
    && _allPatchesEvolvedToGlobalTimestep
  );
  solverState.updateGlobalTimeIntervals(
    _startMaximumLocalTimeInterval,
    _endMaximumLocalTimeInterval,
    _startMinimumLocalTimeInterval,
    _endMinimumLocalTimeInterval
  );
//  solverState.addToTotalNumberOfCellUpdates(_numberOfCellUpdates);
  solverState.updateMinimalTimestep(_minimalTimestep);

  static int counter = 0;
  counter++;

  solverState.setLevelStatisticsForLastGridIteration(_levelStatistics);
  _levelStatistics.clear();

  _averageGlobalTimeInterval = (solverState.getStartMaximumGlobalTimeInterval() + solverState.getEndMaximumGlobalTimeInterval()) / 2.0;

  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void peanoclaw::mappings::SolveTimestep::descend(
  peanoclaw::Cell * const          fineGridCells,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell
) {
  logTraceInWith2Arguments( "descend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "descend(...)" );
}


void peanoclaw::mappings::SolveTimestep::ascend(
  peanoclaw::Cell * const    fineGridCells,
  peanoclaw::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  peanoclaw::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  peanoclaw::Cell&           coarseGridCell
) {
  logTraceInWith2Arguments( "ascend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "ascend(...)" );
}