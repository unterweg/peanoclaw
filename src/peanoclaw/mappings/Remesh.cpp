#include "peanoclaw/mappings/Remesh.h"

#include "peanoclaw/interSubgridCommunication/GridLevelTransfer.h"
#include "peanoclaw/Numerics.h"
#include "peanoclaw/parallel/NeighbourCommunicator.h"
#include "peanoclaw/parallel/MasterWorkerAndForkJoinCommunicator.h"

#include "peano/heap/Heap.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

peanoclaw::records::VertexDescription::IterationParity peanoclaw::mappings::Remesh::_iterationParity
 = peanoclaw::records::VertexDescription::EVEN;

std::map<tarch::la::Vector<DIMENSIONS_PLUS_ONE,double> , peanoclaw::mappings::Remesh::VertexDescription, tarch::la::VectorCompare<DIMENSIONS_PLUS_ONE> >   peanoclaw::mappings::Remesh::_vertexPositionToIndexMap;

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


tarch::logging::Log                peanoclaw::mappings::Remesh::_log( "peanoclaw::mappings::Remesh" ); 


void peanoclaw::mappings::Remesh::dataExchange_receive(
  peanoclaw::Cell&  localCell,
  const peanoclaw::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                       level,
  bool                                      forkOrJoin
) {
      std::cout << "receiving from: " << fromRank << " for " << cellCentre << " and " << cellSize << " on level " << level << std::endl;

      peano::heap::MessageType messageType = forkOrJoin ? peano::heap::ForkOrJoinCommunication : peano::heap::MasterWorkerCommunication;

      std::vector<peanoclaw::records::CellDescription> remoteCellDescriptionVector = peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().receiveData(fromRank, cellCentre, level, messageType);
      CellDescription& remoteCellDescription = remoteCellDescriptionVector.at(0);
 
      if (!peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().isValidIndex(localCell.getCellDescriptionIndex())) {
        // clone existing remote cell Description
        int localCellDescriptionIndex = peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().createData();
        std::vector<CellDescription>& localCellDescriptionVector = peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().getData(localCellDescriptionIndex);

        // copy modified remote CellDescription to local slot
        localCellDescriptionVector.push_back(remoteCellDescription);
 
        localCell.setCellDescriptionIndex(localCellDescriptionIndex);

        CellDescription& localCellDescription = localCellDescriptionVector.at(0);

        // erase array indices from remote Cell to force array allocation
        localCellDescription.setUNewIndex(-1);
        localCellDescription.setUOldIndex(-1);
        localCellDescription.setAuxIndex(-1);
        std::cout << "cloned cell description" << std::endl;
      }

      std::vector<peanoclaw::records::CellDescription>& localCellDescriptionVector = peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex());
      CellDescription& localCellDescription = localCellDescriptionVector.at(0);
      
      // receive data parts if indices != -1
      if (remoteCellDescription.getUNewIndex() != -1) {
        assertion(remoteCellDescription.getUNewIndex() > 0);
        std::vector<Data> remoteMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().receiveData(fromRank, cellCentre, level, messageType);
        
        //create local Vector if needed
        int localIndex = localCellDescription.getUNewIndex();
        if (!peano::heap::Heap<peanoclaw::records::Data>::getInstance().isValidIndex(localIndex)) {
            localIndex = peano::heap::Heap<peanoclaw::records::Data>::getInstance().createData();
            localCellDescription.setUNewIndex(localIndex);
            std::cout << "created uNew vector" << std::endl;
        }

        // get local vector reference
        std::vector<Data>& localMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().getData(localIndex);

        // copy contents:
        std::vector<Data>::iterator it = remoteMatrix.begin();
        localMatrix.assign(it, remoteMatrix.end());

        std::cout << "uNew: copied " << remoteMatrix.size() << " elements" << std::endl;
      }

      if (remoteCellDescription.getUOldIndex() != -1) {
        std::vector<Data> remoteMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().receiveData(fromRank, cellCentre, level, messageType);
        
        //create local Vector if needed
        int localIndex = localCellDescription.getUOldIndex();
        if (!peano::heap::Heap<peanoclaw::records::Data>::getInstance().isValidIndex(localIndex)) {
            localIndex = peano::heap::Heap<peanoclaw::records::Data>::getInstance().createData();
            localCellDescription.setUOldIndex(localIndex);
            std::cout << "created uOld vector" << std::endl;
        }

        // get local vector reference
        std::vector<Data>& localMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().getData(localIndex);

        // copy contents: (TODO: this is wrong most of the time! use assign)
        std::vector<Data>::iterator it = remoteMatrix.begin();
        localMatrix.assign(it, remoteMatrix.end());

        std::cout << "uOld: copied " << remoteMatrix.size() << " elements" << std::endl;
      }

      if (remoteCellDescription.getAuxIndex() != -1) {
        std::vector<Data> remoteMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().receiveData(fromRank, cellCentre, level, messageType);
        
        //create local Vector if needed
        int localIndex = localCellDescription.getAuxIndex();
        if (!peano::heap::Heap<peanoclaw::records::Data>::getInstance().isValidIndex(localIndex)) {
            localIndex = peano::heap::Heap<peanoclaw::records::Data>::getInstance().createData();
            localCellDescription.setAuxIndex(localIndex);
            std::cout << "created Aux vector" << std::endl;
        }

        // get local vector reference
        std::vector<Data>& localMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().getData(localIndex);

        // copy contents: (TODO: this is wrong most of the time! use assign)
        std::vector<Data>::iterator it = remoteMatrix.begin();
        localMatrix.assign(it, remoteMatrix.end());
        
        std::cout << "Aux: copied " << remoteMatrix.size() << " elements" << std::endl;
      }
}

void peanoclaw::mappings::Remesh::dataExchange_send(
  peanoclaw::Cell&  localCell,
  int  toRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level,
  bool                                         forkOrJoin
) {
  peano::heap::MessageType messageType = forkOrJoin ? peano::heap::ForkOrJoinCommunication : peano::heap::MasterWorkerCommunication;

  // get local Cell Description
  assertion3(localCell.getCellDescriptionIndex() > 0, cellCentre, level, toRank);
  CellDescription& cellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0);

  // send Cell description
  std::cout << "sending index to " << toRank << " and position=" << cellCentre << " on level " << level << std::endl;
  peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().sendData(localCell.getCellDescriptionIndex(), toRank, cellCentre, level, messageType);

  Patch localPatch(cellDescription);

  if (cellDescription.getUNewIndex() != -1) {
    int localIndex = cellDescription.getUNewIndex();
    std::vector<Data>& localMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().getData(localIndex);
    std::cout << "sending new data with elements: " << localMatrix.size() << std::endl;
    peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(localIndex, toRank, cellCentre, level, messageType);
  }

  if (cellDescription.getUOldIndex() != -1) {
    int localIndex = cellDescription.getUOldIndex();
    std::vector<Data>& localMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().getData(localIndex);
    std::cout << "sending old data with elements: " << localMatrix.size() << std::endl;
    peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(localIndex, toRank, cellCentre, level, messageType);
  }

  if (cellDescription.getAuxIndex() != -1) {
    int localIndex = cellDescription.getAuxIndex();
    std::vector<Data>& localMatrix = peano::heap::Heap<peanoclaw::records::Data>::getInstance().getData(localIndex);
    std::cout << "sending aux data with elements: " << localMatrix.size() << std::endl;
    peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(localIndex, toRank, cellCentre, level, messageType);
  }
}


peanoclaw::mappings::Remesh::Remesh() {
  logTraceIn( "Remesh()" );
  // @todo Insert your code here
  logTraceOut( "Remesh()" );
}


peanoclaw::mappings::Remesh::~Remesh() {
  logTraceIn( "~Remesh()" );
  // @todo Insert your code here
  logTraceOut( "~Remesh()" );
}


#if defined(SharedMemoryParallelisation)
peanoclaw::mappings::Remesh::Remesh(const Remesh&  masterThread) {
  logTraceIn( "Remesh(Remesh)" );
  // @todo Insert your code here
  logTraceOut( "Remesh(Remesh)" );
}


void peanoclaw::mappings::Remesh::mergeWithWorkerThread(const Remesh& workerThread) {
  logTraceIn( "mergeWithWorkerThread(Remesh)" );
  // @todo Insert your code here
  logTraceOut( "mergeWithWorkerThread(Remesh)" );
}
#endif


void peanoclaw::mappings::Remesh::createHangingVertex(
      peanoclaw::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      peanoclaw::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      peanoclaw::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
 
  //TODO unterweg debug
//  std::cout << "Create hanging vertex" << std::endl;

  fineGridVertex.setShouldRefine(false);
  fineGridVertex.resetSubcellsEraseVeto();

  if(!tarch::la::oneGreater(_domainOffset, fineGridX) && !tarch::la::oneGreater(fineGridX, _domainOffset + _domainSize)) {
    //Project adjacency information down from coarse grid vertex
    _gridLevelTransfer->fillAdjacentPatchIndicesFromCoarseVertices(
      coarseGridVertices,
      coarseGridVerticesEnumerator,
      fineGridVertex,
      fineGridPositionOfVertex
    );

    //Retrieve or create hanging vertex description
    tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> hangingVertexPosition;
    for(int d = 0; d < DIMENSIONS; d++) {
      hangingVertexPosition(d) = fineGridX(d);
    }
    hangingVertexPosition(DIMENSIONS) = coarseGridVerticesEnumerator.getLevel() + 1;

    if( _vertexPositionToIndexMap.find(hangingVertexPosition) == _vertexPositionToIndexMap.end() ) {
      VertexDescription vertexDescription;
      vertexDescription.setTouched(true);
      for(int i = 0; i < TWO_POWER_D; i++) {
        vertexDescription.setIndicesOfAdjacentCellDescriptions(i, -1);
      }
      if(_iterationParity == peanoclaw::records::VertexDescription::EVEN) {
        vertexDescription.setLastUpdateIterationParity(peanoclaw::records::VertexDescription::ODD);
      } else {
        vertexDescription.setLastUpdateIterationParity(peanoclaw::records::VertexDescription::EVEN);
      }
      _vertexPositionToIndexMap[hangingVertexPosition] = vertexDescription;
    }

    VertexDescription& hangingVertexDescription = _vertexPositionToIndexMap[hangingVertexPosition];
    hangingVertexDescription.setTouched(true);

    //Copy indices from coarse level
    _gridLevelTransfer->fillAdjacentPatchIndicesFromCoarseVertices(
      coarseGridVertices,
      coarseGridVerticesEnumerator,
      fineGridVertex,
      fineGridPositionOfVertex
    );

    //Remove deleted indices
    for(int i = 0; i < TWO_POWER_D; i++) {
      //From hanging vertex description
      if(hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i) != -1
          && !peano::heap::Heap<CellDescription>::getInstance().isValidIndex(hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i))) {
        hangingVertexDescription.setIndicesOfAdjacentCellDescriptions(i, -1);
      }
      if(fineGridVertex.getAdjacentCellDescriptionIndex(i) != -1
          && !peano::heap::Heap<CellDescription>::getInstance().isValidIndex(fineGridVertex.getAdjacentCellDescriptionIndex(i))) {
        fineGridVertex.setAdjacentCellDescriptionIndex(i, -1);
      }
    }

    //TODO If the coarse grid vertices are also hanging, a deleted patch index two or more
    //levels coarser than this hanging vertex might not be recognized, yet.
    //Merging adjacency information from stored hanging vertex description and hanging vertex
    //The data stored on the hanging vertex itself must come from the coarser vertex, since
    //the hanging vertex has just been created. So, this data is more recent, when the data
    //in the hanging vertex description describes a patch on a coarser level than the hanging
    //vertex (Should be solved by the check before).
    for(int i = 0; i < TWO_POWER_D; i++) {
      assertion(hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i) >= -1);
      if(hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i) != -1) {
        CellDescription& cellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i)).at(0);
        if(cellDescription.getLevel() == (coarseGridVerticesEnumerator.getLevel() + 1)) {
          fineGridVertex.setAdjacentCellDescriptionIndex(i, hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i));
        }
      }
    }
  }

  logTraceOutWith1Argument( "createHangingVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::Remesh::destroyHangingVertex(
      const peanoclaw::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
 
  //Handle refinement flags
  if(tarch::la::allGreater(fineGridX, _domainOffset) && tarch::la::allGreater(_domainOffset + _domainSize, fineGridX)) {
    _gridLevelTransfer->restrictRefinementFlagsToCoarseVertices(
      coarseGridVertices,
      coarseGridVerticesEnumerator,
      fineGridVertex,
      fineGridPositionOfVertex
    );
  }

  if(!tarch::la::oneGreater(_domainOffset, fineGridX) && !tarch::la::oneGreater(fineGridX, _domainOffset + _domainSize)) {

    //Retrieve hanging vertex description
    tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> hangingVertexPosition;
    for(int d = 0; d < DIMENSIONS; d++) {
      hangingVertexPosition(d) = fineGridX(d);
    }
    hangingVertexPosition(DIMENSIONS) = coarseGridVerticesEnumerator.getLevel() + 1;

    assertionMsg(_vertexPositionToIndexMap.find(hangingVertexPosition) != _vertexPositionToIndexMap.end(), "Hanging vertex description was not created for vertex " << fineGridVertex);

    VertexDescription& hangingVertexDescription = _vertexPositionToIndexMap[hangingVertexPosition];
    hangingVertexDescription.setTouched(true);

    //Copy adjacency information from hanging vertex to hanging vertex description
    for(int i = 0; i < TWO_POWER_D; i++) {
      hangingVertexDescription.setIndicesOfAdjacentCellDescriptions(i, fineGridVertex.getAdjacentCellDescriptionIndex(i));
    }

    //Fill boundary conditions
//    if(hangingVertexDescription.getLastUpdateIterationParity() != _iterationParity) {
      hangingVertexDescription.setLastUpdateIterationParity(_iterationParity);
//    }
  }

  logTraceOutWith1Argument( "destroyHangingVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::Remesh::createInnerVertex(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createInnerVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
 
  fineGridVertex.setShouldRefine(false);
  fineGridVertex.resetSubcellsEraseVeto();

  for(int i = 0; i < TWO_POWER_D; i++) {
    fineGridVertex.setAdjacentCellDescriptionIndex(i, -1);
  }

  //Copy adjacent cell indices from former hanging vertex description, if available.
  tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> vertexPosition;
  for(int d = 0; d < DIMENSIONS; d++) {
    vertexPosition(d) = fineGridX(d);
  }
  vertexPosition(DIMENSIONS) = coarseGridVerticesEnumerator.getLevel() + 1;

  if( _vertexPositionToIndexMap.find(vertexPosition) != _vertexPositionToIndexMap.end() ) {
    VertexDescription& hangingVertexDescription = _vertexPositionToIndexMap[vertexPosition];

    for(int i = 0; i < TWO_POWER_D; i++) {
      fineGridVertex.setAdjacentCellDescriptionIndex(i, hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i));
    }
  } else {
    for(int i = 0; i < TWO_POWER_D; i++) {
      fineGridVertex.setAdjacentCellDescriptionIndex(i, -1);
    }
  }

  logTraceOutWith1Argument( "createInnerVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::Remesh::createBoundaryVertex(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createBoundaryVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
 
  fineGridVertex.setShouldRefine(false);
  fineGridVertex.resetSubcellsEraseVeto();

  for(int i = 0; i < TWO_POWER_D; i++) {
    fineGridVertex.setAdjacentCellDescriptionIndex(i, -1);
  }

  //Copy adjacent cell indices from former hanging vertex description, if available.
  tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> vertexPosition;
  for(int d = 0; d < DIMENSIONS; d++) {
    vertexPosition(d) = fineGridX(d);
  }
  vertexPosition(DIMENSIONS) = coarseGridVerticesEnumerator.getLevel() + 1;

  if( _vertexPositionToIndexMap.find(vertexPosition) != _vertexPositionToIndexMap.end() ) {
    VertexDescription& hangingVertexDescription = _vertexPositionToIndexMap[vertexPosition];

    for(int i = 0; i < TWO_POWER_D; i++) {
      fineGridVertex.setAdjacentCellDescriptionIndex(i, hangingVertexDescription.getIndicesOfAdjacentCellDescriptions(i));
    }
  } else {
    for(int i = 0; i < TWO_POWER_D; i++) {
      fineGridVertex.setAdjacentCellDescriptionIndex(i, -1);
    }
  }

  logTraceOutWith1Argument( "createBoundaryVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::Remesh::destroyVertex(
      const peanoclaw::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );

  //  Retrieve or create hanging vertex description
  tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> hangingVertexPosition;
  for(int d = 0; d < DIMENSIONS; d++) {
    hangingVertexPosition(d) = fineGridX(d);
  }
  hangingVertexPosition(DIMENSIONS) = coarseGridVerticesEnumerator.getLevel() + 1;

  if( _vertexPositionToIndexMap.find(hangingVertexPosition) == _vertexPositionToIndexMap.end() ) {
    VertexDescription vertexDescription;
    if(_iterationParity == peanoclaw::records::VertexDescription::EVEN) {
      vertexDescription.setLastUpdateIterationParity(peanoclaw::records::VertexDescription::ODD);
    } else {
      vertexDescription.setLastUpdateIterationParity(peanoclaw::records::VertexDescription::EVEN);
    }
    _vertexPositionToIndexMap[hangingVertexPosition] = vertexDescription;
  }

  VertexDescription& vertexDescription = _vertexPositionToIndexMap[hangingVertexPosition];
  vertexDescription.setTouched(true);

  //Copy adjacency information from destroyed vertex to hanging vertex description
  for(int i = 0; i < TWO_POWER_D; i++) {
    vertexDescription.setIndicesOfAdjacentCellDescriptions(i, fineGridVertex.getAdjacentCellDescriptionIndex(i));
  }

  logTraceOutWith1Argument( "destroyVertex(...)", fineGridVertex );
}


void peanoclaw::mappings::Remesh::createCell(
      peanoclaw::Cell&                 fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith6Arguments( "createCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, coarseGridVerticesEnumerator.toString(), fineGridPositionOfCell, fineGridVerticesEnumerator.getCellCenter() );

  //TODO unterweg debug
  #ifdef Parallel
//  if(tarch::la::equals(fineGridVerticesEnumerator.getVertexPosition(0)(0), 2.0/3.0)
//    && tarch::la::equals(fineGridVerticesEnumerator.getVertexPosition(0)(1), 23.0/27.0)
//  ) {
    std::cout << "Creating cell at " << fineGridVerticesEnumerator.getVertexPosition(0) << " on level " << fineGridVerticesEnumerator.getLevel() << " with size " << fineGridVerticesEnumerator.getCellSize() << " on rank " << tarch::parallel::Node::getInstance().getRank()
        << " " << fineGridCell.isInside() << std::endl;
//  }
  #endif
 
  //Initialise new Patch
  Patch fineGridPatch(
    fineGridVerticesEnumerator.getVertexPosition(),
    fineGridVerticesEnumerator.getCellSize(),
    _unknownsPerSubcell,
    _auxiliarFieldsPerSubcell,
    _defaultSubdivisionFactor,
    _defaultGhostLayerWidth,
    _initialTimestepSize,
    fineGridVerticesEnumerator.getLevel()
  );
  fineGridCell.setCellDescriptionIndex(fineGridPatch.getCellDescriptionIndex());

  if(fineGridCell.isLeaf()) {
    fineGridPatch.switchToVirtual();
    fineGridPatch.switchToLeaf();
  }

  //Transfer data from coarse to fine patch
  if(!coarseGridCell.isRoot()) {
    assertion4(coarseGridCell.getCellDescriptionIndex() > 0, coarseGridCell.getCellDescriptionIndex(), fineGridVerticesEnumerator.getCellSize(), fineGridVerticesEnumerator.getLevel(), fineGridVerticesEnumerator.getVertexPosition());
    Patch coarseGridPatch(
      coarseGridCell
    );
    assertion1(coarseGridPatch.getTimestepSize() >= 0.0, coarseGridPatch);

    if(!_isInitializing && (coarseGridPatch.isVirtual() || coarseGridPatch.isLeaf())) {
      fineGridPatch.setSkipNextGridIteration(2);

      fineGridPatch.setCurrentTime(coarseGridPatch.getCurrentTime());
      fineGridPatch.setTimestepSize(coarseGridPatch.getTimestepSize());
      fineGridPatch.setEstimatedNextTimestepSize(coarseGridPatch.getEstimatedNextTimestepSize());
      fineGridPatch.updateMinimalNeighborTimeConstraint(
        coarseGridPatch.getMinimalNeighborTimeConstraint(),
        coarseGridPatch.getCellDescriptionIndex()
      );

      //TODO unterweg debug
      #ifdef Parallel
      std::cout << "Interpolating: rank" << tarch::parallel::Node::getInstance().getRank() << " "   << coarseGridCell.isRemote(*_state, true, true) << " " << fineGridCell.isRemote(*_state, true, true) << " level:" << (fineGridVerticesEnumerator.getLevel()-1) << std::endl;
      #endif

      _numerics->interpolate(
        fineGridPatch.getSubdivisionFactor(),
        0,
        coarseGridPatch,
        fineGridPatch,
        false
      );
      _numerics->interpolate(
        fineGridPatch.getSubdivisionFactor(),
        0,
        coarseGridPatch,
        fineGridPatch,
        true
      );
    }
  }

  //Set indices on adjacent vertices and on this cell
  for(int i = 0; i < TWO_POWER_D; i++) {
    fineGridVertices[fineGridVerticesEnumerator(i)].setAdjacentCellDescriptionIndex(i, fineGridCell.getCellDescriptionIndex());
  }


  //TODO unterweg debug
  #ifdef Parallel
  if(tarch::la::equals(fineGridVerticesEnumerator.getVertexPosition(0)(0), 2.0/3.0)
    && tarch::la::equals(fineGridVerticesEnumerator.getVertexPosition(0)(1), 23.0/27.0)
  ) {
    std::cout << "Created cell on rank " << tarch::parallel::Node::getInstance().getRank() << ": " << fineGridPatch
        << " isInitializing: " << _isInitializing << std::endl;
    std::cout << fineGridPatch.toStringUNew() << std::endl;
    std::cout << fineGridPatch.toStringUOldWithGhostLayer() << std::endl;

    Patch coarseGridPatch(
      coarseGridCell
    );
    std::cout << "Coarse patch: " << coarseGridPatch << std::endl;
    std::cout << coarseGridPatch.toStringUNew() << std::endl;
    std::cout << coarseGridPatch.toStringUOldWithGhostLayer() << std::endl;
  }
  #endif

  logTraceOutWith2Arguments( "createCell(...)", fineGridCell, fineGridPatch );
}


void peanoclaw::mappings::Remesh::destroyCell(
      const peanoclaw::Cell&           fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "destroyCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  //TODO unterweg debug
  #ifdef Parallel
  std::cout << "Destroying cell at " << fineGridVerticesEnumerator.getVertexPosition(0) << " on level " << fineGridVerticesEnumerator.getLevel() << " with size " << fineGridVerticesEnumerator.getCellSize() << " on rank " << tarch::parallel::Node::getInstance().getRank()
      << " " << fineGridCell.isInside() << std::endl;
  #endif

  Patch finePatch(
    fineGridCell
  );

  //Fix minimal time patch
  if(_minimalTimePatch.isValid()
      && tarch::la::equals(finePatch.getPosition(), _minimalTimePatch.getPosition())
      && finePatch.getLevel() == _minimalTimePatch.getLevel()) {
    _minimalTimePatch = Patch();
  }
  if(_minimalTimePatchParent.isValid()
      && tarch::la::equals(finePatch.getPosition(), _minimalTimePatchParent.getPosition())
      && finePatch.getLevel() == _minimalTimePatchParent.getLevel()) {
    _minimalTimePatchParent = Patch();
  }

  //Create patch in parent cell if it doesn't exist
//  if(coarseGridVerticesEnumerator.getLevel() > 0) {
  if(!coarseGridCell.isRoot()) {
    CellDescription& coarseCellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(coarseGridCell.getCellDescriptionIndex()).at(0);

    //Fix timestep size
    coarseCellDescription.setTimestepSize(std::max(0.0, coarseCellDescription.getTimestepSize()));

    assertion2(coarseCellDescription.getUNewIndex() != -1, coarseGridCell, finePatch);
    assertion2(coarseCellDescription.getUOldIndex() != -1, coarseGridCell, finePatch);
 
    //Set indices on coarse adjacent vertices and fill adjacent ghostlayers
    for(int i = 0; i < TWO_POWER_D; i++) {
      fineGridVertices[fineGridVerticesEnumerator(i)].setAdjacentCellDescriptionIndex(i, coarseGridCell.getCellDescriptionIndex());
    }

    //Skip update for coarse patch in next grid iteration
    coarseCellDescription.setSkipGridIterations(1);

    //Set demanded mesh width for coarse cell to coarse cell size. Otherwise
    //the coarse patch might get refined immediately.
    coarseCellDescription.setDemandedMeshWidth(coarseGridVerticesEnumerator.getCellSize()(0) / coarseCellDescription.getSubdivisionFactor()(0));
  } else {
    for(int i = 0; i < TWO_POWER_D; i++) {
        fineGridVertices[fineGridVerticesEnumerator(i)].setAdjacentCellDescriptionIndex(i, -1);
    }
  }

  //Delete patch data and description from this cell
  assertion(fineGridCell.getCellDescriptionIndex() != -1);
  CellDescription& fineCellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(fineGridCell.getCellDescriptionIndex()).at(0);
  if(fineCellDescription.getUNewIndex() != -1) {
    //TODO unterweg debug
    peano::heap::Heap<Data>::getInstance().deleteData(fineCellDescription.getUNewIndex());
    fineCellDescription.setUNewIndex(-1);
  }
  if(fineCellDescription.getUOldIndex() != -1) {
    //TODO unterweg debug
    peano::heap::Heap<Data>::getInstance().deleteData(fineCellDescription.getUOldIndex());
    fineCellDescription.setUOldIndex(-1);
  }
  if(fineCellDescription.getAuxIndex() != -1) {
    //TODO unterweg debug
    peano::heap::Heap<Data>::getInstance().deleteData(fineCellDescription.getAuxIndex());
    fineCellDescription.setAuxIndex(-1);
  }
  int cellDescriptionIndex = fineGridCell.getCellDescriptionIndex();
  fineCellDescription.setCellDescriptionIndex(-1);
  //TODO unterweg debug
  assertion(cellDescriptionIndex > 1);
  peano::heap::Heap<CellDescription>::getInstance().deleteData(cellDescriptionIndex);

  logTraceOutWith1Argument( "destroyCell(...)", fineGridCell );
}

#ifdef Parallel
void peanoclaw::mappings::Remesh::mergeWithNeighbour(
  peanoclaw::Vertex&  vertex,
  const peanoclaw::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );

  assertionEquals(vertex.isInside(), neighbour.isInside());
  assertionEquals(vertex.isBoundary(), neighbour.isBoundary());

  int localRank = tarch::parallel::Node::getInstance().getRank();
  tarch::la::Vector<TWO_POWER_D, int> localVertexRanks = vertex.getAdjacentRanks();
  tarch::la::Vector<TWO_POWER_D, int> neighbourVertexRanks = neighbour.getAdjacentRanks();

  for(int i = TWO_POWER_D-1; i >= 0; i--) {
    tarch::la::Vector<DIMENSIONS,double> patchPosition = fineGridX + tarch::la::multiplyComponents(fineGridH, peano::utils::dDelinearised(i, 2).convertScalar<double>() - 1.0);
    peanoclaw::parallel::NeighbourCommunicator communicator(fromRank, patchPosition, level);
    int localAdjacentCellDescriptionIndex = vertex.getAdjacentCellDescriptionIndexInPeanoOrder(i);
    int remoteAdjacentCellDescriptionIndex = neighbour.getAdjacentCellDescriptionIndexInPeanoOrder(i);

    assertion4(
      localAdjacentCellDescriptionIndex != -1
      || vertex.isAdjacentToRemoteRank(),
      localAdjacentCellDescriptionIndex,
      remoteAdjacentCellDescriptionIndex,
      patchPosition,
      level
    );

    if(neighbourVertexRanks(i) == fromRank && remoteAdjacentCellDescriptionIndex != -1) {

      if(localAdjacentCellDescriptionIndex == -1) {
        //Create outside patch
        Patch outsidePatch(
          patchPosition,
          fineGridH,
          _unknownsPerSubcell,
          _auxiliarFieldsPerSubcell,
          _defaultSubdivisionFactor,
          _defaultGhostLayerWidth,
          _initialTimestepSize,
          level
        );
        localAdjacentCellDescriptionIndex = outsidePatch.getCellDescriptionIndex();
        vertex.setAdjacentCellDescriptionIndexInPeanoOrder(i, localAdjacentCellDescriptionIndex);
      }

      communicator.receivePatch(localAdjacentCellDescriptionIndex);
    } else {
      //Receive dummy message
      communicator.receivePaddingPatch();
    }
  }

//  if(tarch::la::equals(fineGridX, 1.0/3.0)
//    && level == 3) {
//    std::cout << "Merged vertex on rank " << tarch::parallel::Node::getInstance().getRank() << " from rank " << fromRank << std::endl;
//    for (int i = 0; i < TWO_POWER_D; i++) {
//      std::cout << " cellDescription(" << i << ")=" << vertex.getAdjacentCellDescriptionIndex(i);
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < TWO_POWER_D; i++) {
//      std::cout << " n.cellDescription(" << i << ")=" << neighbour.getAdjacentCellDescriptionIndex(i);
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < TWO_POWER_D; i++) {
//      std::cout << " n.adjacency(" << i << ")=" << neighbourVertexRanks(i);
//    }
//    std::cout << std::endl;
//  }

  logTraceOut( "mergeWithNeighbour(...)" );
}

void peanoclaw::mappings::Remesh::prepareSendToNeighbour(
  peanoclaw::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
  logTraceInWith3Arguments( "prepareSendToNeighbour(...)", vertex, toRank, level );
  int localRank = tarch::parallel::Node::getInstance().getRank();
  tarch::la::Vector<TWO_POWER_D, int> localVertexRanks = vertex.getAdjacentRanks();

  for(int i = 0; i < TWO_POWER_D; i++) {
    #ifdef Asserts
    tarch::la::Vector<DIMENSIONS,double> patchPosition = x + tarch::la::multiplyComponents(h, peano::utils::dDelinearised(i, 2).convertScalar<double>() - 1.0);
    #else
    tarch::la::Vector<DIMENSIONS,double> patchPosition(0.0);
    #endif
    peanoclaw::parallel::NeighbourCommunicator communicator(toRank, patchPosition, level);

    int adjacentCellDescriptionIndex = vertex.getAdjacentCellDescriptionIndexInPeanoOrder(i);
    if(localVertexRanks(i) == localRank && adjacentCellDescriptionIndex != -1) {
      communicator.sendPatch(adjacentCellDescriptionIndex);
    } else {
      //Send dummy message to avoid problems with heap
      communicator.sendPaddingPatch();
    }
  }

  logTraceOut( "prepareSendToNeighbour(...)" );
}

void peanoclaw::mappings::Remesh::prepareCopyToRemoteNode(
  peanoclaw::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
  logTraceInWith2Arguments( "prepareCopyToRemoteNode(...)", localVertex, toRank);
  // @todo Insert your code here
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void peanoclaw::mappings::Remesh::prepareCopyToRemoteNode(
  peanoclaw::Cell&  localCell,
  int  toRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                       level
) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank, cellCentre, cellSize, level );
  std::cout << "[prepareCopyToRemoteNode] sending to " << toRank << " next item @ " << cellCentre << " on level " << level << std::endl;

//  if(tarch::la::equals(cellSize, 1.0/9.0)
//    && tarch::la::equals(cellCentre(0) - cellSize(0)*0.5, 1.0/3.0)
//    && tarch::la::equals(cellCentre(1) - cellSize(1)*0.5, 0.0)) {
    std::cout << "Copying to remote from " << tarch::parallel::Node::getInstance().getRank()
        << " to " << toRank
        << " at " << (cellCentre - cellSize * 0.5)
        << " on level " << level
        << " isInside " << localCell.isInside() << std::endl;
//  }

  if(localCell.isInside()) {
    peanoclaw::parallel::MasterWorkerAndForkJoinCommunicator communicator(toRank, cellCentre, level, true);
    communicator.sendPatch(localCell.getCellDescriptionIndex());
  }
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void peanoclaw::mappings::Remesh::mergeWithRemoteDataDueToForkOrJoin(
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

void peanoclaw::mappings::Remesh::mergeWithRemoteDataDueToForkOrJoin(
  peanoclaw::Cell&  localCell,
  const peanoclaw::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                       level
) {
  logTraceInWith6Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localCell, masterOrWorkerCell, fromRank, cellCentre, cellSize, level );
  logInfo("", "[mergeWithRemoteDataDueToForkOrJoin] receiving next item @ " << cellCentre << " on level " << level << " localCell.isRemote: " << localCell.isRemote(*_state, false, false));

  //TODO unterweg debug
//  if(tarch::la::equals(cellSize, 1.0/3.0)
//    && tarch::la::equals(cellCentre(0) - cellSize(0)*0.5, 0.0/9.0)
//    && tarch::la::equals(cellCentre(1) - cellSize(1)*0.5, 0.0/9.0)) {
    std::cout << "Merged during fork on rank " << tarch::parallel::Node::getInstance().getRank()
        << " at " << (cellCentre - cellSize*0.5) << " on level " << level << " "
        << " inside:"<< localCell.isInside() << ", remote:" << localCell.isRemote(*_state, false, false) << std::endl;
//  }

  if(localCell.isInside()) {
    peanoclaw::parallel::MasterWorkerAndForkJoinCommunicator communicator(fromRank, cellCentre, level, true);

    if(localCell.isRemote(*_state, false, false)) {
      assertionEquals2(localCell.getCellDescriptionIndex(), 0, cellCentre, level);
      Patch temporaryPatch(
        cellCentre - cellSize * 0.5,
        cellSize,
        _unknownsPerSubcell,
        _auxiliarFieldsPerSubcell,
        _defaultSubdivisionFactor,
        _defaultGhostLayerWidth,
        _initialTimestepSize,
        level
      );

      communicator.receivePatch(temporaryPatch.getCellDescriptionIndex());

      if(temporaryPatch.isLeaf()) {
        temporaryPatch.switchToVirtual();
      }
      if(temporaryPatch.isVirtual()) {
        temporaryPatch.switchToNonVirtual();
      }
      peano::heap::Heap<CellDescription>::getInstance().deleteData(temporaryPatch.getCellDescriptionIndex());
    } else {
      assertion(localCell.getCellDescriptionIndex() != -1);
      assertion2(localCell.getCellDescriptionIndex() != 0, cellCentre, level);

      Patch localPatch(peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0));
      //TODO unterweg: May virtual patches be received during a fork/join
//      _gridLevelTransfer->updatePatchStateDuringForkOrJoin(
//        localPatch,
//        masterOrWorkerCell.getCellDescriptionIndex() != -1
//      );

      assertion2((!localPatch.isLeaf() && !localPatch.isVirtual()) || localPatch.getUNewIndex() > 0, cellCentre, level);

      if(tarch::la::equals(cellSize, 1.0/3.0)
        && tarch::la::equals(cellCentre(0) - cellSize(0)*0.5, 0.0)
        && tarch::la::equals(cellCentre(1) - cellSize(1)*0.5, 0.0)) {
        std::cout << "Merged during fork on rank " << tarch::parallel::Node::getInstance().getRank() << ": " << localPatch << std::endl;
      }

      communicator.receivePatch(localPatch.getCellDescriptionIndex());
      localPatch.loadCellDescription(localCell.getCellDescriptionIndex());

      assertionMsg(!localPatch.isVirtual(), "Receiving virtual patches during a fork is currently not supported!");

      //TODO unterweg debug
      localPatch = Patch(peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0));
      if(tarch::la::equals(cellSize, 1.0/3.0)
        && tarch::la::equals(cellCentre(0) - cellSize(0)*0.5, 0.0)
        && tarch::la::equals(cellCentre(1) - cellSize(1)*0.5, 0.0)) {
        std::cout << "Merged during fork on rank " << tarch::parallel::Node::getInstance().getRank() << ": " << localPatch << std::endl;
        std::cout << localPatch.toStringUNew() << std::endl;
      }
    }
  }

  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void peanoclaw::mappings::Remesh::prepareSendToWorker(
  peanoclaw::Cell&                 fineGridCell,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  logTraceInWith7Arguments( "prepareSendToWorker(...)", fineGridCell, fineGridVerticesEnumerator.toString(), fineGridVerticesEnumerator.getVertexPosition(0), coarseGridCell, coarseGridVerticesEnumerator.toString(), fineGridPositionOfCell, worker );
 
  std::cout << "[prepareSendToWorker] sending to worker " << worker << ": next item @ " << fineGridVerticesEnumerator.getCellCenter() << " on level " << fineGridVerticesEnumerator.getLevel() << " size " << fineGridVerticesEnumerator.getCellSize() << " inside:" << fineGridCell.isInside() << std::endl;

  if(fineGridCell.isInside()){
    peanoclaw::parallel::MasterWorkerAndForkJoinCommunicator communicator(worker, fineGridVerticesEnumerator.getCellCenter(), fineGridVerticesEnumerator.getLevel(), false);
    communicator.sendPatch(fineGridCell.getCellDescriptionIndex());

//    dataExchange_send(fineGridCell, worker, fineGridVerticesEnumerator.getCellCenter(), fineGridVerticesEnumerator.getCellSize(), fineGridVerticesEnumerator.getLevel(), false);

  }

  logTraceOut( "prepareSendToWorker(...)" );
}

void peanoclaw::mappings::Remesh::prepareSendToMaster(
  peanoclaw::Cell&                       localCell,
  peanoclaw::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator,
  const peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const peanoclaw::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
  logTraceInWith3Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString(), verticesEnumerator.getVertexPosition(0) );
  
  std::cout << "[prepareSendToMaster] sending to master: next item @ " << verticesEnumerator.getCellCenter() << " on level " << verticesEnumerator.getLevel() << std::endl;

  //TODO unterweg debug
  if(tarch::la::equals(verticesEnumerator.getVertexPosition(0), tarch::la::Vector<DIMENSIONS, double>(0.0))
    && verticesEnumerator.getLevel() == 1) {
    std::cout << "Send to Master: Updating Patchstate on rank " << tarch::parallel::Node::getInstance().getRank() << std::endl;
  }

  int toRank = tarch::parallel::NodePool::getInstance().getMasterRank();
  if(localCell.isInside()){
    peanoclaw::parallel::MasterWorkerAndForkJoinCommunicator communicator(toRank, verticesEnumerator.getCellCenter(), verticesEnumerator.getLevel(), false);
    communicator.sendPatch(localCell.getCellDescriptionIndex());
  }

  Patch localPatch(
    localCell
  );
  //TODO unterweg debug
//  _gridLevelTransfer->updatePatchStateAfterStepUp(
//    localPatch,
//    localCell.isLeaf()
//  );

  logTraceOut( "prepareSendToMaster(...)" );
}


void peanoclaw::mappings::Remesh::mergeWithMaster(
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
  logTraceInWith7Arguments( "mergeWithMaster(...)", workerGridCell, fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, coarseGridVerticesEnumerator.toString(), fineGridPositionOfCell, worker );

  masterState.updateGlobalTimeIntervals(
        workerState.getStartMaximumGlobalTimeInterval(),
        workerState.getEndMaximumGlobalTimeInterval(),
        workerState.getStartMinimumGlobalTimeInterval(),
        workerState.getEndMinimumGlobalTimeInterval()
  );

  masterState.updateMinimalTimestep(workerState.getMinimalTimestep());

  bool allPatchesEvolvedToGlobalTimestep = workerState.getAllPatchesEvolvedToGlobalTimestep();
  allPatchesEvolvedToGlobalTimestep &= masterState.getAllPatchesEvolvedToGlobalTimestep();

  masterState.setAllPatchesEvolvedToGlobalTimestep(allPatchesEvolvedToGlobalTimestep);

  std::cout << "merging with Master" << std::endl;

  tarch::la::Vector<DIMENSIONS,double>  cellCenter = workerEnumerator.getCellCenter();
  tarch::la::Vector<DIMENSIONS,double>  cellSize = workerEnumerator.getCellSize();
  int level = workerEnumerator.getLevel();

  if(fineGridCell.isInside()) {
    peanoclaw::parallel::MasterWorkerAndForkJoinCommunicator communicator(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getLevel(),
      false
    );
    communicator.receivePatch(fineGridCell.getCellDescriptionIndex());

//    dataExchange_receive(fineGridCell, workerGridCell, worker, cellCenter, cellSize, level, false);
  }

  logTraceOut( "mergeWithMaster(...)" );
}


void peanoclaw::mappings::Remesh::receiveDataFromMaster(
  peanoclaw::Cell&                        receivedCell,
  peanoclaw::Vertex *                     receivedVertices,
  const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
  peanoclaw::Vertex * const               receivedCoarseGridVertices,
  const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
  peanoclaw::Cell&                        receivedCoarseGridCell,
  peanoclaw::Vertex * const               workersCoarseGridVertices,
  const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
  peanoclaw::Cell&                        workersCoarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
  logTraceInWith2Arguments( "receiveDataFromMaster(...)", receivedCell.toString(), receivedVerticesEnumerator.toString() );
  
  //TODO unterweg debug
  std::cout << "Receiving on " << tarch::parallel::Node::getInstance().getRank() << " at remote: " << receivedVerticesEnumerator.getVertexPosition(tarch::la::Vector<DIMENSIONS,int>(0))
      << " size " << receivedVerticesEnumerator.getCellSize()
      << " local: " << workersCoarseGridVerticesEnumerator.getVertexPosition(tarch::la::Vector<DIMENSIONS,int>(0))
      << " size " << workersCoarseGridVerticesEnumerator.getCellSize()
      << " finePosition " << fineGridPositionOfCell
      << std::endl;

  if(receivedCell.isInside()) {
    peanoclaw::parallel::MasterWorkerAndForkJoinCommunicator communicator(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      receivedVerticesEnumerator.getCellCenter(),
      receivedVerticesEnumerator.getLevel(),
      false
    );

    int temporaryCellDescriptionIndex = peano::heap::Heap<CellDescription>::getInstance().createData();
    CellDescription temporaryCellDescription;
    temporaryCellDescription.setUNewIndex(-1);
    temporaryCellDescription.setUOldIndex(-1);
    temporaryCellDescription.setAuxIndex(-1);
    peano::heap::Heap<CellDescription>::getInstance().getData(temporaryCellDescriptionIndex).push_back(temporaryCellDescription);
    receivedCell.setCellDescriptionIndex(temporaryCellDescriptionIndex);

    communicator.receivePatch(temporaryCellDescriptionIndex);
  }

  logTraceOut( "receiveDataFromMaster(...)" );
}


void peanoclaw::mappings::Remesh::mergeWithWorker(
  peanoclaw::Cell&           localCell, 
  const peanoclaw::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );

  //TODO unterweg debug
  if(tarch::la::equals(cellCentre, tarch::la::Vector<DIMENSIONS, double>(0.5)) && level == 1) {
    std::cout << "Merging with Worker: Updating Patchstate on rank " << tarch::parallel::Node::getInstance().getRank() << std::endl;
  }

  if(receivedMasterCell.isInside()) {
    CellDescription& localCellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0);
    CellDescription& masterCellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(receivedMasterCell.getCellDescriptionIndex()).at(0);

    Patch localPatch(localCellDescription);
    //TODO unterweg debug
//    _gridLevelTransfer->updatePatchStateBeforeStepDown(
//      localPatch,
//      masterCellDescription.getUNewIndex() != -1
//    );

    //Delete current content of patch
    if(masterCellDescription.getUNewIndex() != -1) {
      assertion(localPatch.isVirtual() || localPatch.isLeaf());
      peano::heap::Heap<Data>::getInstance().deleteData(localCellDescription.getUNewIndex());
      peano::heap::Heap<Data>::getInstance().deleteData(localCellDescription.getUOldIndex());
      if(localCellDescription.getAuxIndex() != -1) {
        peano::heap::Heap<Data>::getInstance().deleteData(localCellDescription.getAuxIndex());
      }

      localCellDescription.setUNewIndex(masterCellDescription.getUNewIndex());
      localCellDescription.setUOldIndex(masterCellDescription.getUOldIndex());
      localCellDescription.setAuxIndex(masterCellDescription.getAuxIndex());
    }

    peano::heap::Heap<CellDescription>::getInstance().deleteData(receivedMasterCell.getCellDescriptionIndex());
  } else {
    CellDescription& localCellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0);
    Patch localPatch(localCellDescription);
    //TODO unterweg debug
//    _gridLevelTransfer->updatePatchStateBeforeStepDown(
//      localPatch,
//      false
//    );
  }

  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}


void peanoclaw::mappings::Remesh::mergeWithWorker(
  peanoclaw::Vertex&        localVertex,
  const peanoclaw::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localVertex.toString(), receivedMasterVertex.toString() );
  // @todo Insert your code here
  logTraceOutWith1Argument( "mergeWithWorker(...)", localVertex.toString() );
}
#endif

void peanoclaw::mappings::Remesh::touchVertexFirstTime(
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


void peanoclaw::mappings::Remesh::touchVertexLastTime(
      peanoclaw::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexLastTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
 
  //Regain 2-irregularity if necessary
  if(fineGridVertex.getRefinementControl() == peanoclaw::Vertex::Records::Refined
      || fineGridVertex.getRefinementControl() == peanoclaw::Vertex::Records::Refining) {
    tarch::la::Vector<DIMENSIONS, int> coarseGridPositionOfVertex(0);
    for(int d = 0; d < DIMENSIONS; d++) {
      if(fineGridPositionOfVertex(d) > 1) {
        coarseGridPositionOfVertex(d) = 1;
      }
    }

    peanoclaw::Vertex& coarseVertex = coarseGridVertices[coarseGridVerticesEnumerator(coarseGridPositionOfVertex)];
    if(coarseVertex.getRefinementControl() == peanoclaw::Vertex::Records::Unrefined
        && !coarseVertex.isHangingNode()) {
      coarseVertex.refine();
    }
  }

  logTraceOutWith1Argument( "touchVertexLastTime(...)", fineGridVertex );
}


void peanoclaw::mappings::Remesh::enterCell(
      peanoclaw::Cell&                 fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith5Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, coarseGridVerticesEnumerator.toString(), fineGridPositionOfCell );

  //TODO unterweg debug
  #ifdef Parallel
  if(tarch::la::equals(fineGridVerticesEnumerator.getVertexPosition(0), tarch::la::Vector<DIMENSIONS, double>(0.0))
    && fineGridVerticesEnumerator.getLevel() == 1) {
    std::cout << "enterCell: Updating Patchstate on rank " << tarch::parallel::Node::getInstance().getRank() << std::endl;
  }
  #endif

  Patch patch(
    fineGridCell
  );

  assertionEquals1(patch.getLevel(), fineGridVerticesEnumerator.getLevel(), patch.toString());

  _gridLevelTransfer->updatePatchStateBeforeStepDown(
    patch,
    fineGridVertices,
    fineGridVerticesEnumerator,
    _isInitializing,
    _state->isInvolvedInJoinOrFork()
  );

  _gridLevelTransfer->stepDown(
    coarseGridCell.isRoot() ? -1 : coarseGridCell.getCellDescriptionIndex(),
    patch,
    fineGridVertices,
    fineGridVerticesEnumerator
  );

#ifdef Asserts
if(patch.isLeaf() && !fineGridCell.isLeaf()) {
  bool isRefining = false;
  for(int i = 0; i < TWO_POWER_D; i++) {
    if(fineGridVertices[fineGridVerticesEnumerator(i)].getRefinementControl() == peanoclaw::records::Vertex::Refining) {
      isRefining = true;
    }
  }
  assertion1(isRefining, patch);
}
#endif

  logTraceOutWith2Arguments( "enterCell(...)", fineGridCell, patch );
}


void peanoclaw::mappings::Remesh::leaveCell(
      peanoclaw::Cell&           fineGridCell,
      peanoclaw::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "leaveCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  //TODO unterweg debug
  #ifdef Parallel
  if(tarch::la::equals(fineGridVerticesEnumerator.getVertexPosition(0), tarch::la::Vector<DIMENSIONS, double>(0.0))
    && fineGridVerticesEnumerator.getLevel() == 1) {
    std::cout << "leaveCell: Updating Patchstate on rank " << tarch::parallel::Node::getInstance().getRank() << std::endl;
  }
  #endif

  assertion(fineGridCell.isInside());
  assertion(fineGridCell.getCellDescriptionIndex() != -1);
  assertion(coarseGridCell.getCellDescriptionIndex() != -1);

  Patch finePatch(
    fineGridCell
  );

  if(!coarseGridCell.isRoot()) {
    Patch coarsePatch(
      coarseGridCell
    );

    //Copy vertices
    peanoclaw::Vertex vertices[TWO_POWER_D];
    for(int i = 0; i < TWO_POWER_D; i++) {
      vertices[i] = fineGridVertices[fineGridVerticesEnumerator(i)];
    }
  }

  _gridLevelTransfer->stepUp(
    coarseGridCell.getCellDescriptionIndex(),
    finePatch,
    fineGridCell.isLeaf(),
    fineGridVertices,
    fineGridVerticesEnumerator
  );

  _gridLevelTransfer->updatePatchStateAfterStepUp(
    finePatch,
    fineGridCell.isLeaf()
  );

  assertionEquals1(finePatch.isLeaf(), fineGridCell.isLeaf(), finePatch);
  assertionEquals1(finePatch.getLevel(), fineGridVerticesEnumerator.getLevel(), finePatch.toString());

  //Todo unterweg debug
  if(finePatch.isLeaf() && (finePatch.getCurrentTime() + finePatch.getTimestepSize() < _minimalPatchTime)) {
    _minimalPatchTime = finePatch.getCurrentTime() + finePatch.getTimestepSize();
    _minimalTimePatch = finePatch;
    _minimalPatchCoarsening = peano::grid::aspects::VertexStateAnalysis::doesOneVertexCarryRefinementFlag
                            (
                              coarseGridVertices,
                              coarseGridVerticesEnumerator,
                              peanoclaw::records::Vertex::Erasing
                            );
    _minimalPatchIsAllowedToAdvanceInTime = finePatch.isAllowedToAdvanceInTime();
    _minimalPatchShouldSkipGridIteration = finePatch.shouldSkipNextGridIteration();

    if(coarseGridCell.getCellDescriptionIndex() > 0) {
      Patch coarsePatch(coarseGridCell);
      _minimalTimePatchParent = coarsePatch;
    }
  }

  for(int i = 0; i < TWO_POWER_D; i++) {
    fineGridVertices[fineGridVerticesEnumerator(i)].setAdjacentCellDescriptionIndex(i, fineGridCell.getCellDescriptionIndex());
  }

  logTraceOutWith1Argument( "leaveCell(...)", fineGridCell );
}


void peanoclaw::mappings::Remesh::beginIteration(
  peanoclaw::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );
 
  _unknownsPerSubcell       = solverState.getUnknownsPerSubcell();

  _auxiliarFieldsPerSubcell = solverState.getAuxiliarFieldsPerSubcell();

  _defaultSubdivisionFactor = solverState.getDefaultSubdivisionFactor();

  _defaultGhostLayerWidth   = solverState.getDefaultGhostLayerWidth();

  _initialTimestepSize      = solverState.getInitialTimestepSize();

  _numerics                   = solverState.getNumerics();

  _domainOffset             = solverState.getDomainOffset();

  _domainSize               = solverState.getDomainSize();

  if(_iterationParity == peanoclaw::records::VertexDescription::EVEN) {
    _iterationParity = peanoclaw::records::VertexDescription::ODD;
  } else {
    _iterationParity = peanoclaw::records::VertexDescription::EVEN;
  }

  _gridLevelTransfer = new peanoclaw::interSubgridCommunication::GridLevelTransfer(solverState.useDimensionalSplitting(), *_numerics);

  _initialMinimalMeshWidth = solverState.getInitialMinimalMeshWidth();

  _additionalLevelsForPredefinedRefinement = solverState.getAdditionalLevelsForPredefinedRefinement();
  _isInitializing = solverState.getIsInitializing();
  _averageGlobalTimeInterval = (solverState.getStartMaximumGlobalTimeInterval() + solverState.getEndMaximumGlobalTimeInterval()) / 2.0;

  //TODO unterweg debug
  _minimalPatchTime = std::numeric_limits<double>::max();
  _minimalTimePatch = Patch();
  _minimalTimePatchParent = Patch();

  _useDimensionalSplitting = solverState.useDimensionalSplitting();

  //Reset touched for all hanging vertex descriptions
  std::map<tarch::la::Vector<DIMENSIONS_PLUS_ONE,double> , VertexDescription, tarch::la::VectorCompare<DIMENSIONS_PLUS_ONE> >::iterator i = _vertexPositionToIndexMap.begin();
  while(i != _vertexPositionToIndexMap.end()) {
    if(i->second.getTouched()) {
      i->second.setTouched(false);
      i++;
    } else {
      _vertexPositionToIndexMap.erase(i++);
    }
  }
 
  peano::heap::Heap<peanoclaw::records::Data>::getInstance().startToSendOrReceiveHeapData(solverState.isTraversalInverted());
  peano::heap::Heap<CellDescription>::getInstance().startToSendOrReceiveHeapData(solverState.isTraversalInverted());

  _state = &solverState;

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void peanoclaw::mappings::Remesh::endIteration(
  peanoclaw::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
 
  delete _gridLevelTransfer;

  //Todo unterweg debug
  _minimalTimePatch.reloadCellDescription();
  _minimalTimePatchParent.reloadCellDescription();
  if(_minimalTimePatch.isValid()) {
//    std::cout << "Minimal time patch: " << _minimalTimePatch << std::endl;
  //  std::cout << "Minimal time patch parent: " << _minimalTimePatchParent << std::endl;

//    if(_minimalTimePatch.getConstrainingNeighborIndex() != -1) {
//      Patch constrainingPatch(peano::heap::Heap<CellDescription>::getInstance().getData(_minimalTimePatch.getConstrainingNeighborIndex()).at(0));
//      std::cout << "Contrained by " << constrainingPatch << std::endl;
//    }
  }

  peano::heap::Heap<peanoclaw::records::Data>::getInstance().finishedToSendOrReceiveHeapData();
  peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().finishedToSendOrReceiveHeapData();

  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void peanoclaw::mappings::Remesh::descend(
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


void peanoclaw::mappings::Remesh::ascend(
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
