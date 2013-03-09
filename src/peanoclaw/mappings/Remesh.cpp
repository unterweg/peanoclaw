#include "peanoclaw/mappings/Remesh.h"

#include "peanoclaw/interSubgridCommunication/GridLevelTransfer.h"
#include "peanoclaw/Numerics.h"

#include "peano/heap/Heap.h"

peanoclaw::records::VertexDescription::IterationParity peanoclaw::mappings::Remesh::_iterationParity
 = peanoclaw::records::VertexDescription::EVEN;

std::map<tarch::la::Vector<DIMENSIONS_PLUS_ONE,double> , peanoclaw::mappings::Remesh::VertexDescription, tarch::la::VectorCompare<DIMENSIONS_PLUS_ONE> >   peanoclaw::mappings::Remesh::_vertexPositionToIndexMap;

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   peanoclaw::mappings::Remesh::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WHOLE_TREE,false,false);
}


tarch::logging::Log                peanoclaw::mappings::Remesh::_log( "peanoclaw::mappings::Remesh" ); 


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
  logTraceInWith4Arguments( "createCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
 
  //Initialise Cell Description
  fineGridCell.setCellDescriptionIndex(peano::heap::Heap<CellDescription>::getInstance().createData());
  std::vector<CellDescription>& descriptions = peano::heap::Heap<CellDescription>::getInstance().getData(fineGridCell.getCellDescriptionIndex());
  descriptions.push_back(CellDescription());
  assertionEquals(peano::heap::Heap<CellDescription>::getInstance().getData(fineGridCell.getCellDescriptionIndex()).size(), 1);

  CellDescription& cellDescription = descriptions.at(0);
  cellDescription.setCellDescriptionIndex(fineGridCell.getCellDescriptionIndex());
  cellDescription.setGhostLayerWidth(_defaultGhostLayerWidth);
  cellDescription.setSubdivisionFactor(_defaultSubdivisionFactor);
  cellDescription.setUnknownsPerSubcell(_unknownsPerSubcell);
  cellDescription.setAuxiliarFieldsPerSubcell(_auxiliarFieldsPerSubcell);
  cellDescription.setTime(0.0);
  cellDescription.setTimestepSize(0.0);
  cellDescription.setEstimatedNextTimestepSize(_initialTimestepSize);
  cellDescription.setMinimalNeighborTime(std::numeric_limits<double>::max());
  cellDescription.setMinimalNeighborTimeConstraint(std::numeric_limits<double>::max());
  cellDescription.setLevel(fineGridVerticesEnumerator.getLevel());
  cellDescription.setPosition(fineGridVerticesEnumerator.getVertexPosition(0));
  cellDescription.setSize(fineGridVerticesEnumerator.getCellSize());
  cellDescription.setIsVirtual(false);
  cellDescription.setSkipNextGridIteration(false);
  cellDescription.setDemandedMeshWidth(fineGridVerticesEnumerator.getCellSize()(0) / _defaultSubdivisionFactor(0) * 3.0);
  cellDescription.setAgeInGridIterations(0);
  cellDescription.setUOldIndex(-1);
  cellDescription.setUNewIndex(-1);
  cellDescription.setAuxIndex(-1);
  cellDescription.setRestrictionLowerBounds(std::numeric_limits<double>::max());
  cellDescription.setRestrictionUpperBounds(-std::numeric_limits<double>::max());

  Patch fineGridPatch(
    fineGridVerticesEnumerator.getVertexPosition(0),
    fineGridVerticesEnumerator.getCellSize(),
    fineGridCell
  );

  if(fineGridCell.isLeaf()) {
    fineGridPatch.switchToVirtual();
    fineGridPatch.switchToLeaf();
  }

//  double minimalNeighborTime = std::numeric_limits<double>::max();
  //Transfer data from coarse to fine patch
  //TODO Kann man die Root-Zelle richtig setzen? Dann kann man diese Abfrage rausschmeißen.
  if(coarseGridVerticesEnumerator.getLevel() > 0) {
    Patch coarseGridPatch(
      coarseGridVerticesEnumerator.getVertexPosition(0),
      coarseGridVerticesEnumerator.getCellSize(),
      coarseGridCell
    );

    //if(coarseGridPatch.isLeaf() && !_isInitializing) {
    if(!_isInitializing) {
      cellDescription.setSkipNextGridIteration(true);
      cellDescription.setEstimatedNextTimestepSize(coarseGridPatch.getEstimatedNextTimestepSize() / 6.0);

      cellDescription.setTime(coarseGridPatch.getCurrentTime());
      cellDescription.setTimestepSize(coarseGridPatch.getTimestepSize());
      cellDescription.setEstimatedNextTimestepSize(coarseGridPatch.getEstimatedNextTimestepSize());
      cellDescription.setMinimalNeighborTimeConstraint(coarseGridPatch.getMinimalNeighborTimeConstraint());

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

  //Create patch in parent cell if it doesn't exist
  if(coarseGridVerticesEnumerator.getLevel() > 0) {
    CellDescription& coarseCellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(coarseGridCell.getCellDescriptionIndex()).at(0);
    Patch finePatch(
      fineGridVerticesEnumerator.getVertexPosition(0),
      fineGridVerticesEnumerator.getCellSize(),
      fineGridCell
    );

    //Fix timestep size
    coarseCellDescription.setTimestepSize(std::max(0.0, coarseCellDescription.getTimestepSize()));

    assertion2(coarseCellDescription.getUNewIndex() != -1, coarseGridCell, finePatch);
    assertion2(coarseCellDescription.getUOldIndex() != -1, coarseGridCell, finePatch);
 
    //Set indices on coarse adjacent vertices and fill adjacent ghostlayers
    for(int i = 0; i < TWO_POWER_D; i++) {
      fineGridVertices[fineGridVerticesEnumerator(i)].setAdjacentCellDescriptionIndex(i, coarseGridCell.getCellDescriptionIndex());
    }

    //Skip update for coarse patch in next grid iteration
    coarseCellDescription.setSkipNextGridIteration(true);

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
    peano::heap::Heap<Data>::getInstance().deleteData(fineCellDescription.getUNewIndex());
    fineCellDescription.setUNewIndex(-1);
  }
  if(fineCellDescription.getUOldIndex() != -1) {
    peano::heap::Heap<Data>::getInstance().deleteData(fineCellDescription.getUOldIndex());
    fineCellDescription.setUOldIndex(-1);
  }
  if(fineCellDescription.getAuxIndex() != -1) {
    peano::heap::Heap<Data>::getInstance().deleteData(fineCellDescription.getAuxIndex());
    fineCellDescription.setAuxIndex(-1);
  }
  int cellDescriptionIndex = fineGridCell.getCellDescriptionIndex();
  fineCellDescription.setCellDescriptionIndex(-1);
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
  // @todo Insert your code here
  logTraceOut( "mergeWithNeighbour(...)" );
}

void peanoclaw::mappings::Remesh::prepareSendToNeighbour(
  peanoclaw::Vertex&  vertex,
  int  toRank,
  int  level
) {
  logTraceInWith3Arguments( "prepareSendToNeighbour(...)", vertex, toRank, level );
  // @todo Insert your code here
  logTraceOut( "prepareSendToNeighbour(...)" );
}

void peanoclaw::mappings::Remesh::prepareCopyToRemoteNode(
  peanoclaw::Vertex&  localVertex,
  int  toRank
) {
  logTraceInWith2Arguments( "prepareCopyToRemoteNode(...)", localVertex, toRank );
  // @todo Insert your code here
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void peanoclaw::mappings::Remesh::prepareCopyToRemoteNode(
  peanoclaw::Cell&  localCell,
  int  toRank
) {
  logTraceInWith2Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank );
  // @todo Insert your code here
  if (!tarch::parallel::Node::getInstance().isGlobalMaster()) {
      //CellDescription* _cellDescription = &peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0);
      /*Patch patch(*_cellDescription);
      
      peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().sendData(localCell.getCellDescriptionIndex(), toRank, patch.getPosition(), patch.getLevel(), true);
     
      if (patch.getUNewIndex() != -1) {
         peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(patch.getUNewIndex(), toRank, patch.getPosition(), patch.getLevel(), true);
      }

      if (patch.getUOldIndex() != -1) {
         peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(patch.getUOldIndex(), toRank, patch.getPosition(), patch.getLevel(), true);
      }

      if (patch.getAuxIndex() != -1) {
          peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(patch.getAuxIndex(), toRank, patch.getPosition(), patch.getLevel(), true);
      }*/
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
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
  logTraceInWith3Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localCell, masterOrWorkerCell, fromRank );
  // @todo Insert your code here
  /*const int globalMasterRank = tarch::parallel::Node::getInstance().getGlobalMasterRank();
  if (fromRank != globalMasterRank) {
      peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().receiveData(fromRank, x, level, true);
      peano::heap::Heap<peanoclaw::records::Data>::getInstance().receiveData(fromRank, x, level, true);

      localCell.setCellDescriptionIndex(peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().createData());
      CellDescription* _localCellDescription = &peano::heap::Heap<CellDescription>::getInstance().getData(localCell.getCellDescriptionIndex()).at(0);
      CellDescription* _remoteCellDescription = &peano::heap::Heap<CellDescription>::getInstance().getData(masterOrWorkerCell.getCellDescriptionIndex()).at(0);
      *_localCellDescription = *_remoteCellDescription;
  }*/
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
  logTraceIn( "prepareSendToWorker(...)" );
  // @todo Insert your code here
#if 0
  //if (!tarch::parallel::Node::getInstance().isGlobalMaster()) {
      /*Patch finePatch(
        fineGridVerticesEnumerator.getVertexPosition(0),
        fineGridVerticesEnumerator.getCellSize(),
        fineGridCell
      );*/

      Patch coarsePatch(
        coarseGridVerticesEnumerator.getVertexPosition(0),
        coarseGridVerticesEnumerator.getCellSize(),
        coarseGridCell
      );


      /*if (finePatch.getUNewIndex() != -1) {
         peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(finePatch.getUNewIndex(), worker, fineGridVerticesEnumerator.getVertexPosition(0), fineGridVerticesEnumerator.getLevel(), true);
      }

      if (finePatch.getUOldIndex() != -1) {
         peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(finePatch.getUOldIndex(), worker, fineGridVerticesEnumerator.getVertexPosition(0), fineGridVerticesEnumerator.getLevel(), true);
      }

      if (finePatch.getAuxIndex() != -1) {
          peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(finePatch.getAuxIndex(), worker, fineGridVerticesEnumerator.getVertexPosition(0), fineGridVerticesEnumerator.getLevel(), true);
      }*/

      if (coarsePatch.getUNewIndex() != -1) {
         peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(coarsePatch.getUNewIndex(), worker, coarseGridVerticesEnumerator.getVertexPosition(0), coarseGridVerticesEnumerator.getLevel(), true);
      }

      if (coarsePatch.getUOldIndex() != -1) {
         peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(coarsePatch.getUOldIndex(), worker, coarseGridVerticesEnumerator.getVertexPosition(0), coarseGridVerticesEnumerator.getLevel(), true);
      }

      if (coarsePatch.getAuxIndex() != -1) {
          peano::heap::Heap<peanoclaw::records::Data>::getInstance().sendData(coarsePatch.getAuxIndex(), worker, coarseGridVerticesEnumerator.getVertexPosition(0), coarseGridVerticesEnumerator.getLevel(), true);
      }


      peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().sendData(coarseGridCell.getCellDescriptionIndex(), worker, coarseGridVerticesEnumerator.getVertexPosition(0), coarseGridVerticesEnumerator.getLevel(), true);

      //peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().sendData(fineGridCell.getCellDescriptionIndex(), worker, fineGridVerticesEnumerator.getVertexPosition(0), fineGridVerticesEnumerator.getLevel(), true);
  //}
#endif
  logTraceOut( "prepareSendToWorker(...)" );
}

void peanoclaw::mappings::Remesh::prepareSendToMaster(
  peanoclaw::Cell&     localCell,
  peanoclaw::Vertex *  vertices,
  const peano::grid::VertexEnumerator&  verticesEnumerator
) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );
  // @todo Insert your code here
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
  logTraceIn( "mergeWithMaster(...)" );
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

  logTraceOut( "mergeWithMaster(...)" );
}


void peanoclaw::mappings::Remesh::receiveDataFromMaster(
  peanoclaw::Cell&                    receivedCell, 
  peanoclaw::Vertex *                 receivedVertices,
  const peano::grid::VertexEnumerator&    verticesEnumerator
) {
  logTraceInWith2Arguments( "receiveDataFromMaster(...)", receivedCell.toString(), verticesEnumerator.toString() );
  // @todo Insert your code here
#if 0
  const int masterRank = tarch::parallel::NodePool::getInstance().getMasterRank();
  const int globalMasterRank = tarch::parallel::Node::getInstance().getGlobalMasterRank();

  if (masterRank != globalMasterRank) {
      peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().receiveData(masterRank, verticesEnumerator.getVertexPosition(), verticesEnumerator.getLevel(), true);
      peano::heap::Heap<peanoclaw::records::Data>::getInstance().receiveData(masterRank, verticesEnumerator.getVertexPosition(), verticesEnumerator.getLevel(), true);
      //Patch patch(verticesEnumerator.getVertexPosition(), verticesEnumerator.getCellSize(), receivedCell);
  }
#endif
  logTraceOut( "receiveDataFromMaster(...)" );
}


void peanoclaw::mappings::Remesh::mergeWithWorker(
  peanoclaw::Cell&           localCell, 
  const peanoclaw::Cell&     receivedMasterCell
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );
  // @todo Insert your code here
  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}


void peanoclaw::mappings::Remesh::mergeWithWorker(
  peanoclaw::Vertex&        localVertex,
  const peanoclaw::Vertex&  receivedMasterVertex
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
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  
  Patch patch(
    fineGridVerticesEnumerator.getVertexPosition(0),
    fineGridVerticesEnumerator.getCellSize(),
    fineGridCell
  );

  assertionEquals1(patch.getLevel(), fineGridVerticesEnumerator.getLevel(), patch.toString());

  _gridLevelTransfer->stepDown(
    coarseGridCell.getCellDescriptionIndex(),
    patch,
    fineGridVertices,
    fineGridVerticesEnumerator,
    _isInitializing
  );

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

  assertion(fineGridCell.getCellDescriptionIndex() != -1);
  assertion(coarseGridCell.getCellDescriptionIndex() != -1);

  Patch finePatch(
    fineGridVerticesEnumerator.getVertexPosition(0),
    fineGridVerticesEnumerator.getCellSize(),
    fineGridCell
  );

//  if(!fineGridCell.isLeaf()) {
//    finePatch.updateTimeIntervalFromMinimalFineGridTimeInterval();
//  }

  if(coarseGridVerticesEnumerator.getLevel() > 0) {
    Patch coarsePatch(coarseGridVerticesEnumerator.getVertexPosition(0), coarseGridVerticesEnumerator.getCellSize(), coarseGridCell);

    //Copy vertices
    peanoclaw::Vertex vertices[TWO_POWER_D];
    for(int i = 0; i < TWO_POWER_D; i++) {
      vertices[i] = fineGridVertices[fineGridVerticesEnumerator(i)];
    }
  }

  _gridLevelTransfer->stepUp(
    coarseGridCell.getCellDescriptionIndex(),
    finePatch, fineGridCell.isLeaf(),
    fineGridVertices,
    fineGridVerticesEnumerator
  );

  assertionEquals1(finePatch.getLevel(), fineGridVerticesEnumerator.getLevel(), finePatch.toString());

  //Todo unterweg debug
  if(finePatch.isLeaf() && (finePatch.getCurrentTime() + finePatch.getTimestepSize() < _minimalPatchTime)) {
    _minimalPatchTime = finePatch.getCurrentTime() + finePatch.getTimestepSize();
    _minimalTimePatch = finePatch;
  }

//  CellDescription& cellDescription = peano::heap::Heap<CellDescription>::getInstance().getData(fineGridCell.getCellDescriptionIndex()).at(0);

  // Delete indices on adjacent vertices, if the cell is not a leaf cell and, thus, contains no real patch.
//  if(!fineGridCell.isLeaf()) {
//    if(cellDescription.getUNewIndex() != -1) {
//      logDebug("enterCell(...)", "Deleting uNew for non-leaf cell at " << fineGridVerticesEnumerator.getCellCenter() << ", level " << fineGridVerticesEnumerator.getLevel() << ".");
//      peano::heap::Heap<Data>::getInstance().deleteData(cellDescription.getUNewIndex());
//      cellDescription.setUNewIndex(-1);
//    }
//    if(cellDescription.getUOldIndex() != -1) {
//      logDebug("enterCell(...)", "Deleting uOld for non-leaf cell at " << fineGridVerticesEnumerator.getCellCenter() << ", level " << fineGridVerticesEnumerator.getLevel() << ".");
//      peano::heap::Heap<Data>::getInstance().deleteData(cellDescription.getUOldIndex());
//      cellDescription.setUOldIndex(-1);
//    }
//  }

  for(int i = 0; i < TWO_POWER_D; i++) {
//    assertion(fineGridCell.getCellDescriptionIndex() == -1 || fineGridVertices[fineGridVerticesEnumerator(i)].getAdjacentCellDescriptionIndex(i) == -1);
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
 
  //peano::heap::Heap<peanoclaw::records::Data>::getInstance().startToSendOrReceiveHeapData(solverState.isTraversalInverted());
  //peano::heap::Heap<CellDescription>::getInstance().startToSendOrReceiveHeapData(solverState.isTraversalInverted());

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void peanoclaw::mappings::Remesh::endIteration(
  peanoclaw::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
 
  delete _gridLevelTransfer;

  //Todo unterweg debug
  /*std::cout << "Minimal time patch: " << _minimalTimePatch.toString() << std::endl;
  if(_minimalTimePatch.isValid()) {
    std::cout << "is allowed to advance: " << _minimalTimePatch.isAllowedToAdvanceInTime() << std::endl
        << "should skip next grid iteration: " << _minimalTimePatch.shouldSkipNextGridIteration() << std::endl;
  }*/
  
  //peano::heap::Heap<peanoclaw::records::Data>::getInstance().finishedToSendOrReceiveHeapData();
  //peano::heap::Heap<peanoclaw::records::CellDescription>::getInstance().finishedToSendOrReceiveHeapData();

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
