/*
 * VTKSubgridPlotter.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: kristof
 */
#include "peanoclaw/grid/plotter/VTKSubgridPlotter.h"

#include "peanoclaw/Patch.h"

peanoclaw::grid::plotter::VTKSubgridPlotter::VTKSubgridPlotter(
  std::string fileName,
  int unknownsPerSubcell,
  int parametersWithoutGhostlayerPerSubcell,
  int parametersWithGhostlayerPerSubcell,
  std::set<int> plotQ = std::set<int>(),
  std::set<int> plotParameterWithoutGhostlayer = std::set<int>(),
  std::set<int> plotParameterWithGhostlayer = std::set<int>(),
  bool plotMetainformation = true
) : _fileName(fileName),
    _vtkWriter(),
    _gap(0.015),
    _plotQ(plotQ),
    _plotParameterWithoutGhostlayer(plotParameterWithoutGhostlayer),
    _plotParameterWithGhostlayer(plotParameterWithGhostlayer),
    _plotMetainformation(plotMetainformation) {
  _vertex2IndexMap.clear();
  _vertexWriter     = _vtkWriter.createVertexWriter();
  _cellWriter       = _vtkWriter.createCellWriter();

  for(int i = 0; i < unknownsPerSubcell; i++) {
    if(_plotQ.empty() || _plotQ.find(i) != _plotQ.end()) {
      std::stringstream stringstream;
      stringstream << "q" << i;
      _cellQWriter.push_back( _vtkWriter.createCellDataWriter(stringstream.str(), 1) );
    }
  }
  for(int i = 0; i < parametersWithoutGhostlayerPerSubcell; i++) {
    if(_plotParameterWithoutGhostlayer.empty() || _plotParameterWithoutGhostlayer.find(i) != _plotParameterWithoutGhostlayer.end()) {
      std::stringstream stringstream;
      stringstream << "pWoG" << i;
      _cellParameterWithoutGhostlayerWriter.push_back( _vtkWriter.createCellDataWriter(stringstream.str(), 1) );
    }
  }
  for(int i = 0; i < parametersWithGhostlayerPerSubcell; i++) {
    if(_plotParameterWithGhostlayer.empty() || _plotParameterWithGhostlayer.find(i) != _plotParameterWithGhostlayer.end()) {
      std::stringstream stringstream;
      stringstream << "pWG" << i;
      _cellParameterWithGhostlayerWriter.push_back( _vtkWriter.createCellDataWriter(stringstream.str(), 1) );
    }
  }

  if(_plotMetainformation) {
    _cellSubdivisionFactorWriter = _vtkWriter.createCellDataWriter("SubdivisionFactor",1);
    _cellGhostLayerWidthWriter   = _vtkWriter.createCellDataWriter("GhostlayerWidth",1);
    _cellTimeOldWriter           = _vtkWriter.createCellDataWriter("timeOld", 1);
    _cellTimeNewWriter           = _vtkWriter.createCellDataWriter("timeNew", 1);
    _cellDemandedMeshWidthWriter = _vtkWriter.createCellDataWriter("demandedMeshWidth", 1);
    _cellAgeWriter               = _vtkWriter.createCellDataWriter("age", 1);
    #ifdef Parallel
    _cellRankWriter               = _vtkWriter.createCellDataWriter("rank", 1);
    #endif
  }
}

peanoclaw::grid::plotter::VTKSubgridPlotter::~VTKSubgridPlotter() {
  close();

  _vtkWriter.writeToFile( _fileName );

  delete _vertexWriter;
  delete _cellWriter;
  for(unsigned int i = 0; i < _cellQWriter.size(); i++) {
    delete _cellQWriter.at(i);
  }
  for(unsigned int i = 0; i < _cellParameterWithoutGhostlayerWriter.size(); i++) {
    delete _cellParameterWithoutGhostlayerWriter.at(i);
  }
  for(unsigned int i = 0; i < _cellParameterWithGhostlayerWriter.size(); i++) {
    delete _cellParameterWithGhostlayerWriter.at(i);
  }

  if(_plotMetainformation) {
    delete _cellSubdivisionFactorWriter;
    delete _cellGhostLayerWidthWriter;
    delete _cellTimeOldWriter;
    delete _cellTimeNewWriter;
    delete _cellDemandedMeshWidthWriter;
    delete _cellAgeWriter;
    #ifdef Parallel
    delete _cellRankWriter;
    #endif
  }
}

void peanoclaw::grid::plotter::VTKSubgridPlotter::plotSubcell(
  const Patch&                            patch,
  const peanoclaw::grid::SubgridAccessor& accessor,
  tarch::la::Vector<DIMENSIONS, int>      subcellIndex//,
//  peanoclaw::Vertex * const               vertices,
//  const peano::grid::VertexEnumerator&    enumerator
) {
  double localGap = _gap * patch.getLevel();
  tarch::la::Vector<DIMENSIONS, double> subcellSize = patch.getSubcellSize() / (1.0 + localGap);
  tarch::la::Vector<DIMENSIONS, double> x = patch.getPosition() + patch.getSize() * localGap/2.0 + tarch::la::multiplyComponents(subcellIndex.convertScalar<double>(), subcellSize);

  //Retreive adjacent vertices.
  int vertexIndices[TWO_POWER_D];
  dfor2(vertexIndex)
    tarch::la::Vector<DIMENSIONS,double> currentVertexPosition = x + tarch::la::multiplyComponents(vertexIndex.convertScalar<double>(), subcellSize);
    assertion1 ( _vertex2IndexMap.find(currentVertexPosition) != _vertex2IndexMap.end(), currentVertexPosition );
    vertexIndices[vertexIndexScalar] = _vertex2IndexMap[currentVertexPosition];
  enddforx

  int number = -1;
#ifdef Dim2
    number = _cellWriter->plotQuadrangle(vertexIndices);
#elif Dim3
    number = _cellWriter->plotHexahedron(vertexIndices);
#endif

    int qWriterIndex = 0;
    //for(int i = 0; i < (int)_cellQWriter.size(); i++) {
    for(int i = 0; i < patch.getUnknownsPerSubcell(); i++) {
      if(_plotQ.empty() || _plotQ.find(i) != _plotQ.end()) {
        _cellQWriter[qWriterIndex]->plotCell(number, accessor.getValueUNew(subcellIndex, i));
        qWriterIndex++;
      }
    }

    int parameterWithoutGhostlayerWriterIndex = 0;
    for(int i = 0; i < patch.getNumberOfParametersWithoutGhostlayerPerSubcell(); i++) {
      if(_plotParameterWithoutGhostlayer.empty() || _plotParameterWithoutGhostlayer.find(i) != _plotParameterWithoutGhostlayer.end()) {
        _cellParameterWithoutGhostlayerWriter[parameterWithoutGhostlayerWriterIndex]->plotCell(number, accessor.getParameterWithoutGhostlayer(subcellIndex, i));
        parameterWithoutGhostlayerWriterIndex++;
      }
    }

    int parameterWithGhostlayerWriterIndex = 0;
    for(int i = 0; i < patch.getNumberOfParametersWithGhostlayerPerSubcell(); i++) {
      if(_plotParameterWithGhostlayer.empty() || _plotParameterWithGhostlayer.find(i) != _plotParameterWithGhostlayer.end()) {
        _cellParameterWithGhostlayerWriter[parameterWithGhostlayerWriterIndex]->plotCell(number, accessor.getParameterWithGhostlayer(subcellIndex, i));
        parameterWithGhostlayerWriterIndex++;
      }
    }

    if(_plotMetainformation) {
      _cellSubdivisionFactorWriter->plotCell(number, patch.getSubdivisionFactor()(0));
      _cellGhostLayerWidthWriter->plotCell(number, patch.getGhostlayerWidth());
      _cellTimeOldWriter->plotCell(number, patch.getTimeIntervals().getCurrentTime());
      _cellTimeNewWriter->plotCell(number, patch.getTimeIntervals().getCurrentTime() + patch.getTimeIntervals().getTimestepSize());
      _cellDemandedMeshWidthWriter->plotCell(number, tarch::la::norm1(patch.getDemandedMeshWidth()) / DIMENSIONS);
      _cellAgeWriter->plotCell(number, patch.getAge());

      #ifdef Parallel
      _cellRankWriter->plotCell(number, tarch::parallel::Node::getInstance().getRank());
      #endif
    }
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::grid::plotter::VTKSubgridPlotter::computeGradient(
  const Patch&                            patch,
  const peanoclaw::grid::SubgridAccessor& accessor,
  tarch::la::Vector<DIMENSIONS, int>      subcellIndex,
  int                                     unknown
) {
  tarch::la::Vector<DIMENSIONS, double> gradient;

  for(int d = 0; d < DIMENSIONS; d++) {
    tarch::la::Vector<DIMENSIONS, int> neighborIndex = subcellIndex;
    neighborIndex(d)++;
    gradient(d) = (accessor.getValueUOld(neighborIndex, unknown) - accessor.getValueUOld(subcellIndex, unknown)) / patch.getSubcellSize()(d);
  }

  return gradient;
}

void peanoclaw::grid::plotter::VTKSubgridPlotter::plotSubgrid(
  const Patch& subgrid
) {
  double localGap = _gap * subgrid.getLevel();
  tarch::la::Vector<DIMENSIONS, double> subcellSize = subgrid.getSubcellSize() / (1.0 + localGap);

  peanoclaw::grid::SubgridAccessor accessor = subgrid.getAccessor();

  // Plot vertices
  dfor( vertexIndex, subgrid.getSubdivisionFactor()+1) {
    tarch::la::Vector<DIMENSIONS, double> x = subgrid.getPosition() + subgrid.getSize() * localGap/2.0 + tarch::la::multiplyComponents(vertexIndex.convertScalar<double>(), subcellSize);
    if ( _vertex2IndexMap.find(x) == _vertex2IndexMap.end() ) {
      assertion(_vertexWriter!=NULL);
      #if defined(Dim2) || defined(Dim3)
      int index = _vertexWriter->plotVertex(x);
      _vertex2IndexMap[x] = index;
      #else
      _vertex2IndexMap[x] = _vertexWriter->plotVertex(tarch::la::Vector<3,double>(x.data()));
      #endif
    }
  }

  // Plot cells
  dfor(subcellIndex, subgrid.getSubdivisionFactor()) {
//    if(tarch::la::norm2(computeGradient(patch, subcellIndex, 0)) > 1.0) {
      plotSubcell(
        subgrid,
        accessor,
        subcellIndex//,
//        vertices,
//        enumerator
      );
//    }
  }
}


void peanoclaw::grid::plotter::VTKSubgridPlotter::close() {
  _vertexWriter->close();
  _cellWriter->close();
  for(std::vector<tarch::plotter::griddata::Writer::CellDataWriter*>::iterator i = _cellQWriter.begin(); i != _cellQWriter.end(); i++) {
    (*i)->close();
  }
  for(std::vector<tarch::plotter::griddata::Writer::CellDataWriter*>::iterator i = _cellParameterWithoutGhostlayerWriter.begin(); i != _cellParameterWithoutGhostlayerWriter.end(); i++) {
    (*i)->close();
  }
  for(std::vector<tarch::plotter::griddata::Writer::CellDataWriter*>::iterator i = _cellParameterWithGhostlayerWriter.begin(); i != _cellParameterWithGhostlayerWriter.end(); i++) {
    (*i)->close();
  }
  if(_plotMetainformation) {
    _cellSubdivisionFactorWriter->close();
    _cellGhostLayerWidthWriter->close();
    _cellTimeOldWriter->close();
    _cellTimeNewWriter->close();
    _cellDemandedMeshWidthWriter->close();
    _cellAgeWriter->close();
    #ifdef Parallel
    _cellRankWriter->close();
    #endif
  }
}



