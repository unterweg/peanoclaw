/*
 * TimeIntervals.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: kristof
 */
#include "peanoclaw/grid/TimeIntervals.h"

#include "peanoclaw/Patch.h"

tarch::logging::Log peanoclaw::grid::TimeIntervals::_log("peanoclaw::grid::TimeIntervals");

peanoclaw::grid::TimeIntervals::TimeIntervals()
  : _cellDescription(0)
{
}

peanoclaw::grid::TimeIntervals::TimeIntervals(
  CellDescription* cellDescription
) : _cellDescription(cellDescription)
{
}

bool peanoclaw::grid::TimeIntervals::isAllowedToAdvanceInTime(
  double                                   maximumTimestepDueToGlobalTimestep,
  peanoclaw::Vertex * const                fineGridVertices,
  const peano::grid::VertexEnumerator&     fineGridVerticesEnumerator,
  peanoclaw::Vertex * const                coarseGridVertices,
  const peano::grid::VertexEnumerator&     coarseGridVerticesEnumerator
) {
  //Does Patch coarsen?
  bool patchCoarsening =
                          (peano::grid::aspects::VertexStateAnalysis::doesOneVertexCarryRefinementFlag
                            (
                              coarseGridVertices,
                              coarseGridVerticesEnumerator,
                              peanoclaw::records::Vertex::Erasing
                            )
                          ||
                          peano::grid::aspects::VertexStateAnalysis::doesOneVertexCarryRefinementFlag
                            (
                              coarseGridVertices,
                              coarseGridVerticesEnumerator,
                              peanoclaw::records::Vertex::EraseTriggered
                            )
                          )
                          &&
                          !peano::grid::aspects::VertexStateAnalysis::doesOneVertexCarryRefinementFlag
                            (
                              coarseGridVertices,
                              coarseGridVerticesEnumerator,
                              peanoclaw::records::Vertex::Refined
                            )
                          &&
                          !peano::grid::aspects::VertexStateAnalysis::doesOneVertexCarryRefinementFlag
                            (
                              coarseGridVertices,
                              coarseGridVerticesEnumerator,
                              peanoclaw::records::Vertex::Refining
                            )
                          &&
                          !peano::grid::aspects::VertexStateAnalysis::doesOneVertexCarryRefinementFlag
                            (
                              coarseGridVertices,
                              coarseGridVerticesEnumerator,
                              peanoclaw::records::Vertex::RefinementTriggered
                            );

  //Are all adjacent vertices set correctly?
  bool allAdjacentVerticesValid = true;
  for(int i = 0; i < TWO_POWER_D; i++) {
    peanoclaw::Vertex& vertex = fineGridVertices[fineGridVerticesEnumerator(i)];
    allAdjacentVerticesValid &= (vertex.getAgeInGridIterations() > 0) || vertex.isHangingNode();
  }

  //Statistics
//  if(!tarch::la::greater(maximumTimestepDueToGlobalTimestep, 0.0)) {
//    _subgridStatistics.addBlockedPatchDueToGlobalTimestep(patch);
//  } else if (!patch.getTimeIntervals().isAllowedToAdvanceInTime()) {
//    _subgridStatistics.addBlockedPatchDueToNeighborTimeConstraint(patch);
//  } else if (patch.shouldSkipNextGridIteration()) {
//    _subgridStatistics.addBlockedPatchDueToSkipIteration(patch);
//  } else if (patchCoarsening) {
//    _subgridStatistics.addBlockedPatchDueToCoarsening(patch);
//  }

  return
    tarch::la::greater(maximumTimestepDueToGlobalTimestep, 0.0)
    && !isBlockedByNeighbors()
    && !_cellDescription->getSkipGridIterations() > 0
    && !patchCoarsening
    && allAdjacentVerticesValid;
}

bool peanoclaw::grid::TimeIntervals::isBlockedByNeighbors() const {
  return tarch::la::greater(
        _cellDescription->getTime() + _cellDescription->getTimestepSize(),
        _cellDescription->getMinimalNeighborTimeConstraint());
}


double peanoclaw::grid::TimeIntervals::getCurrentTime() const {
  return _cellDescription->getTime();
}

double peanoclaw::grid::TimeIntervals::getTimeUOld() const {
  assertion1(Patch::isLeaf(_cellDescription) || Patch::isVirtual(_cellDescription), toString());
  if (Patch::isLeaf(_cellDescription)) {
    return _cellDescription->getTime();
  } else {
    return _cellDescription->getMinimalNeighborTime();
  }
}

double peanoclaw::grid::TimeIntervals::getTimeUNew() const {
  assertion1(Patch::isLeaf(_cellDescription) || Patch::isVirtual(_cellDescription), toString());
  if (Patch::isLeaf(_cellDescription)) {
    return _cellDescription->getTime() + _cellDescription->getTimestepSize();
  } else {
    return _cellDescription->getMinimalNeighborTime()
        + _cellDescription->getMaximalNeighborTimestep();
  }
}

void peanoclaw::grid::TimeIntervals::advanceInTime() {
  _cellDescription->setTime(
      _cellDescription->getTime() + _cellDescription->getTimestepSize());
}

double peanoclaw::grid::TimeIntervals::getTimestepSize() const {
  return _cellDescription->getTimestepSize();
}

void peanoclaw::grid::TimeIntervals::setTimestepSize(double timestepSize) {
  _cellDescription->setTimestepSize(timestepSize);
}

double peanoclaw::grid::TimeIntervals::getEstimatedNextTimestepSize() const {
  return _cellDescription->getEstimatedNextTimestepSize();
}

void peanoclaw::grid::TimeIntervals::setEstimatedNextTimestepSize(
    double estimatedNextTimestepSize) {
  _cellDescription->setEstimatedNextTimestepSize(estimatedNextTimestepSize);
}

double peanoclaw::grid::TimeIntervals::getMinimalNeighborTimeConstraint() const {
  return _cellDescription->getMinimalNeighborTimeConstraint();
}

double peanoclaw::grid::TimeIntervals::getMinimalLeafNeighborTimeConstraint() const {
  return _cellDescription->getMinimalLeafNeighborTimeConstraint();
}

void peanoclaw::grid::TimeIntervals::updateMinimalNeighborTimeConstraint(
    double neighborTimeConstraint, int neighborIndex) {
  if (tarch::la::smaller(neighborTimeConstraint,
      _cellDescription->getMinimalNeighborTimeConstraint())) {
    _cellDescription->setMinimalNeighborTimeConstraint(neighborTimeConstraint);
    _cellDescription->setConstrainingNeighborIndex(neighborIndex);
  }
}

void peanoclaw::grid::TimeIntervals::updateMinimalLeafNeighborTimeConstraint(
    double leafNeighborTime) {
  if (tarch::la::smaller(leafNeighborTime,
      _cellDescription->getMinimalLeafNeighborTimeConstraint())) {
    _cellDescription->setMinimalLeafNeighborTimeConstraint(leafNeighborTime);
  }
}

void peanoclaw::grid::TimeIntervals::resetMinimalNeighborTimeConstraint() {
  _cellDescription->setMinimalNeighborTimeConstraint(
      std::numeric_limits<double>::max());
  _cellDescription->setMinimalLeafNeighborTimeConstraint(
      std::numeric_limits<double>::max());
  _cellDescription->setConstrainingNeighborIndex(-1);
}

int peanoclaw::grid::TimeIntervals::getConstrainingNeighborIndex() const {
  return _cellDescription->getConstrainingNeighborIndex();
}

void peanoclaw::grid::TimeIntervals::resetMaximalNeighborTimeInterval() {
  _cellDescription->setMinimalNeighborTime(std::numeric_limits<double>::max());
  _cellDescription->setMaximalNeighborTimestep(-std::numeric_limits<double>::max());
}

void peanoclaw::grid::TimeIntervals::updateMaximalNeighborTimeInterval(double neighborTime,
    double neighborTimestepSize) {
  if (neighborTime + neighborTimestepSize
      > _cellDescription->getMinimalNeighborTime()
          + _cellDescription->getMaximalNeighborTimestep()) {
    _cellDescription->setMaximalNeighborTimestep(
        neighborTime + neighborTimestepSize
            - _cellDescription->getMinimalNeighborTime());
  }

  if (neighborTime < _cellDescription->getMinimalNeighborTime()) {
    _cellDescription->setMaximalNeighborTimestep(
        _cellDescription->getMaximalNeighborTimestep()
            + _cellDescription->getMinimalNeighborTime() - neighborTime);
    _cellDescription->setMinimalNeighborTime(neighborTime);
  }
}

void peanoclaw::grid::TimeIntervals::resetMinimalFineGridTimeInterval() {
  _cellDescription->setMaximumFineGridTime(-1.0);
  _cellDescription->setMinimumFineGridTimestep(
      std::numeric_limits<double>::max());
}

void peanoclaw::grid::TimeIntervals::updateMinimalFineGridTimeInterval(double fineGridTime,
    double fineGridTimestepSize) {
  if (fineGridTime > _cellDescription->getMaximumFineGridTime()) {
    _cellDescription->setMinimumFineGridTimestep(
        _cellDescription->getMinimumFineGridTimestep()
            - (fineGridTime - _cellDescription->getMaximumFineGridTime()));
    _cellDescription->setMaximumFineGridTime(fineGridTime);
  }
  _cellDescription->setMinimumFineGridTimestep(
      std::min(_cellDescription->getMinimumFineGridTimestep(),
          (fineGridTime + fineGridTimestepSize
              - _cellDescription->getMaximumFineGridTime())));
}

double peanoclaw::grid::TimeIntervals::getTimeConstraint() const {
  return _cellDescription->getTime() + _cellDescription->getTimestepSize();
}

void peanoclaw::grid::TimeIntervals::setFineGridsSynchronize(bool synchronizeFineGrids) {
  _cellDescription->setSynchronizeFineGrids(synchronizeFineGrids);
}

bool peanoclaw::grid::TimeIntervals::shouldFineGridsSynchronize() const {
  return _cellDescription->getSynchronizeFineGrids();
}

std::string peanoclaw::grid::TimeIntervals::toString() const {
  std::stringstream str;
  if (_cellDescription != 0) {
    str << "currentTime=" << _cellDescription->getTime() << ",timestepSize="
    << _cellDescription->getTimestepSize()
    << ",minimalNeighborTimeConstraint="
    << _cellDescription->getMinimalNeighborTimeConstraint();
  } else {
    str << "null";
  }
  return str.str();
}
