/*
 * NeighbourCommunicator.cpp
 *
 *  Created on: Mar 15, 2013
 *      Author: unterweg
 */
#include "NeighbourCommunicator.h"

#include "peanoclaw/Heap.h"
#include "peanoclaw/ParallelSubgrid.h"
#include "peanoclaw/records/CellDescription.h"
#include "peanoclaw/records/Data.h"

#include "peano/utils/Loop.h"

tarch::logging::Log peanoclaw::parallel::NeighbourCommunicator::_log("peanoclaw::parallel::NeighbourCommunicator");

void peanoclaw::parallel::NeighbourCommunicator::sendPatch(
  const Patch& transferedSubgrid
) {
  logTraceInWith3Arguments("sendPatch", cellDescriptionIndex, _position, _level);
  #ifdef Parallel

  #ifdef Asserts
  assertion3(tarch::la::allGreater(transferedSubgrid.getSubdivisionFactor(), 0), transferedSubgrid.toString(), _position, _level);
  assertion1(!transferedSubgrid.isRemote(), transferedSubgrid);

  //Check for zeros in transfered patch
  if(transferedSubgrid.isValid() && transferedSubgrid.isLeaf()) {
    dfor(subcellIndex, transferedSubgrid.getSubdivisionFactor()) {
      assertion3(tarch::la::greater(transferedSubgrid.getValueUNew(subcellIndex, 0), 0.0), subcellIndex, transferedSubgrid, transferedSubgrid.toStringUNew());
      assertion3(tarch::la::greater(transferedSubgrid.getValueUOld(subcellIndex, 0), 0.0), subcellIndex, transferedSubgrid, transferedSubgrid.toStringUOldWithGhostLayer());
    }
  }
  #endif

  _statistics.sentNeighborData();
  _subgridCommunicator.sendCellDescription(transferedSubgrid.getCellDescriptionIndex());

  if(transferedSubgrid.getUNewIndex() != -1) {
    _subgridCommunicator.sendDataArray(transferedSubgrid.getUNewIndex());
  } else {
    _subgridCommunicator.sendPaddingDataArray();
  }
  #endif
  logTraceOut("sendPatch");
}

void peanoclaw::parallel::NeighbourCommunicator::sendPaddingPatch(
  const tarch::la::Vector<DIMENSIONS, double>& position,
  int                                          level,
  const tarch::la::Vector<DIMENSIONS, double>& subgridSize
) {
  logTraceInWith2Arguments("sendPaddingPatch", _position, _level);

  //TODO unterweg debug
  logDebug("", "Sending padding patch to " << _remoteRank << " at " << position << " on level " << level);

  _statistics.sentPaddingNeighborData();
  _subgridCommunicator.sendPaddingCellDescription(
    position,
    level,
    subgridSize
  );
  _subgridCommunicator.sendPaddingDataArray(); //UNew
  logTraceOut("sendPaddingPatch");
}

void peanoclaw::parallel::NeighbourCommunicator::receivePatch(Patch& localSubgrid) {
  #ifdef Parallel
  logTraceInWith3Arguments("receivePatch", localCellDescriptionIndex, _position, _level);

  std::vector<CellDescription> remoteCellDescriptionVector = CellDescriptionHeap::getInstance().receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);

  logDebug("", "Receiving patch from " << _remoteRank << " at " << localCellDescription.getPosition() << " on level " << localCellDescription.getLevel());

  if(!_onlySendSubgridsAfterChange || remoteCellDescriptionVector.size() > 0) {
    assertionEquals4(remoteCellDescriptionVector.size(), 1, _position, _level, _remoteRank, localSubgrid);
    CellDescription remoteCellDescription = remoteCellDescriptionVector[0];
    assertionEquals2(localSubgrid.getLevel(), _level, localSubgrid, tarch::parallel::Node::getInstance().getRank());

    #ifdef Asserts
    logDebug("", "Receiving patch: " << remoteCellDescription.toString());
    assertion6(!remoteCellDescription.getIsRemote(), localSubgrid, remoteCellDescription.toString(), _position, _level, _remoteRank, tarch::parallel::Node::getInstance().getRank());
    assertionNumericalEquals6(localSubgrid.getPosition(), remoteCellDescription.getPosition(),
        localSubgrid, remoteCellDescription.toString(), tarch::parallel::Node::getInstance().getRank(), _remoteRank, _position, _level);
    assertionNumericalEquals4(localSubgrid.getSize(), remoteCellDescription.getSize(),
        localSubgrid, remoteCellDescription.toString(), tarch::parallel::Node::getInstance().getRank(), _remoteRank);
    assertionEquals4(localSubgrid.getLevel(), remoteCellDescription.getLevel(),
        localSubgrid, remoteCellDescription.toString(), tarch::parallel::Node::getInstance().getRank(), _remoteRank);
    #endif

    _statistics.receivedNeighborData();
    if(remoteCellDescription.getUIndex() != -1) {
      remoteCellDescription.setUIndex(_subgridCommunicator.receiveDataArray());
    } else {
      _subgridCommunicator.receiveDataArray();
    }

    //Copy remote cell description to local cell description
    _subgridCommunicator.deleteArraysFromSubgrid(localSubgrid);
    remoteCellDescription.setCellDescriptionIndex(localSubgrid.getCellDescriptionIndex());
    remoteCellDescription.setIsRemote(true); //TODO unterweg: Remote patches are currently never destroyed.
    CellDescriptionHeap::getInstance().getData(localSubgrid.getCellDescriptionIndex()).at(0) = remoteCellDescription;
    assertionEquals(CellDescriptionHeap::getInstance().getData(localSubgrid.getCellDescriptionIndex()).size(), 1);

    //Initialize non-parallel fields
    Patch remotePatch(remoteCellDescription);
    remotePatch.initializeNonParallelFields();

    //Check for zeros in transfered patch
    #ifdef Asserts
    if(remotePatch.isValid() && remotePatch.isLeaf()) {
      dfor(subcellIndex, remotePatch.getSubdivisionFactor()) {
        assertion3(tarch::la::greater(remotePatch.getValueUNew(subcellIndex, 0), 0.0), subcellIndex, remotePatch, remotePatch.toStringUNew());
        assertion3(tarch::la::greater(remotePatch.getValueUOld(subcellIndex, 0), 0.0), subcellIndex, remotePatch, remotePatch.toStringUOldWithGhostLayer());
      }
    }
    assertionEquals(CellDescriptionHeap::getInstance().getData(localSubgrid.getCellDescriptionIndex()).at(0).getCellDescriptionIndex(), localSubgrid.getCellDescriptionIndex());
    assertion1(!tarch::la::smaller(remotePatch.getTimeIntervals().getTimestepSize(), 0.0) || !remotePatch.isLeaf(), remotePatch);
    #endif
  } else {
    //Padding patch received -> receive padding data
    DataHeap::getInstance().receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  }

  logTraceOut("receivePatch");
  #endif
}

void peanoclaw::parallel::NeighbourCommunicator::receivePaddingPatch() {
  logTraceInWith2Arguments("receivePaddingPatch", _position, _level);
  #ifdef Parallel
  _statistics.receivedPaddingNeighborData();

  logDebug("", "Receiving padding patch from " << _remoteRank << " at " << _position << " on level " << _level);

  //Receive padding patch
  std::vector<CellDescription> remoteCellDescriptionVector = CellDescriptionHeap::getInstance().receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  assertionEquals4(remoteCellDescriptionVector.size(), 0, _position, _level, _remoteRank, remoteCellDescriptionVector[0].toString());

  //UNew
  DataHeap::getInstance().receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  #endif
  logTraceOut("receivePaddingPatch");
}

void peanoclaw::parallel::NeighbourCommunicator::createOrFindRemoteSubgrid(
  Vertex& localVertex,
  int     adjacentSubgridIndexInPeanoOrder,
  const tarch::la::Vector<DIMENSIONS, double>& subgridSize
) {
  tarch::la::Vector<DIMENSIONS, double> subgridPosition
    = _position + tarch::la::multiplyComponents(subgridSize, peano::utils::dDelinearised(adjacentSubgridIndexInPeanoOrder, 2).convertScalar<double>() - 1.0);
  tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> key = createRemoteSubgridKey(subgridPosition, _level);
   if(_remoteSubgridMap.find(key) != _remoteSubgridMap.end()) {
     localVertex.setAdjacentCellDescriptionIndexInPeanoOrder(adjacentSubgridIndexInPeanoOrder, _remoteSubgridMap[key]);
   } else {
    Patch outsidePatch(
      subgridPosition,
      subgridSize,
      0,
      0,
      1,
      1,
      1.0,
      _level
    );
    localVertex.setAdjacentCellDescriptionIndexInPeanoOrder(adjacentSubgridIndexInPeanoOrder, outsidePatch.getCellDescriptionIndex());
    _remoteSubgridMap[key] = outsidePatch.getCellDescriptionIndex();
   }
}

tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> peanoclaw::parallel::NeighbourCommunicator::createRemoteSubgridKey(
  const tarch::la::Vector<DIMENSIONS, double> subgridPosition,
  int                                         level
) const {
  tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> key;
  for(int d = 0; d < DIMENSIONS; d++) {
    key(d) = subgridPosition(d);
  }
  key(DIMENSIONS) = level;
  return key;
}

peanoclaw::parallel::NeighbourCommunicator::NeighbourCommunicator(
  int                                         remoteRank,
  const tarch::la::Vector<DIMENSIONS,double>& position,
  int                                         level,
  const tarch::la::Vector<DIMENSIONS,double>& subgridSize,
  RemoteSubgridMap&                           remoteSubgridMap,
  peanoclaw::statistics::ParallelStatistics&  statistics
) : _subgridCommunicator(remoteRank, position, level, peano::heap::NeighbourCommunication),
    _remoteRank(remoteRank),
    _position(position),
    _level(level),
    _subgridSize(subgridSize),
    _remoteSubgridMap(remoteSubgridMap),
    _statistics(statistics),
    //En-/Disable optimizations
    _avoidMultipleTransferOfSubgridsIfPossible(true),
    _onlySendSubgridsAfterChange(true) {
  logTraceInWith3Arguments("NeighbourCommunicator", remoteRank, position, level);

  logTraceOut("NeighbourCommunicator");
}

void peanoclaw::parallel::NeighbourCommunicator::sendSubgridsForVertex(
  peanoclaw::Vertex&                           vertex,
  const tarch::la::Vector<DIMENSIONS, double>& vertexPosition,
  const tarch::la::Vector<DIMENSIONS, double>& adjacentSubgridSize,
  int                                          level,
  const State&                                 state
) {
  #ifdef Parallel

  //TODO unterweg debug
  logDebug("sendSubgridsForVertex", "Sending subgrids to rank " << _remoteRank
      << " for vertex " << vertexPosition << " on level " << level << " adj:" << vertex.getAdjacentRanks()
      << " joining:" << state.isJoiningWithMaster() << " " << state.isJoinTriggered() << " " << state.isJoinWithMasterTriggered()
      << " " << state.isJoiningWithWorker() << " " << state.isInvolvedInJoinOrFork()
      << " master:" << (tarch::parallel::Node::getInstance().isGlobalMaster() ? -1 : tarch::parallel::NodePool::getInstance().getMasterRank()));

  if(!tarch::parallel::Node::getInstance().isGlobalMaster() && _remoteRank != 0) {
    int localRank = tarch::parallel::Node::getInstance().getRank();

    for(int i = 0; i < TWO_POWER_D; i++) {
      int localSubgridRank = vertex.getAdjacentRanks()(i);
      if(state.isJoinWithMasterTriggered() && localSubgridRank == localRank) {
        localSubgridRank = tarch::parallel::NodePool::getInstance().getMasterRank();
      }

      bool sentSubgrid = false;
      if(localSubgridRank == localRank && vertex.getAdjacentCellDescriptionIndexInPeanoOrder(i) != -1) {
        Patch           localSubgrid(vertex.getAdjacentCellDescriptionIndexInPeanoOrder(i));
        ParallelSubgrid localParallelSubgrid(vertex.getAdjacentCellDescriptionIndexInPeanoOrder(i));

        if(
            !localParallelSubgrid.wasCurrentStateSent()
            || vertex.wereAdjacentRanksChanged()
            || vertex.getAgeInGridIterations() <= 1
            //TODO unterweg dissertion
            //Warum kann ich nicht einfach den Status der Feingitter restringieren, ob die sich geändert haben?
            //Gibt es da einen Zusammenhang zwischen Zeitschritten der benachbarten Grobgitter und den Zeitschritten
            //der Feingitter, den ich hier nicht erkenne?
            //TODO unterweg debug
            || !localSubgrid.isLeaf()
            || !_onlySendSubgridsAfterChange) {

          //Find adjacent rank
          bool isSingleTransferOfSubgridToRank = false;
          tarch::la::Vector<THREE_POWER_D_MINUS_ONE, int> adjacentRanks = localParallelSubgrid.getAdjacentRanks();
          for(int i = 0; i < THREE_POWER_D_DIVIDED_BY_THREE; i++) {
            if(adjacentRanks(i) == _remoteRank) {
              if(localParallelSubgrid.getNumberOfSharedAdjacentVertices()(i) == 1) {
                isSingleTransferOfSubgridToRank = true;
              }
              break;
            }
          }

//          logInfo("sendSubgridsForVertex", "Sending subgrid to rank " << _remoteRank << " at " << localSubgrid.getPosition() << " on level " << localSubgrid.getLevel()
//              << " adjacentRanks: " << localParallelSubgrid.getAdjacentRanks() << " shared: " << localParallelSubgrid.getNumberOfSharedAdjacentVertices()
//              << " isSingleTransfer: " << isSingleTransferOfSubgridToRank << " transferToSkip: " << localParallelSubgrid.getNumberOfTransfersToBeSkipped());

          if(
            isSingleTransferOfSubgridToRank
            || !_avoidMultipleTransferOfSubgridsIfPossible
          ) {
            logDebug("sendSubgridsForVertex", "Sending subgrid to rank " << _remoteRank << ": " << localSubgrid << " for vertex " << _position);

            sendPatch(localSubgrid);

            sentSubgrid = true;
          } else {
            logDebug("sendSubgridsForVertex", "(Skipped) Sending subgrid to rank " << _remoteRank << ": " << localSubgrid
                << " for vertex " << _position << ", getNumberOfSharedAdjacentVertices=" << localParallelSubgrid.getNumberOfSharedAdjacentVertices());
          }
          localParallelSubgrid.decreaseNumberOfSharedAdjacentVertices(_remoteRank);
        }
      }

      if(!sentSubgrid) {
        logDebug("sendSubgridsForVertex", "Sending padding subgrid to rank " << _remoteRank
            << " at " << (_position + _subgridSize * (peano::utils::dDelinearised(i, 2).convertScalar<double>()*2.0-1.0)) << " on level " << _level);

        sendPaddingPatch(_position, _level, _subgridSize);
      }
    }
  }
  #endif
}

void peanoclaw::parallel::NeighbourCommunicator::receiveSubgridsForVertex(
  peanoclaw::Vertex&                           localVertex,
  const peanoclaw::Vertex&                     remoteVertex,
  const tarch::la::Vector<DIMENSIONS, double>& vertexPosition,
  const tarch::la::Vector<DIMENSIONS, double>& adjacentSubgridSize,
  int                                          level
) {
  #ifdef Parallel

  //TODO unterweg debug
  logDebug("receiveSubgridsForVertex", "Receiving subgrids from rank " << _remoteRank
      << " for vertex " << vertexPosition << " on level " << level << " adj(local):" << localVertex.getAdjacentRanks() << " adj(remote):" << remoteVertex.getAdjacentRanks());

  if(!tarch::parallel::Node::getInstance().isGlobalMaster() && _remoteRank != 0) {
    for(int i = TWO_POWER_D - 1; i >= 0; i--) {
      int localAdjacentCellDescriptionIndex = localVertex.getAdjacentCellDescriptionIndexInPeanoOrder(i);
      int remoteAdjacentCellDescriptionIndex = remoteVertex.getAdjacentCellDescriptionIndexInPeanoOrder(i);

      bool receivedSubgrid = false;
      if(remoteVertex.getAdjacentRanks()(i) == _remoteRank && remoteAdjacentCellDescriptionIndex != -1) {

          //Create remote subgrid if necessary
          if(localVertex.getAdjacentCellDescriptionIndexInPeanoOrder(i) == -1 && localVertex.getAdjacentRanks()(i) != 0) {
            createOrFindRemoteSubgrid(localVertex, i, adjacentSubgridSize);
            localAdjacentCellDescriptionIndex = localVertex.getAdjacentCellDescriptionIndexInPeanoOrder(i);
          }

          assertion(localAdjacentCellDescriptionIndex != -1);
          Patch           localSubgrid(localAdjacentCellDescriptionIndex);
          ParallelSubgrid localParallelSubgrid(localAdjacentCellDescriptionIndex);

//          logInfo("receiveSubgridsForVertex", "Receiving subgrid from rank " << _remoteRank << " at " << localSubgrid.getPosition() << " on level " << localSubgrid.getLevel()
//              << " transfersToSkip: " << localParallelSubgrid.getNumberOfTransfersToBeSkipped() << " adjacentRanks: " << localParallelSubgrid.getAdjacentRanks()
//              << " transfersToSkip: " << localParallelSubgrid.getNumberOfTransfersToBeSkipped()
//              );

          //Only receive if the incoming transfer was not skipped to avoid multiple send of the subgrid.
          if(localParallelSubgrid.getNumberOfTransfersToBeSkipped() == 0 || !_avoidMultipleTransferOfSubgridsIfPossible) {

            receivePatch(localSubgrid);
            receivedSubgrid = true;

            localSubgrid.reloadCellDescription();
            logDebug("receiveSubgridsForVertex", "Received subgrid from rank " << _remoteRank << ": " << localSubgrid
                << " for vertex " << _position);
          } else {

            assertion1(localSubgrid.isRemote(), localSubgrid);

            localParallelSubgrid.decreaseNumberOfTransfersToBeSkipped();

            logDebug("receiveSubgridsForVertex", "(Skipped) Received subgrid from rank " << _remoteRank << ": " << localSubgrid
                << " for vertex " << _position << ",getNumberOfTransfersToBeSkipped=" << localParallelSubgrid.getNumberOfTransfersToBeSkipped());
          }
      }


      if(!receivedSubgrid) {
        //Receive dummy message
        receivePaddingPatch();

        logDebug("receiveSubgridsForVertex", "Received padding subgrid from rank " << _remoteRank
            << " at "  << (_position + _subgridSize * (peano::utils::dDelinearised(i, 2).convertScalar<double>()*2.0-1.0)) << " on level " << _level);
      }
    }
  }
  #endif
}

void peanoclaw::parallel::NeighbourCommunicator::switchToRemote(peanoclaw::Patch& subgrid) {
  #ifdef Parallel
  subgrid.setIsRemote(true);
  tarch::la::Vector<DIMENSIONS_PLUS_ONE, double> key = createRemoteSubgridKey(subgrid.getPosition(), subgrid.getLevel());

  _remoteSubgridMap[key] = subgrid.getCellDescriptionIndex();
  #endif
}
