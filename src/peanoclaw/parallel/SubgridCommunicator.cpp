/*
 * SubgridCommunicator.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: kristof
 */
#include "peanoclaw/parallel/SubgridCommunicator.h"

#include "peanoclaw/geometry/Region.h"
#include "peanoclaw/ParallelSubgrid.h"

#include "tarch/Assertions.h"

#if defined(Parallel) && defined(UseBlockedMeshCommunication)
#include "peano/parallel/Serialization.h"
#include "peano/parallel/SerializationMap.h"
#endif

tarch::logging::Log peanoclaw::parallel::SubgridCommunicator::_log("peanoclaw::parallel::SubgridCommunicator");

std::vector<peanoclaw::records::CellDescription> peanoclaw::parallel::SubgridCommunicator::receiveCellDescription() {
  logTraceIn("receiveCellDescription()");

  std::vector<peanoclaw::records::CellDescription> remoteCellDescriptionVector;
  if(_packCommunication && _messageType == peano::heap::NeighbourCommunication) {
    #if defined(Parallel) && defined(UseBlockedMeshCommunication)
    Serialization::ReceiveBuffer& recvbuffer = peano::parallel::SerializationMap::getInstance().getReceiveBuffer(_remoteRank)[1];
    assertion1(recvbuffer.isBlockAvailable(), "cannot read heap data from Serialization Buffer - not enough blocks");

    Serialization::Block block = recvbuffer.nextBlock();

    assertion4(block.size() >= 1, "cell description block is too small",
        block.size(), tarch::parallel::Node::getInstance().getRank(), _remoteRank);

    unsigned char BlockType;
    block >> BlockType;

    assertion2(BlockType == 0x10, "expected cell description but got something else", BlockType);

    int cellDescriptionSize;
    MPI_Pack_size(1, CellDescription::Packed::Datatype, MPI_COMM_WORLD, &cellDescriptionSize);

    size_t numberOfCellDescriptions = (block.size() - 1) / cellDescriptionSize;
    assertion1(block.size() - 1 == 0 || block.size() - 1 == (size_t)cellDescriptionSize, block.size());

    //assertion1(numberOfCellDescriptions > 0, "no cell descriptions, huh? we always send one");

    remoteCellDescriptionVector.resize(numberOfCellDescriptions);

    int block_position = 1;
    //std::cout << " ||||||| unpacking " << numberOfCellDescriptions << " cell descriptions " << std::endl;
    for (size_t i=0; i < numberOfCellDescriptions; i++) {
        CellDescription::Packed packed;
        MPI_Unpack(block.data(), block.size(), &block_position, &packed, 1, CellDescription::Packed::Datatype, MPI_COMM_WORLD );
        remoteCellDescriptionVector[i] = packed.convert();
    }
    #endif
  } else {
    remoteCellDescriptionVector = CellDescriptionHeap::getInstance().receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  }

  logTraceOut("receiveCellDescription()");
  return remoteCellDescriptionVector;
}

peanoclaw::parallel::SubgridCommunicator::SubgridCommunicator(
  int                                         remoteRank,
  const tarch::la::Vector<DIMENSIONS,double>& position,
  int                                         level,
  peano::heap::MessageType                    messageType,
  bool                                        onlySendOverlappedCells,
  bool                                        packCommunication
) : _remoteRank(remoteRank),
    _position(position),
    _level(level),
    _messageType(messageType),
    _onlySendOverlappedCells(onlySendOverlappedCells),
    _packCommunication(packCommunication) {

}

void peanoclaw::parallel::SubgridCommunicator::sendSubgrid(Patch& subgrid) {
  logTraceInWith1Argument("sendSubgrid(Subgrid)", subgrid);

  if (!_packCommunication) {
      sendCellDescription(subgrid.getCellDescriptionIndex());
  }

  if(subgrid.getUIndex() != -1) {
    #if defined(Asserts) && defined(AssertForPositiveValues)
    if(subgrid.isValid() && subgrid.isLeaf() && !_onlySendOverlappedCells) {
      dfor(subcellIndex, subgrid.getSubdivisionFactor()) {
        assertion4(tarch::la::greater(subgrid.getAccessor().getValueUNew(subcellIndex, 0), 0.0), tarch::parallel::Node::getInstance().getRank(), subcellIndex, subgrid, subgrid.toStringUNew());
        assertion4(tarch::la::greater(subgrid.getAccessor().getValueUOld(subcellIndex, 0), 0.0), tarch::parallel::Node::getInstance().getRank(), subcellIndex, subgrid, subgrid.toStringUOldWithGhostLayer());
      }
    }
    assertionEquals(CellDescriptionHeap::getInstance().getData(subgrid.getCellDescriptionIndex()).at(0).getCellDescriptionIndex(), subgrid.getCellDescriptionIndex());
    assertion1(!tarch::la::smaller(subgrid.getTimeIntervals().getTimestepSize(), 0.0) || !subgrid.isLeaf(), subgrid);
    #endif

    if(_onlySendOverlappedCells) {
      sendOverlappedCells(subgrid);
    } else {
      sendDataArray(subgrid.getUIndex());
    }
  } else if(_messageType == peano::heap::NeighbourCommunication) {
    sendPaddingDataArray();
  } else {
    //Don't send anything
  }
 
  if (_packCommunication) {
    sendCellDescription(subgrid.getCellDescriptionIndex());
  }
  logTraceOut("sendSubgrid(Subgrid)");
}

void peanoclaw::parallel::SubgridCommunicator::sendPaddingSubgrid() {
  logTraceIn("sendPaddingSubgrid()");

  if (_packCommunication && _messageType == peano::heap::NeighbourCommunication) {
    sendPaddingDataArray();
    sendPaddingCellDescription();
  } else {
    sendPaddingCellDescription();
    sendPaddingDataArray();
  }

  logTraceOut("sendPaddingSubgrid()");
}

void peanoclaw::parallel::SubgridCommunicator::sendCellDescription(int cellDescriptionIndex)
{
  logTraceInWith1Argument("sendCellDescription", cellDescriptionIndex);
  #if defined(Asserts) && defined(Parallel)
  CellDescription& cellDescription = CellDescriptionHeap::getInstance().getData(cellDescriptionIndex).at(0);
  assertion1(!cellDescription.getIsPaddingSubgrid(), cellDescription.toString());
  #endif

  if(_packCommunication && _messageType == peano::heap::NeighbourCommunication) {
    #if defined(Parallel) && defined(UseBlockedMeshCommunication)
    Serialization::SendBuffer& sendbuffer = peano::parallel::SerializationMap::getInstance().getSendBuffer(_remoteRank)[1];

    std::vector<CellDescription>& localCellDescriptionVector = CellDescriptionHeap::getInstance().getData(cellDescriptionIndex);

    size_t numberOfCellDescriptions = localCellDescriptionVector.size();
    assertionEquals(numberOfCellDescriptions, 1);
    int cellDescriptionSize;// = sizeof(CellDescription::Packed);
    MPI_Pack_size(1, CellDescription::Packed::Datatype, MPI_COMM_WORLD, &cellDescriptionSize);
    Serialization::Block block = sendbuffer.reserveBlock(1+cellDescriptionSize*numberOfCellDescriptions);
    unsigned char CellDescriptionType = 0x10;
    block << CellDescriptionType;

    int block_position = 1;
    //std::cout << " ||||||| packing " << numberOfCellDescriptions << " cell descriptions " << std::endl;
    for (size_t i=0; i < numberOfCellDescriptions; i++) {
        CellDescription::Packed packed = localCellDescriptionVector[i].convert();
        MPI_Pack(&packed, 1, CellDescription::Packed::Datatype, block.data(), block.size(), &block_position, MPI_COMM_WORLD );
        assertionEquals((size_t)block_position, (1+cellDescriptionSize*numberOfCellDescriptions));
    }
    #endif
  } else {
    CellDescriptionHeap::getInstance().sendData(cellDescriptionIndex, _remoteRank, _position, _level, _messageType);
  }
  logTraceOut("sendCellDescription");
}

void peanoclaw::parallel::SubgridCommunicator::sendPaddingCellDescription() {
  logTraceIn("sendPaddingCellDescription");
  int cellDescriptionIndex = CellDescriptionHeap::getInstance().createData();

  if(_packCommunication && _messageType == peano::heap::NeighbourCommunication) {
    #if defined(Parallel) && defined(UseBlockedMeshCommunication)
    Serialization::SendBuffer& sendbuffer = peano::parallel::SerializationMap::getInstance().getSendBuffer(_remoteRank)[1];
    Serialization::Block block = sendbuffer.reserveBlock(1);
    unsigned char CellDescriptionType = 0x10;
    block << CellDescriptionType;
    #endif
  } else {
    CellDescriptionHeap::getInstance().sendData(cellDescriptionIndex, _remoteRank, _position, _level, _messageType);
  }

  CellDescriptionHeap::getInstance().deleteData(cellDescriptionIndex);
  logTraceOut("sendPaddingCellDescription");
}

void peanoclaw::parallel::SubgridCommunicator::sendDataArray(int index) {
  logTraceInWith3Arguments("sendDataArray", index, _position, _level);
 
  if ( _packCommunication && _messageType == peano::heap::NeighbourCommunication ) {
    #if defined(Parallel) && defined(UseBlockedMeshCommunication)
    Serialization::SendBuffer& sendbuffer = peano::parallel::SerializationMap::getInstance().getSendBuffer(_remoteRank)[1];

    std::vector<Data>& localDataVector = DataHeap::getInstance().getData(index);

    //assertion1(localDataVector.size() <= 27*27*3+29*29+1, "tooo large");

    size_t numberOfDataElements = localDataVector.size();
    int dataSize = sizeof(Data::Packed);
    Serialization::Block block = sendbuffer.reserveBlock(dataSize*numberOfDataElements);

    int block_position = 0;
    //std::cout << " ||||||| packing " << numberOfDataElements << " data elements" << std::endl;
    for (size_t i=0; i < numberOfDataElements; i++) {
        Data::Packed packed = localDataVector[i].convert();
        MPI_Pack(&packed, 1, Data::Packed::Datatype, block.data(), block.size(), &block_position, MPI_COMM_WORLD );
    }
    #endif
  } else {
 
      DataHeap::getInstance().sendData(index, _remoteRank, _position, _level, _messageType);

  }

  logTraceOut("sendDataArray");
}

void peanoclaw::parallel::SubgridCommunicator::sendOverlappedCells(
  Patch& subgrid
) {
  logTraceInWith1Argument("sendOverlappedCellsOfDataArray(...)", subgrid);

  ParallelSubgrid parallelSubgrid(subgrid);
  peanoclaw::grid::SubgridAccessor& subgridAccessor = subgrid.getAccessor();

  peanoclaw::geometry::Region regions[THREE_POWER_D_MINUS_ONE];
  int numberOfRegions = peanoclaw::geometry::Region::getRegionsOverlappedByRemoteGhostlayers(
    parallelSubgrid.getAdjacentRanks(),
    parallelSubgrid.getOverlapOfRemoteGhostlayers(),
    subgrid.getSubdivisionFactor(),
    _remoteRank,
    regions
  );

  int numberOfCells = 0;
  for(int i = 0; i < numberOfRegions; i++) {
    numberOfCells += tarch::la::volume(regions[i]._size);
  }

  //Allocate data array
  int numberOfEntries = numberOfCells * subgrid.getUnknownsPerSubcell() * 2;
  int temporaryIndex = DataHeap::getInstance().createData(numberOfEntries, numberOfEntries);
  std::vector<Data>& temporaryArray = DataHeap::getInstance().getData(temporaryIndex);

  int entry = 0;
  for(int i = 0; i < numberOfRegions; i++) {
    peanoclaw::geometry::Region& region = regions[i];

    //U new
    dfor(subcellIndex, region._size) {
      int linearIndex = subgridAccessor.getLinearIndexUNew(region._offset + subcellIndex);
      for(int unknown = 0; unknown < subgrid.getUnknownsPerSubcell(); unknown++) {
        temporaryArray[entry++] = subgridAccessor.getValueUNew(linearIndex, unknown);
      }

      #if defined(AssertForPositiveValues) && defined(Asserts)
      if(subgrid.isLeaf()) {
        assertion6(
          tarch::la::greater(subgridAccessor.getValueUNew(linearIndex, 0), 0.0),
          subgrid,
          subcellIndex,
          region._offset,
          region._size,
          subgridAccessor.getValueUNew(linearIndex, 0),
          subgrid.toStringUNew()
        );
      }
      #endif
    }

    //U old
    dfor(subcellIndex, region._size) {
      int linearIndex = subgridAccessor.getLinearIndexUOld(region._offset + subcellIndex);
      for(int unknown = 0; unknown < subgrid.getUnknownsPerSubcell(); unknown++) {
        temporaryArray[entry++] = subgridAccessor.getValueUOld(linearIndex, unknown);
      }
    }
  }

  //Send
  sendDataArray(temporaryIndex);

  assertion1(DataHeap::getInstance().isValidIndex(temporaryIndex), temporaryIndex);

  DataHeap::getInstance().deleteData(temporaryIndex);

  logTraceOut("sendOverlappedCellsOfDataArray(...)");
}

void peanoclaw::parallel::SubgridCommunicator::sendPaddingDataArray() {
  logTraceInWith2Arguments("sendPaddingDataArray", _position, _level);
  int index = DataHeap::getInstance().createData();
  sendDataArray(index);
  DataHeap::getInstance().deleteData(index);

  logTraceOut("sendPaddingDataArray");
}

void peanoclaw::parallel::SubgridCommunicator::receivePaddingSubgrid() {
  logTraceIn("receivePaddingSubgrid()");
  #ifdef Parallel
  logDebug("", "Receiving padding patch from " << _remoteRank << " at " << _position << " on level " << _level);

  std::vector<CellDescription> remoteCellDescriptionVector = receiveCellDescription();

  assertionEquals4(remoteCellDescriptionVector.size(), 0, _position, _level, _remoteRank, remoteCellDescriptionVector[0].toString());

  //UNew
  receivePaddingDataArray();

  #endif
  logTraceOut("receivePaddingSubgrid()");
}

int peanoclaw::parallel::SubgridCommunicator::receiveDataArray() {
  logTraceIn("receiveDataArray");

  int localIndex = DataHeap::getInstance().createData();
  
  if ( _packCommunication && _messageType == peano::heap::NeighbourCommunication ) {
    #if defined(Parallel) && defined(UseBlockedMeshCommunication)
    Serialization::ReceiveBuffer& recvbuffer = peano::parallel::SerializationMap::getInstance().getReceiveBuffer(_remoteRank)[1];
    assertion1(recvbuffer.isBlockAvailable(), "cannot read heap data from Serialization Buffer - not enough blocks");

    Serialization::Block block = recvbuffer.nextBlock();

    size_t numberOfDataElements = block.size() / sizeof(Data::Packed);

    //assertion1(numberOfCellDescriptions > 0, "no cell descriptions, huh? we always send one");
 
    std::vector<Data>& remoteDataVector = DataHeap::getInstance().getData(localIndex);
    remoteDataVector.resize(numberOfDataElements);

    int block_position = 0;
    //std::cout << " ||||||| unpacking " << numberOfCellDescriptions << " cell descriptions " << std::endl;
    for (size_t i=0; i < numberOfDataElements; i++) {
        Data::Packed packed;
        MPI_Unpack(block.data(), block.size(), &block_position, &packed, 1, Data::Packed::Datatype, MPI_COMM_WORLD );
        remoteDataVector[i] = packed.convert();
    }
    #endif
  } else {

      DataHeap::getInstance().receiveData(
        localIndex,
        _remoteRank,
        _position,
        _level,
        _messageType
      );

  }

  logTraceOut("receiveDataArray");
  return localIndex;
}

void peanoclaw::parallel::SubgridCommunicator::receivePaddingDataArray() {
  int paddingArrayIndex = receiveDataArray();
  DataHeap::getInstance().deleteData(paddingArrayIndex);
}

void peanoclaw::parallel::SubgridCommunicator::receiveOverlappedCells(
  const CellDescription& remoteCellDescription,
  Patch&                 subgrid
) {
  logTraceInWith1Argument("receiveOverlappedCells", subgrid);
  #ifdef Parallel
  peanoclaw::grid::SubgridAccessor& subgridAccessor = subgrid.getAccessor();
  peanoclaw::geometry::Region regions[THREE_POWER_D_MINUS_ONE];
  int numberOfRegions = peanoclaw::geometry::Region::getRegionsOverlappedByRemoteGhostlayers(
    remoteCellDescription.getAdjacentRanks(),
    remoteCellDescription.getOverlapByRemoteGhostlayer(),
    remoteCellDescription.getSubdivisionFactor(),
    tarch::parallel::Node::getInstance().getRank(),
    regions
  );

  //TODO unterweg debug
//  std::cout << "Receiving " << numberOfRegions << " regions from " << _remoteRank << " on " << tarch::parallel::Node::getInstance().getRank()
//      << " for subgrid " << subgrid.getPosition() << ", " << subgrid.getSize()
//      << ", adj:" << remoteCellDescription.getAdjacentRanks() << ", overlap:" << remoteCellDescription.getOverlapByRemoteGhostlayer() << ": " << std::endl;
//  for(int i = 0; i < numberOfRegions; i++) {
//    std::cout << "\t" << regions[i]._offset << ", " << regions[i]._size << std::endl;
//  }


  #ifdef Asserts
  int numberOfCells = 0;
  for(int i = 0; i < numberOfRegions; i++) {
    numberOfCells += tarch::la::volume(regions[i]._size);
  }
  #endif

  std::vector<Data> remoteData;

  if(_packCommunication && _messageType == peano::heap::NeighbourCommunication) {
    #if defined(Parallel) && defined(UseBlockedMeshCommunication)
    Serialization::ReceiveBuffer& recvbuffer = peano::parallel::SerializationMap::getInstance().getReceiveBuffer(_remoteRank)[1];
    assertion1(recvbuffer.isBlockAvailable(), "cannot read heap data from Serialization Buffer - not enough blocks");

    Serialization::Block block = recvbuffer.nextBlock();

    size_t numberOfDataElements = block.size() / sizeof(Data::Packed);

    //assertion1(numberOfCellDescriptions > 0, "no cell descriptions, huh? we always send one");
 
    remoteData.resize(numberOfDataElements);

    int block_position = 0;
    //std::cout << " ||||||| unpacking " << numberOfCellDescriptions << " cell descriptions " << std::endl;
    for (size_t i=0; i < numberOfDataElements; i++) {
        Data::Packed packed;
        MPI_Unpack(block.data(), block.size(), &block_position, &packed, 1, Data::Packed::Datatype, MPI_COMM_WORLD );
        remoteData[i] = packed.convert();
    }
    #endif
  } else {

      //Allocate data array
      remoteData = DataHeap::getInstance().receiveData(
                                        _remoteRank,
                                        _position,
                                        _level,
                                        _messageType
                                      );

  }

  int entry = 0;
  for(int i = 0; i < numberOfRegions; i++) {
    peanoclaw::geometry::Region& region = regions[i];

    //U new
    dfor(subcellIndex, region._size) {
      int linearIndex = subgridAccessor.getLinearIndexUNew(region._offset + subcellIndex);
      for(int unknown = 0; unknown < remoteCellDescription.getUnknownsPerSubcell(); unknown++) {
        subgridAccessor.setValueUNew(linearIndex, unknown, remoteData[entry++].getU());
      }

      #if defined(AssertForPositiveValues) && defined(Asserts)
      if(subgrid.isLeaf()) {
        std::stringstream s;
        for(size_t i = 0; i < remoteData.size(); i++) s << remoteData[i].getU() << " ";
        assertion9(
          tarch::la::greater(subgrid.getAccessor().getValueUNew(linearIndex, 0), 0.0),
          subgrid,
          subcellIndex,
          region._offset,
          region._size,
          subgrid.getAccessor().getValueUNew(linearIndex, 0),
          subgrid.toStringUNew(),
          entry,
          s.str(),
          tarch::parallel::Node::getInstance().getRank()
        );
      }
      #endif
    }

    //U old
    dfor(subcellIndex, region._size) {
      int linearIndex = subgridAccessor.getLinearIndexUOld(region._offset + subcellIndex);
      for(int unknown = 0; unknown < remoteCellDescription.getUnknownsPerSubcell(); unknown++) {
//        subgridAccessor.setValueUOldAndResize(linearIndex, unknown, remoteData[entry++].getU());
        subgridAccessor.setValueUOld(linearIndex, unknown, remoteData[entry++].getU());
      }

      #if defined(AssertForPositiveValues) && defined(Asserts)
      if(subgrid.isLeaf()) {
        assertion3(tarch::la::greater(subgrid.getAccessor().getValueUOld(linearIndex, 0), 0.0), subgrid, subcellIndex, subgrid.getAccessor().getValueUOld(linearIndex, 0));
      }
      #endif
    }
  }

//  assertion4(subgrid.getUSize() >= entry, entry, subgrid.getUSize(), subgrid, remoteCellDescription.getSubdivisionFactor());
//  assertion2(DataHeap::getInstance().getData(subgrid.getUIndex()).size() >= entry, DataHeap::getInstance().getData(subgrid.getUIndex()).size(), entry);

  #endif
  logTraceOut("receiveOverlappedCells");
}

void peanoclaw::parallel::SubgridCommunicator::deleteArraysFromSubgrid(
  Patch& subgrid
) {
  logTraceInWith1Argument("deleteArraysFromSubgrid", subgrid);
  if(subgrid.getUIndex() != -1) {
    DataHeap::getInstance().deleteData(subgrid.getUIndex());
  }
  logTraceOut("deleteArraysFromSubgrid");
}
