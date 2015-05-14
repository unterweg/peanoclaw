/*
 * TreeAwareNodePoolStrategy.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kristof
 */
#ifdef Parallel
#include <mpi.h>
#endif

#include "peanoclaw/parallel/TreeAwareNodePoolStrategy.h"

#include "peano/utils/Dimensions.h"
#include "peano/utils/Globals.h"

#include <sstream>

tarch::logging::Log peanoclaw::parallel::TreeAwareNodePoolStrategy::_log("peanoclaw::parallel::TreeAwareNodePoolStrategy");

peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::Node()
  : _state(Unregistered),
    _rank(-1),
    _level(-100)
{
}

peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::Node(int rank)
  : _state(Registered),
    _rank(rank),
    _level(-100)
{
}

bool peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::isIdle() const {
  return _state == Idle;
}

bool peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::isRegistered() const {
  return _state != Unregistered;
}

bool peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::isWorking() const {
  return _state == Working;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::setIdle() {
  _state = Idle;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::setWorking() {
  _state = Working;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::setRegistered() {
  _state = Registered;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::setUnregistered() {
  _state = Unregistered;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::setLevel(int level) {
  _level = level;
}

int peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::getLevel() const {
  return _level;
}

int peanoclaw::parallel::TreeAwareNodePoolStrategy::Node::getRank() const {
  return _rank;
}

int peanoclaw::parallel::TreeAwareNodePoolStrategy::getRanksPerSubTree(int level) const {
  int workerInSubtree = 0;
  if(level < _numberOfLevels - 2) {
    for(int i = level + 1; i < _numberOfLevels - 1; i++) {
      workerInSubtree += pow(3, DIMENSIONS * (i - level - 1));
    }
    workerInSubtree += pow(3, DIMENSIONS * (_numberOfLevels - 2 - level - 1)) * _workerPerMasterOnLastLevel;
  } else {
    workerInSubtree = 1;
  }
  //TODO unterweg debug
//  std::cout << "masterInSubtree=" << workerInSubtree << " wPMOLL=" << _workerPerMasterOnLastLevel << " level=" << level << " workerInSubtree=" << workerInSubtree << std::endl;
  return workerInSubtree;
}

int peanoclaw::parallel::TreeAwareNodePoolStrategy::getWorkerRankForMaster(int masterRank, int masterLevel) const {
//  if(masterLevel >= _numberOfLevels - 1) {
//    return -1;
//  }

  if(masterRank == 0) {
    return _nodes.size() > 1 ? 1 : -1;
  } else {
    int workerRank = masterRank + 1;
    int workerInSubtree = getRanksPerSubTree(masterLevel);

    for(int i = 0; i < THREE_POWER_D; i++) {
      //TODO unterweg debug
//      std::cout << "worker=" << workerRank << " workerInSubtree=" << workerInSubtree << " masterLevel=" << masterLevel << " nodes.size=" << _nodes.size()
//          << " isIdle=" << _nodes[workerRank].isIdle() << " isRegistered=" << _nodes[workerRank].isRegistered() << " isWorking=" << _nodes[workerRank].isWorking() << std::endl;
      if(workerRank < int(_nodes.size()) && _nodes[workerRank].isIdle()) {
        assertion2(workerRank >= 0 && workerRank < int(_nodes.size()), workerRank, _nodes.size());
        return workerRank;
      }
      workerRank += workerInSubtree;
    }
    return -1;
  }
}

peanoclaw::parallel::TreeAwareNodePoolStrategy::TreeAwareNodePoolStrategy()
  : _tag(-100),
    _numberOfNodes(0),
    _numberOfIdleNodes(0),
    _nodes()
{
  int numberOfRanks;
  #ifdef Parallel
  MPI_Comm_size (MPI_COMM_WORLD, &numberOfRanks);
  #else
  numberOfRanks = 1;
  #endif
  _nodes.reserve(numberOfRanks);
  for(int i = 0; i < numberOfRanks; i++) {
    _nodes.push_back(Node(i));
  }
  _nodes[0].setLevel(-1);

  //Get number of levels
  _numberOfLevels = 0;
  int potentialRanksInWorkerTree = 1;
  while(numberOfRanks > potentialRanksInWorkerTree) {
    potentialRanksInWorkerTree += pow(3, DIMENSIONS * _numberOfLevels);
    _numberOfLevels++;
  }

  //Get number of worker per master that is located on the second last level of the spacetree.
  int numberOfMasters = 1;
  for(int level = 0; level < _numberOfLevels - 1; level++) {
    numberOfMasters += pow(3, DIMENSIONS * level);
  }
  numberOfMasters -= pow(3, DIMENSIONS * (_numberOfLevels - 2));
  assertion3(numberOfRanks-1 >= numberOfMasters, numberOfRanks, numberOfMasters, _numberOfLevels);

  _workerPerMasterOnLastLevel = (numberOfRanks - 1 - numberOfMasters) / pow(3, DIMENSIONS * (_numberOfLevels - 2));

  //TODO unterweg debug
  logInfo("TreeAwareNodePoolStrategy()", "number of ranks: " << numberOfRanks << " number of levels: " << _numberOfLevels
      << " numberOfMasters=" << numberOfMasters
      << " worker per master on last level: " << _workerPerMasterOnLastLevel);
}

peanoclaw::parallel::TreeAwareNodePoolStrategy::~TreeAwareNodePoolStrategy() {
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::setNodePoolTag(int tag) {
  _tag = tag;
}

tarch::parallel::messages::WorkerRequestMessage peanoclaw::parallel::TreeAwareNodePoolStrategy::extractElementFromRequestQueue(
  tarch::parallel::NodePoolStrategy::RequestQueue& queue
) {
  assertion( !queue.empty() );
  tarch::parallel::messages::WorkerRequestMessage result = queue.front();
  queue.pop_front();
  return result;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::fillWorkerRequestQueue(RequestQueue& queue) {
  #ifdef Parallel
  assertion( _tag >= 0 );
  while ( tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_tag, true) ) {
    tarch::parallel::messages::WorkerRequestMessage message;
    message.receive(MPI_ANY_SOURCE,_tag, true, SendAndReceiveLoadBalancingMessagesBlocking);
    queue.push_back( message );
  }
  #endif
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::addNode(
  const tarch::parallel::messages::RegisterAtNodePoolMessage& node
) {
  #ifdef Parallel
  if(!_nodes[node.getSenderRank()].isRegistered()) {
    _nodes[node.getSenderRank()].setRegistered();
  }
  _numberOfNodes++;
  #endif
}


void peanoclaw::parallel::TreeAwareNodePoolStrategy::removeNode( int rank ) {
  logInfo("removeNode(rank)", "Removing node " << rank << " (" << (_numberOfNodes-1) << " nodes remaining");

  if(isIdleNode(rank)) {
    _numberOfIdleNodes--;
  }
  _nodes[rank].setUnregistered();
  _nodes[rank].setLevel(-100);
  _numberOfNodes--;
}


int peanoclaw::parallel::TreeAwareNodePoolStrategy::getNumberOfIdleNodes() const {
  return _numberOfIdleNodes;
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::setNodeIdle( int rank ) {
  logInfo("setNodeIdle(rank)", "Setting node " << rank << " idle");
  if(!_nodes[rank].isIdle()) {
    _numberOfIdleNodes++;
  }
  _nodes[rank].setIdle();

  //TODO unterweg debug
//  std::cout << " this=" << this << std::endl;
//  for(int i = 0; i < _nodes.size(); i++) {
//    std::cout << _nodes[i].isIdle() << " ";
//  }
//  std::cout << std::endl;

  assertion(_nodes[rank].isIdle());
}


int peanoclaw::parallel::TreeAwareNodePoolStrategy::reserveNode(int forMaster) {
  logInfo("reserveNode(master)", "Reserving node for master " << forMaster);

  int masterLevel = _nodes[forMaster].getLevel();
  int workerRank = getWorkerRankForMaster(forMaster, masterLevel);

  //TODO unterweg debug
//  std::cout << "Reserving worker " << workerRank << " for master " << forMaster << std::endl;

  if(workerRank == -1) {
    return -1;
  } else {
    _nodes[workerRank].setLevel(masterLevel + 1);
    _nodes[workerRank].setWorking();
    _numberOfIdleNodes--;
    return workerRank;
  }
}

void peanoclaw::parallel::TreeAwareNodePoolStrategy::reserveParticularNode(int rank) {
  logInfo("reserveParticularNode(rank)", "Reserving node " << rank);
  assertion(_nodes[rank].isIdle());
  _nodes[rank].setWorking();
  _numberOfIdleNodes--;
}

bool peanoclaw::parallel::TreeAwareNodePoolStrategy::isRegisteredNode(int rank) const {
  if(rank >= int(_nodes.size())) {
    return false;
  }

  return _nodes[rank].isRegistered();
}


bool peanoclaw::parallel::TreeAwareNodePoolStrategy::isIdleNode(int rank) const {
  if(!isRegisteredNode(rank)) {
    return false;
  }
  return _nodes[rank].isIdle();
}

int peanoclaw::parallel::TreeAwareNodePoolStrategy::getNumberOfRegisteredNodes() const {
  return _numberOfNodes;
}

std::string peanoclaw::parallel::TreeAwareNodePoolStrategy::toString() const {
  return "TreeAwareNodePoolStrategy";
}

bool peanoclaw::parallel::TreeAwareNodePoolStrategy::hasIdleNode(int forMaster) const {
  //During termination Peano asks forMaster=AnyMaster
  if(forMaster == tarch::parallel::NodePoolStrategy::AnyMaster) {
    return _numberOfNodes > 0;
  } else {
    int masterLevel = _nodes[forMaster].getLevel();

    //TODO unterweg debug
//    std::cout << "Idle worker " << getWorkerRankForMaster(forMaster, masterLevel) << " for master " << forMaster << " this=" << this << std::endl;
//    for(int i = 0; i < _nodes.size(); i++) {
//      std::cout << _nodes[i].isIdle() << " ";
//    }
//    std::cout << std::endl;

    return getWorkerRankForMaster(forMaster, masterLevel) != -1;
  }
}

int peanoclaw::parallel::TreeAwareNodePoolStrategy::removeNextIdleNode() {
  for(size_t rank = 1; rank < _nodes.size(); rank++) {
    if(_nodes[rank].isRegistered()) {
      removeNode(rank);
      return rank;
    }
  }

  assertionFail("No next idle node to remove... number of registered nodes=" << _numberOfNodes);
  return -1;
}





