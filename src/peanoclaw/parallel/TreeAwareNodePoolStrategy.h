/*
 * TreeAwareNodePoolStrategy.h
 *
 *  Created on: Apr 1, 2015
 *      Author: kristof
 */

#ifndef PEANOCLAW_PARALLEL_TREEAWARENODEPOOLSTRATEGY_H_
#define PEANOCLAW_PARALLEL_TREEAWARENODEPOOLSTRATEGY_H_

namespace peanoclaw {
  namespace parallel {
    class TreeAwareNodePoolStrategy;
  }
}

#include "tarch/parallel/NodePoolStrategy.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#include <vector>
#include <map>

class peanoclaw::parallel::TreeAwareNodePoolStrategy : public tarch::parallel::NodePoolStrategy {

  class Node {
    private:
      enum State {
        Unregistered,
        Registered,
        Idle,
        Working
      };
      State _state;

      int   _rank;
      int   _level;

    public:
      Node();
      Node(int rank);

      //State
      bool isIdle() const;
      bool isRegistered() const;
      bool isWorking() const;
      void setWorking();
      void setIdle();
      void setRegistered();
      void setUnregistered();

      //Level
      void setLevel(int level);
      int getLevel() const;

      int getRank() const;

      bool operator<(const Node& right) const;
  };

  /**
   * Logging Device
   */
  static tarch::logging::Log _log;

  /**
   * Tag on which the node pool works
   */
  int _tag;

  int _numberOfNodes;

  int _numberOfIdleNodes;

  std::vector<Node> _nodes;

  int _numberOfLevels;

  int _workerPerMasterOnLastLevel;

  int getRanksPerSubTree(int level) const;

  int getWorkerRankForMaster(int masterRank, int masterLevel) const;

public:
  /**
   * Constructor
   *
   * Construct all the attributes.
   */
  TreeAwareNodePoolStrategy();
  virtual ~TreeAwareNodePoolStrategy();

  virtual void setNodePoolTag(int tag);

  virtual tarch::parallel::messages::WorkerRequestMessage extractElementFromRequestQueue(tarch::parallel::NodePoolStrategy::RequestQueue& queue);

  virtual void fillWorkerRequestQueue(RequestQueue& queue);

  virtual void addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node );

  virtual void removeNode( int rank );

  virtual int getNumberOfIdleNodes() const;

  virtual void setNodeIdle( int rank );

  virtual int reserveNode(int forMaster);

  virtual void reserveParticularNode(int rank);

  virtual bool isRegisteredNode(int rank) const;

  virtual bool isIdleNode(int rank) const;

  virtual int getNumberOfRegisteredNodes() const;

  virtual std::string toString() const;

  virtual bool hasIdleNode(int forMaster) const;

  virtual int removeNextIdleNode();
};



#endif /* PEANOCLAW_PARALLEL_TREEAWARENODEPOOLSTRATEGY_H_ */
