#ifndef _PEANO_PARALLEL_MESHCOMMUNICATION_H_
#define _PEANO_PARALLEL_MESHCOMMUNICATION_H_

#include <set>

#include "peano/parallel/messages/LoadBalancingMessage.h"
#include "peano/grid/SingleLevelEnumerator.h"
#include "peano/grid/SingleElementVertexEnumerator.h"

#if defined(Parallel)

class MeshCommunication {
    private:
      static tarch::logging::Log _log;

    public:
        template<typename VertexData, typename CellData> 
        struct MasterWorkerForkJoin {
            // TODO: we need true MPI buffers
            CellData _cellData;
            VertexData _vertexData[TWO_POWER_D];
#if DIMENSIONS==2
            unsigned char adjacencyBitSet; // encode TWO_POWER_D = 4 entries
#elif DIMENSIONS==3
            unsigned char adjacencyBitset; // encode TWO_POWER_D = 8 entries
#else
            #error not yet implemented for current number of dimensions
#endif
        } __attribute__((packed));

        template<typename VertexData, typename CellData> 
        struct DataExchangeAscend {
            CellData _cellData;
            VertexData _vertexData[TWO_POWER_D];
        } __attribute__((packed));

        template<typename VertexData, typename CellData> 
        struct DataExchangeDescend {
            CellData _cellData[2];
            VertexData _vertexData[TWO_POWER_D];
        } __attribute__((packed));


        template<typename Vertex, typename Cell> 
        static 
#if defined(Debug)
        int 
#else
        void
#endif
        receiveInitialMeshDataFromMaster(
          std::bitset<NUMBER_OF_VERTICES_PER_ELEMENT>& adjacencyBitset, 
          Cell& receivedCell, 
          Vertex* receivedVertices
        );
 
        template<typename Vertex, typename Cell> 
        static void sendInitialMeshDataToWorker(
          int currentWorker,
          const std::bitset<NUMBER_OF_VERTICES_PER_ELEMENT>& adjacencyBitset, 
          Cell& fineGridCell, 
          Vertex* fineGridVertices,
          const peano::grid::SingleLevelEnumerator&              fineGridVerticesEnumerator
        );

        template<typename Vertex, typename Cell> 
        static 
#if defined(Debug)
        int 
#else
        void
#endif
         receiveFinalMeshDataFromWorker(
          int currentWorker,
          std::bitset<NUMBER_OF_VERTICES_PER_ELEMENT>& adjacencyBitset, 
          Cell& receivedCell, 
          Vertex* receivedVertices
        );

        template<typename Vertex, typename Cell> 
        static void sendFinalMeshDataToMaster(
          const std::bitset<NUMBER_OF_VERTICES_PER_ELEMENT>& adjacencyBitset, 
          Cell& fineGridCell, 
          Vertex* fineGridVertices,
          const peano::grid::SingleLevelEnumerator&              fineGridVerticesEnumerator
        );

        template<typename Vertex, typename Cell> 
        static void receiveMeshDataFromWorker(
          int currentWorker,
          Cell& receivedCell, 
          Vertex* receivedVertices
        );

        template<typename Vertex, typename Cell> 
        static void sendMeshDataToMaster(
          Cell& receivedCell, 
          Vertex* receivedVertices,
          const peano::grid::SingleLevelEnumerator&              fineGridVerticesEnumerator
        );
     
        template<typename Vertex, typename Cell>
        static void sendMeshDataToWorker(
           int currentWorker,
           Cell& coarseGridCell,
           Vertex* coarseGridVertices,
           const peano::grid::SingleLevelEnumerator&              coarseGridVerticesEnumerator,
           Cell& fineGridCell,
           Vertex* fineGridVertices,
           const peano::grid::SingleLevelEnumerator&              fineGridVerticesEnumerator
        );

        template<typename Vertex, typename Cell>
        static void receiveMeshDataFromMaster(
           Cell& coarseGridCell,
           Vertex* coarseGridVertices,
           Cell& fineGridCell,
           Vertex* fineGridVertices
        );

        template<typename State> 
        static void receiveStateDataFromWorker(
          int currentWorker,
          State& workerState
        );

        template<typename State>
        static void sendStateDataToMaster(
           State& state
        );

        template<typename State> 
        static void sendStateDataToWorker(
          int currentWorker,
          peano::parallel::messages::LoadBalancingMessage& loadBalancingMessage,
          State& workerState
        );

        template<typename State>
        static void receiveStateDataFromMaster(
           peano::parallel::messages::LoadBalancingMessage& loadBalancingMessage,
           State& state
        );

        template<typename Vertex>
        static void receiveVertexFromNeighbour(
            int neighbourRank,
            Vertex& vertex
        );
 
        template<typename Vertex>
        static void sendVertexToNeighbour(
            int neighbourRank,
            const Vertex& vertex
        );

        template<typename Data>
        static void sendAsynchronousAndBlocking(
                Data& data,
                int destination, 
                int tag, 
                bool exchangeOnlyAttributesMarkedWithParallelise,
                bool exchangeDataBlocking
        );
 
        template<typename Data>
        static void receiveAsynchronousAndBlocking(
                Data& data,
                int destination, 
                int tag, 
                bool exchangeOnlyAttributesMarkedWithParallelise,
                bool exchangeDataBlocking
        );

        static void sendAsynchronousAndBlocking(
                void *data,
                int size,
                const MPI_Datatype& datatype,
                int destination, 
                int tag,
                bool exchangeDataBlocking
        );
 
        static void receiveAsynchronousAndBlocking(
                void *data,
                int size,
                const MPI_Datatype& datatype,
                int source, 
                int tag,
                bool exchangeDataBlocking
        );
};

#include "peano/parallel/MeshCommunication.cpph"

#endif

#endif // _PEANO_PARALLEL_MESHCOMMUNICATION_H_
