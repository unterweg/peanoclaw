#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/Assertions.h"

#include "SerializationMap.h"

peano::parallel::SerializationMap& peano::parallel::SerializationMap::getInstance() {
    static SerializationMap singleton;
    return singleton;
}

peano::parallel::SerializationMap::~SerializationMap() {

}

peano::parallel::SerializationMap::SerializationMap() {

}

peano::parallel::SerializationMap::ValueType& peano::parallel::SerializationMap::operator[](int rank) {
    BufferMap::iterator it = _bufferMap.find(rank);
    if (it == _bufferMap.end()) {
      // object does not exist yet, insert a new set of buffers

      Serialization::SendBuffer sendbuffer_temp(true);
      Serialization::ReceiveBuffer recvbuffer_temp(true);

      peano::parallel::SerializationMap::SendBufferArrayType sendbuffer_array = {sendbuffer_temp, sendbuffer_temp};
      peano::parallel::SerializationMap::ReceiveBufferArrayType recvbuffer_array = {recvbuffer_temp, recvbuffer_temp};

      peano::parallel::SerializationMap::ValueType maptype(
          sendbuffer_array, recvbuffer_array
      );

      std::pair<const int, peano::parallel::SerializationMap::ValueType> temp_value(rank, maptype);
      std::pair<BufferMap::iterator, bool> result = _bufferMap.insert(temp_value);
      it = result.first;
    }
    return it->second;
}

peano::parallel::SerializationMap::SendBufferArrayType& peano::parallel::SerializationMap::getSendBuffer(int rank) {
    return this->operator[](rank).first;
}

peano::parallel::SerializationMap::ReceiveBufferArrayType& peano::parallel::SerializationMap::getReceiveBuffer(int rank) {
    return this->operator[](rank).second;
}

void peano::parallel::SerializationMap::exchangeDataWithNeighbors() {
#ifdef Parallel
    // first determine maximum number of buffer pairs (some might have no data, though)
    // and allocate a MPI_RequestBuffer of twice the size (one for each send AND receive buffer)
    std::vector<MPI_Request> requestBuffer(2*numberOfBufferPairs(), MPI_REQUEST_NULL);

    // initiate MPI calls for deterministic buffers
    size_t requestBufferPosition = 0;
    for (BufferMap::iterator it = begin(); it != end(); it++)
    {
        int neighbourRank = it->first;
		ValueType& bufferPair = it->second;

        // deterministic buffers: we are able to determine the size for the required receive buffers 
        // independently from other processes
        Serialization::SendBuffer& deterministicSendBuffer = bufferPair.first[0];
        Serialization::ReceiveBuffer& deterministicReceiveBuffer = bufferPair.second[0];

        assertion(deterministicSendBuffer.verifyBlocks());

        // nondeterministic buffers: the opposite of deterministic buffers, hence we need some
        // additional information to set an appropriate buffer size
        // in our case: we will put the size of the non deterministic buffer as payload inside the 
        // deterministic buffer stream
        Serialization::SendBuffer& nondeterministicSendBuffer = bufferPair.first[1];
        Serialization::ReceiveBuffer& nondeterministicReceiveBuffer = bufferPair.second[1];

        assertion(nondeterministicSendBuffer.verifyBlocks());
        
        Serialization::Block block = deterministicSendBuffer.reserveBlock(sizeof(size_t));
        const size_t nondeterministicSendBufferSize = nondeterministicSendBuffer.size();
        block << nondeterministicSendBufferSize;


        deterministicReceiveBuffer.reset(deterministicSendBuffer.size());

        //std::cout << "my rank " << tarch::parallel::Node::getInstance().getRank() << " ------------ SENDING DETERMINISTIC TO " << neighbourRank << " size " << deterministicSendBuffer.size() << std::endl;
        //std::cout << "my rank " << tarch::parallel::Node::getInstance().getRank() << " ------------ EXPECTING DETERMINISTIC FROM " << neighbourRank << " size " << deterministicSendBuffer.size() << std::endl;

        if (deterministicSendBuffer.size() > 0) {
            MPI_Isend(deterministicSendBuffer.data(), deterministicSendBuffer.size(), MPI_CHAR, neighbourRank, 1337, MPI_COMM_WORLD, requestBuffer.data()+2*requestBufferPosition);
            MPI_Irecv(deterministicReceiveBuffer.data(), deterministicSendBuffer.size(), MPI_CHAR, neighbourRank, 1337, MPI_COMM_WORLD, requestBuffer.data()+2*requestBufferPosition+1);
        }

        requestBufferPosition++;
    }

    // wait for MPI actions to finish data exchange of deterministic buffers
    MPI_Waitall(2*numberOfBufferPairs(), requestBuffer.data(), MPI_STATUSES_IGNORE);

    #ifdef Asserts
    for (BufferMap::iterator it = begin(); it != end(); it++)
    {
      int neighbourRank = it->first;
      ValueType& bufferPair = it->second;
      Serialization::ReceiveBuffer& deterministicReceiveBuffer = bufferPair.second[0];
      assertion(deterministicReceiveBuffer.verifyBlocks());
    }
    #endif

    // clear deterministic send buffer, keep determinisitc receive buffers, 
    // and initiate non deterministic MPI actions
    requestBufferPosition = 0;
    for (BufferMap::iterator it = begin(); it != end(); it++) {
        int neighbourRank = it->first;
		ValueType& bufferPair = it->second;

        Serialization::SendBuffer& deterministicSendBuffer = bufferPair.first[0];
        Serialization::ReceiveBuffer& deterministicReceiveBuffer = bufferPair.second[0];

        Serialization::SendBuffer& nondeterministicSendBuffer = bufferPair.first[1];
        Serialization::ReceiveBuffer& nondeterministicReceiveBuffer = bufferPair.second[1];

        deterministicSendBuffer.reset();

        // extract size for nondeterministic buffer from deterministic reveive buffer
        assertion1(deterministicReceiveBuffer.isBlockAvailable(), "block for nondeterministic buffer size is missing!");
        Serialization::Block block = deterministicReceiveBuffer.nextBlock();

        assertion2(block.size() == sizeof(size_t), "block for nondeterministic buffer size has wrong size", block.size());

        size_t nondeterministicReceiveBufferSize = 0;
        block >> nondeterministicReceiveBufferSize;

        nondeterministicReceiveBuffer.reset(nondeterministicReceiveBufferSize);

        if (nondeterministicSendBuffer.size() > 0) {
            //std::cout << "my rank " << tarch::parallel::Node::getInstance().getRank() << " ------------ SENDING NONDETERMINISTIC TO " << neighbourRank << " size " << nondeterministicSendBuffer.size() << std::endl;
            MPI_Isend(nondeterministicSendBuffer.data(), nondeterministicSendBuffer.size(), MPI_CHAR, neighbourRank, 1337, MPI_COMM_WORLD, requestBuffer.data()+2*requestBufferPosition);
        } else {
            //std::cout << "my rank " << tarch::parallel::Node::getInstance().getRank() << " ------------ NOT SENDING NONDETERMINISTIC TO " << neighbourRank << std::endl;
            requestBuffer[2*requestBufferPosition] = MPI_REQUEST_NULL;
        }

        if (nondeterministicReceiveBuffer.size() > 0) {
            //std::cout << "my rank " << tarch::parallel::Node::getInstance().getRank() << " ------------ EXPECTING NONDETERMINISTIC FROM " << neighbourRank << " size " << nondeterministicReceiveBuffer.size() << std::endl;
            MPI_Irecv(nondeterministicReceiveBuffer.data(), nondeterministicReceiveBuffer.size(), MPI_CHAR, neighbourRank, 1337, MPI_COMM_WORLD, requestBuffer.data()+2*requestBufferPosition+1);
        } else {
            //std::cout << "my rank " << tarch::parallel::Node::getInstance().getRank() << " ------------ NOT EXPECTING NONDETERMINISTIC FROM " << neighbourRank << std::endl;
            requestBuffer[2*requestBufferPosition+1] = MPI_REQUEST_NULL;
        }
 
        requestBufferPosition++;
    }

    // wait for MPI actions to finish data exchange of nondeterministic buffers
    MPI_Waitall(2*numberOfBufferPairs(), requestBuffer.data(), MPI_STATUSES_IGNORE);
 
    // clear nondeterministic send buffer
    requestBufferPosition = 0;
    for (BufferMap::iterator it = begin(); it != end(); it++) {
		peano::parallel::SerializationMap::ValueType& bufferPair = it->second;

        Serialization::SendBuffer& nondeterministicSendBuffer = bufferPair.first[1];
        nondeterministicSendBuffer.reset();

        #ifdef Asserts
        Serialization::ReceiveBuffer& nondeterministicReceiveBuffer = bufferPair.second[1];
        assertion(nondeterministicReceiveBuffer.verifyBlocks());
        #endif
    }
#endif
}
