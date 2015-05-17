#ifndef _PEANO_PARALLEL_SERIALIZATION_MAP_H_
#define _PEANO_PARALLEL_SERIALIZATION_MAP_H_

#if __cplusplus >= 201103L
    #include <array>
#else
    #warning compiler does not use C++11 support. Trying C++ TR1 support
    #include <tr1/array>
#endif

#include <map>
#include <utility>

#include "peano/parallel/Serialization.h"

namespace peano {
    namespace parallel {
        class SerializationMap;
    }
}

class peano::parallel::SerializationMap {
    public:
#if __cplusplus >= 201103L
	    typedef std::array<Serialization::SendBuffer, 2> SendBufferArrayType;
	    typedef std::array<Serialization::ReceiveBuffer, 2> ReceiveBufferArrayType;
#else
	    typedef std::tr1::array<Serialization::SendBuffer, 2> SendBufferArrayType;
	    typedef std::tr1::array<Serialization::ReceiveBuffer, 2> ReceiveBufferArrayType;
#endif

        typedef std::pair<SendBufferArrayType, ReceiveBufferArrayType> ValueType;
        typedef std::map<int, ValueType> BufferMap;


        static SerializationMap& getInstance();

        virtual ~SerializationMap();

        ValueType& operator[](int rank);
       	SendBufferArrayType& getSendBuffer(int rank);
		ReceiveBufferArrayType& getReceiveBuffer(int rank);
 
        inline size_t numberOfBufferPairs() {
            return _bufferMap.size();
        }

        inline BufferMap::iterator begin() {
            return _bufferMap.begin();
        }

        inline BufferMap::iterator end() {
            return _bufferMap.end();
        }

        void exchangeDataWithNeighbors();

    private:
        SerializationMap();

        BufferMap _bufferMap;
};

#endif // _PEANO_PARALLEL_SERIALIZATION_MAP_H_
