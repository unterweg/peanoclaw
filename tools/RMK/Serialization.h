#ifndef __SERIALIZATION_H__
#define __SERIALIZATION_H__

#include <vector>
#include <cstring> // memcpy
#include <cstddef> // size_t

namespace Serialization {
    class Block;
    class SendBuffer;
    class ReceiveBuffer;
}

// keep this as generic if possible
// TODO: split into read and write block?
class Serialization::Block {
    public:
        Block(size_t blocksize, unsigned char * const block_ptr);
        virtual ~Block();

        void pack(size_t count, const unsigned char* buffer);
        void unpack(size_t count, unsigned char* buffer);

        inline size_t size() {
            return blocksize;
        }

        inline unsigned char* data() {
            return block_ptr;
        }

    private:
        const size_t blocksize;
        unsigned char * const block_ptr;
 
        size_t position;
};

class Serialization::SendBuffer {
    public:
        // stack = true: put size at end of block
        // stack = false: put size at beginning of block
        SendBuffer(bool stack);
        virtual ~SendBuffer();

        // every section is built of payload data + size information
        Serialization::Block reserveBlock(size_t blocksize);

        inline size_t size() const {
            return storage.size();
        }

        void reset();

        inline unsigned char* data() {
            return storage.data();
        }

        bool verifyBlocks() const;

    private:
        bool stack;
        size_t position;
        std::vector<unsigned char> storage;
};

class Serialization::ReceiveBuffer {
    public:
        ReceiveBuffer(bool stack);
        virtual ~ReceiveBuffer();

        inline size_t size() const {
            return storage.size();
        }

        bool isBlockAvailable();
        Serialization::Block nextBlock();
        void reset(size_t storageSize);
 
        inline unsigned char* data() {
            return storage.data();
        }

        bool verifyBlocks() const;

    private:
        bool stack;
        size_t position;
        std::vector<unsigned char> storage;
};

template<typename Datatype>
Serialization::Block& operator<<(Serialization::Block& block, const Datatype& data) {
    block.pack(sizeof(Datatype), reinterpret_cast<const unsigned char*>(&data));
    return block;
}

template<typename Datatype>
Serialization::Block& operator>>(Serialization::Block& block, Datatype& data) {
    block.unpack(sizeof(Datatype), reinterpret_cast<unsigned char*>(&data));
    return block;
}

#endif // __SERIALIZATION_H__
