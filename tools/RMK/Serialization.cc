#include "Serialization.h"

#include "tarch/Assertions.h"

Serialization::Block::Block(size_t blocksize, unsigned char * const block_ptr) 
:   blocksize(blocksize), 
    block_ptr(block_ptr),
    position(0)
{

}

Serialization::Block::~Block() {

}

void Serialization::Block::pack(size_t count, const unsigned char* buffer) {
    // TODO: assertion: position < blocksize
    memcpy(block_ptr+position, buffer, count);
    position += count;
}

void Serialization::Block::unpack(size_t count, unsigned char* buffer) {
    // TODO: assertion: position < blocksize
    memcpy(buffer, block_ptr+position, count);
    position += count;
}

/////////////////////////////////////////////////////////////////


// stack = true: put size at end of block
// stack = false: put size at beginning of block
Serialization::SendBuffer::SendBuffer(bool stack)
    : stack(stack)
{
    reset();
}

Serialization::SendBuffer::~SendBuffer() {

}

// every section is built of payload data + size information
Serialization::Block Serialization::SendBuffer::reserveBlock(size_t blocksize) {
    const size_t total_blocksize = blocksize + sizeof(size_t);
    const size_t current_storageSize = storage.size();
    const size_t required_storageSize = position + total_blocksize;
    if (required_storageSize > current_storageSize) {
        storage.resize(required_storageSize);
    }

    unsigned char* const totalblock_ptr = storage.data() + position;
    unsigned char* result_ptr;
    size_t* size_ptr;

    if (stack) {
      result_ptr = totalblock_ptr;
      size_ptr = reinterpret_cast<size_t*>(totalblock_ptr + blocksize);
    } else {
      assertionFail("No stack.");
      size_ptr = reinterpret_cast<size_t*>(totalblock_ptr);
      result_ptr = totalblock_ptr + sizeof(size_t);
    }
    *size_ptr = blocksize;
    
    position += total_blocksize;
    return Serialization::Block(blocksize, result_ptr);
}

void Serialization::SendBuffer::reset() {
    position = 0;
    storage.clear();
}

bool Serialization::SendBuffer::verifyBlocks() const {
  if(stack) {
    size_t pos = position;
    while(pos >= sizeof(size_t)) {
      //TODO unterweg debug
//      std::cout << "pos=" << pos << " storage.size=" << storage.size() << std::endl;
      assertion3(pos <= position, pos, position, size());
      size_t blockSize = *(reinterpret_cast<const size_t*>(storage.data() + pos - sizeof(size_t)));
//      std::cout << "\tblockSize=" << blockSize << std::endl;
      pos -= (blockSize + sizeof(size_t));
    }
    assertionEquals(pos, 0);
    return pos == 0;
  } else {
    assertionFail("Not implemented...");
    return false;
  }
}

/////////////////////////////////////////////////////////////////

Serialization::ReceiveBuffer::ReceiveBuffer(bool stack)
    : stack(stack),
      position(0)
{
}

Serialization::ReceiveBuffer::~ReceiveBuffer() {

}

bool Serialization::ReceiveBuffer::isBlockAvailable() {
    bool result = !storage.empty();
    if (stack) {
        result &= (position != 0);
    } else {
        result &= (position != size());
    }
    return result;
}

Serialization::Block Serialization::ReceiveBuffer::nextBlock() { 
    unsigned char * const totalblock_ptr = storage.data() + position;
    size_t* size_ptr;
    unsigned char* result_ptr;

    size_t blocksize = 0;
    if (stack) {
        size_ptr = reinterpret_cast<size_t*>(totalblock_ptr - sizeof(size_t));
        blocksize = *size_ptr;
        result_ptr = totalblock_ptr - (sizeof(size_t) + blocksize);

        position -= (sizeof(size_t) + blocksize);
    } else {
        assertionFail("No stack.");
        size_ptr = reinterpret_cast<size_t*>(totalblock_ptr);
        blocksize = *size_ptr;
        result_ptr = totalblock_ptr + sizeof(size_t);

        position += (sizeof(size_t) + blocksize);
    }
    return Serialization::Block(blocksize, result_ptr);
}

void Serialization::ReceiveBuffer::reset(size_t storageSize) {
    storage.clear();
    storage.resize(storageSize);

    if (stack) {
        position = size();
    } else {
        position = 0;
    }
}

bool Serialization::ReceiveBuffer::verifyBlocks() const {
  if(stack) {
    size_t pos = size();
    while(pos >= sizeof(size_t)) {
      //TODO unterweg debug
//      std::cout << "pos=" << pos << " storage.size=" << storage.size() << std::endl;
      assertion2(pos <= size(), pos, size());
      size_t blockSize = *(reinterpret_cast<const size_t*>(storage.data() + pos - sizeof(size_t)));
//      std::cout << "\tblockSize=" << blockSize << std::endl;
      pos -= (blockSize + sizeof(size_t));
    }
    assertionEquals(pos, 0);
    return pos == 0;
  } else {
    assertionFail("Not implemented...");
    return false;
  }
}


