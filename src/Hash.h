#ifndef HASH_INTERFACE_H_
#define HASH_INTERFACE_H_

#include "../util/pair.h"
template<class T>
class Hash {
  public:
    Hash(void) = default;
    ~Hash(void) = default;
    virtual int Insert(T, Value_t) = 0;
    virtual bool Delete(T) = 0;
    virtual Value_t Get(T) = 0;
};

#endif  // _HASH_INTERFACE_H_