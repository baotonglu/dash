#ifndef HASH_INTERFACE_H_
#define HASH_INTERFACE_H_

#include "../util/pair.h"
#ifdef PMEM
#include <libpmemobj.h>
#endif

template<class T>
class Hash {
  public:
    Hash(void) = default;
    ~Hash(void) = default;
    virtual void Insert(T, Value_t) = 0;
    virtual void Insert(T, Value_t, bool) = 0;
    virtual void Insert(T, Value_t, int){

    };
    virtual bool Delete(T) = 0;
    virtual bool Delete(T, bool) = 0;
    virtual Value_t Get(T) = 0;
    virtual Value_t Get(T key, bool is_in_epoch) = 0;
    virtual void Recovery() = 0;
    virtual void getNumber() = 0;
};

#endif  // _HASH_INTERFACE_H_