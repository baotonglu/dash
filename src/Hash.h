
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
#ifndef HASH_INTERFACE_H_
#define HASH_INTERFACE_H_

#include "../util/pair.h"
#ifdef PMEM
#include <libpmemobj.h>
#endif

/*
* Parent function of all hash indexes
* Used to define the interface of the hash indexes
*/

template <class T>
class Hash {
 public:
  Hash(void) = default;
  ~Hash(void) = default;
  /*0 means success insert, -1 means this key already exist, directory return*/
  virtual int Insert(T, Value_t) = 0;
  virtual int Insert(T, Value_t, bool) = 0;
 
  virtual void bootRestore(){

  };
  virtual void reportRestore(){

  };
  virtual bool Delete(T) = 0;
  virtual bool Delete(T, bool) = 0;
  virtual bool Get(T, Value_t*) = 0;
  virtual bool Get(T key, Value_t*, bool is_in_epoch) = 0;
  virtual void getNumber() = 0;
};

#endif  // _HASH_INTERFACE_H_
