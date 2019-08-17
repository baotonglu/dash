#ifndef CCEH_H_
#define CCEH_H_

#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <cmath>
#include <thread>
#include <shared_mutex>
#include <bitset>
#include <cassert>
#include <unordered_map>
#include "../../util/hash.h"
#include "../../util/pair.h"
#include "../../util/persist.h"
#include "../allocator.h"

#ifdef PMEM
#include <libpmemobj.h>
#endif

//#define PERSISTENT_LOCK 1

#define INPLACE 1
template<class T>
struct _Pair{
  T key;
  Value_t value;
};

//const size_t kCacheLineSize = 64;
constexpr size_t kSegmentBits = 8;
constexpr size_t kMask = (1 << kSegmentBits)-1;
constexpr size_t kShift = kSegmentBits;
constexpr size_t kSegmentSize = (1 << kSegmentBits) * 16 * 4;
constexpr size_t kNumPairPerCacheLine = kCacheLineSize/16;
constexpr size_t kNumCacheLine = 4;

uint64_t clflushCount;

template<class T>
struct Segment {
  static const size_t kNumSlot = kSegmentSize/sizeof(_Pair<T>);

  Segment(void)
  : local_depth{0}, sema{0}, count{0}, seg_lock{0}, mutex(), rwlock()
  { 
    memset((void*)&_[0],255,sizeof(_Pair<T>)*kNumSlot);
  }

  Segment(size_t depth)
  :local_depth{depth}, sema{0}, count{0}, seg_lock{0}, mutex(), rwlock()
  {
    memset((void*)&_[0],255,sizeof(_Pair<T>)*kNumSlot);
  }

  static void New(Segment **seg, size_t depth){
#ifdef PMEM
  auto callback = [](PMEMobjpool *pool, void *ptr, void *arg){
    auto value_ptr = reinterpret_cast<size_t*>(arg);
    auto seg_ptr = reinterpret_cast<Segment*>(ptr);
    seg_ptr->local_depth = *value_ptr;
    seg_ptr->sema = 0;
    seg_ptr->count = 0;
    seg_ptr->seg_lock = 0;
    //new &(seg_ptr->mutex) std::shared_mutex();
    memset((void *)&seg_ptr->mutex, 0, sizeof(std::shared_mutex));
    memset((void *)&seg_ptr->rwlock, 0, sizeof(PMEMrwlock));
    memset((void*)&seg_ptr->_[0],255,sizeof(_Pair<T>)*kNumSlot);
    return 0;
  };
  Allocator::Allocate((void**)seg, kCacheLineSize, sizeof(Segment),
                      callback, reinterpret_cast<void*>(&depth));
#else
  Allocator::ZAllocate((void **)seg, kCacheLineSize, sizeof(Segment));
  new (*seg) Segment(depth);
#endif
  }

  ~Segment(void) {}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

  int Insert(PMEMobjpool *,T, Value_t, size_t, size_t);
  int Insert4split(T, Value_t, size_t);
  bool Put(T, Value_t, size_t);
  Segment* Split(PMEMobjpool *);
  
  void get_lock(PMEMobjpool* pop) {
    #ifdef PERSISTENT_LOCK
    pmemobj_rwlock_wrlock(pop, &rwlock);
    #else
    uint64_t temp = 0;
    while(!CAS(&seg_lock, &temp, 1)){
      temp = 0;
    }
    #endif
  }

  void release_lock(PMEMobjpool* pop) {
    #ifdef PERSISTENT_LOCK
    pmemobj_rwlock_unlock(pop, &rwlock);
    #else
    uint64_t temp = 1;
    while(!CAS(&seg_lock, &temp, 0)){
      temp = 1;
    }
    #endif
  }

  void get_rd_lock(PMEMobjpool* pop){
    #ifdef PERSISTENT_LOCK
    pmemobj_rwlock_rdlock(pop, &rwlock);
    #else
    uint64_t temp = 0;
    while(!CAS(&seg_lock, &temp, 1)){
      temp = 0;
    }
    #endif
  }

  void release_rd_lock(PMEMobjpool* pop){
    #ifdef PERSISTENT_LOCK
    pmemobj_rwlock_unlock(pop, &rwlock);
    #else
    uint64_t temp = 1;
    while(!CAS(&seg_lock, &temp, 0)){
      temp = 1;
    }
    #endif
  }

  bool try_get_lock(PMEMobjpool* pop){
    #ifdef PERSISTENT_LOCK   
    if (pmemobj_rwlock_trywrlock(pop, &rwlock) == 0)
    {
      return true;
    }
    return false;
    #else
    uint64_t temp = 0;
    return CAS(&seg_lock, &temp, 1);
    #endif
  }

  bool try_get_rd_lock(PMEMobjpool* pop){
    #ifdef PERSISTENT_LOCK
    if (pmemobj_rwlock_tryrdlock(pop, &rwlock) == 0)
    {
      return true;
    }
    return false;
    #else
    uint64_t temp = 0;
    return CAS(&seg_lock, &temp, 1);
    #endif
  }

  _Pair<T> _[kNumSlot];
  size_t local_depth;
  int64_t sema = 0;
  size_t pattern = 0;
  int count=0;
  std::shared_mutex mutex;
  uint64_t seg_lock;
  PMEMrwlock rwlock;
};

template<class T>
struct Seg_array{
  Segment<T> **_; /*physical addr*/
  size_t global_depth;
  Seg_array(size_t capacity){
    global_depth = static_cast<size_t>(log2(capacity));
    _ = new Segment<T>*[capacity];
  }

  static void New(Seg_array **sa, size_t capacity){
    Segment<T> **seg_array{nullptr};
    Allocator::ZAllocate((void**)&seg_array, kCacheLineSize,
                           sizeof(uint64_t) * capacity);
#ifdef PMEM
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg){
      auto value_ptr = reinterpret_cast<std::pair<Segment<T>**, size_t> *>(arg);
      auto sa_ptr = reinterpret_cast<Seg_array *>(ptr);
      sa_ptr->_ = value_ptr->first;
      sa_ptr->global_depth = static_cast<size_t>(log2(value_ptr->second));
      return 0;
    };
    auto call_args = std::make_pair(seg_array, capacity);
    Allocator::Allocate((void**)sa, kCacheLineSize, sizeof(Seg_array), callback,
                          reinterpret_cast<void*>(&call_args));
#else
    Allocator::ZAllocate((void**)sa, kCacheLineSize, sizeof(Seg_array));
    new (*sa) Seg_array(capacity);
#endif
  }

};

template<class T>
struct Directory {
  static const size_t kDefaultDirectorySize = 1024;
  Seg_array<T>* sa;
  Seg_array<T> *new_sa;
  size_t capacity;
  bool lock;
  int sema = 0;

  Directory(Seg_array<T> *_sa) {
    capacity = kDefaultDirectorySize;
    sa = _sa;
    new_sa = nullptr;
    lock = false;
    sema = 0;
  }

  Directory(size_t size, Seg_array<T> *_sa) {
    capacity = size;
    //sa = new Seg_array<T>(capacity);
    sa = _sa;
    new_sa = nullptr;
    lock = false;
    sema = 0;
  }

  static void New(Directory** dir, size_t capacity){
    Seg_array<T>* temp_sa;
    Seg_array<T>::New(&temp_sa, capacity);
#ifdef PMEM
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg){
      auto value_ptr 
          = reinterpret_cast<std::pair<size_t, Seg_array<T> *> *>(arg);
      auto dir_ptr = reinterpret_cast<Directory *>(ptr);
      dir_ptr->capacity = value_ptr->first;
      dir_ptr->sa = value_ptr->second;
      dir_ptr->lock = false;
      dir_ptr->sema = 0; 
      dir_ptr = nullptr;
      return 0;
    };

    auto call_args = std::make_pair(capacity, temp_sa);
    Allocator::Allocate((void **)dir, kCacheLineSize, sizeof(Directory), 
                        callback, reinterpret_cast<void*>(&call_args));
#else
    Allocator::ZAllocate((void **)dir, kCacheLineSize, sizeof(Directory));
    new (*dir) Directory(capacity, temp_sa); 
#endif
  }

  ~Directory(void) {
    
  }

  void get_item_num(){
    /*
    size_t count = 0;
    size_t seg_num = 0;
    Seg_array<T> *seg = sa;
    Segment<T>** dir_entry = seg->_;
    Segment<T> *ss;
    auto global_depth = seg->global_depth;
    size_t depth_diff;
    for (int i = 0; i < capacity;)
    {
      ss = dir_entry[i];
      depth_diff = global_depth - ss->local_depth;
      count += ss->count;
      seg_num++;
      i += pow(2, depth_diff);
    }
    printf("#items: %lld\n", count);
    printf("load_factor: %f\n", (double)count/(seg_num*256*4));*/
  }

  bool Acquire(void) {
    bool unlocked = false;
    return CAS(&lock, &unlocked, true);
  }
  
  bool Release(void) {
    bool locked = true;
    return CAS(&lock, &locked, false);
  }
  
  void SanityCheck(void*);
};

template<class T>
class CCEH{
  public:
    CCEH(void);
    CCEH(int);
    ~CCEH(void);
    void Insert(T key, Value_t value);
    bool InsertOnly(T, Value_t);
    bool Delete(T);
    Value_t Get(T);
    Value_t FindAnyway(T);
    double Utilization(void);
    size_t Capacity(void);
    bool Recovery(void);
    void Directory_Doubling(int x, Segment<T> *s0, Segment<T> *s1);
    void Directory_Update(int x, Segment<T> *s0, Segment<T> *s1);
    void Lock_Directory();
    void Unlock_Directory();
    void getNumber(){
      dir->get_item_num();
    }

    Directory<T> *dir;
    int seg_num;
    int restart;
#ifdef PMEM
    PMEMobjpool *pool_addr;
#endif
};
#endif  // EXTENDIBLE_PTR_H_

template<class T>
int Segment<T>::Insert(PMEMobjpool *pool_addr, T key, Value_t value, size_t loc, size_t key_hash) {
#ifdef INPLACE
  if (sema == -1) {
    return 2;
  };
  get_lock(pool_addr);
  if ((key_hash >> (8*sizeof(key_hash)-local_depth)) != pattern || sema == -1){
    release_lock(pool_addr);
    return 2;
  } 
  int ret = 1;
  T LOCK = (T)INVALID;

  /*uniqueness check*/
  auto slot = loc;
  for (unsigned i = 0; i < kNumCacheLine*kNumPairPerCacheLine; ++i)
  {
    slot = (loc+i)%kNumSlot;
    if constexpr (std::is_pointer_v<T>){
      if (_[slot].key != (T)INVALID && (strcmp(key, _[slot].key) == 0))
      {
        release_lock(pool_addr);
        return 0;
      }
    }else{
      if (_[slot].key == key)
      {
        release_lock(pool_addr);
        return 0;
      }
    }
  }

  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    slot = (loc + i) % kNumSlot;
    if constexpr (std::is_pointer_v<T>){
      if ((_[slot].key != (T)INVALID) && ((h(_[slot].key,strlen(_[slot].key)) >> (8*sizeof(key_hash)-local_depth)) != pattern)) {
        _[slot].key = (T)INVALID;
      }
      if (CAS(&_[slot].key, &LOCK, SENTINEL)) {
        //T new_key = (T)malloc(strlen(key)+1);
        //strcpy(new_key, key);
        //clflush((char*)new_key, strlen(key)+1);
        _[slot].value = value;
        mfence();       
        _[slot].key = key;
        //clflush((char*)&_[slot],sizeof(_Pair<T>));
#ifdef PMEM
        Allocator::Persist(&_[slot], sizeof(_Pair<T>));
#endif
        count++;
        ret = 0;
        break;
      } else {
        LOCK = (T)INVALID;
      }
    }else{
      if ((h(&_[slot].key,sizeof(Key_t)) >> (8*sizeof(key_hash)-local_depth)) != pattern) {
        _[slot].key = INVALID;
      }
      if (CAS(&_[slot].key, &LOCK, SENTINEL)) {
        _[slot].value = value;
        mfence();
        _[slot].key = key;
        //clflush((char*)&_[slot],sizeof(_Pair<T>));
        Allocator::Persist(&_[slot], sizeof(_Pair<T>));
        count++;
        ret = 0;
        break;
      } else {
        LOCK = INVALID;
      }
    }
  }
  release_lock(pool_addr);
  return ret;
#else
  if (sema == -1) return 2;
  if ((key_hash >> (8*sizeof(key_hash)-local_depth)) != pattern) return 2;
  auto lock = sema;
  int ret = 1;
  while (!CAS(&sema, &lock, lock+1)) {
    lock = sema;
  }
  Key_t LOCK = INVALID;
  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (loc + i) % kNumSlot;
    if (CAS(&_[slot].key, &LOCK, SENTINEL)) {
      _[slot].value = value;
      mfence();
      _[slot].key = key;
      clflush((char*)&_[slot],sizeof(_Pair<T>));
      ret = 0;
      break;
    } else {
      LOCK = INVALID;
    }
  }
  lock = sema;
  while (!CAS(&sema, &lock, lock-1)) {
    lock = sema;
  }
  return ret;
#endif
}

template<class T>
int Segment<T>::Insert4split(T key, Value_t value, size_t loc) {
  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (loc+i) % kNumSlot;
    if (_[slot].key == (T)INVALID) {
      _[slot].key = key;
      _[slot].value = value;
      return 0;
    }
  }
  return -1;
}

template<class T>
Segment<T>* Segment<T>::Split(PMEMobjpool *pool_addr){
  using namespace std;
#ifndef INPLACE
  int64_t lock = 0;
  if (!CAS(&sema, &lock, -1)) return nullptr;
#else
  if(!try_get_lock(pool_addr))
  {
    return nullptr;
  }
  sema = -1;
#endif

#ifdef INPLACE
  /* How to get the pool handler?*/
  size_t new_pattern = (pattern << 1) + 1;
  size_t old_pattern = pattern << 1;
  //Segment<T> *split = new Segment<T>(local_depth+1);
  Segment<T> *split;
  Segment::New(&split, local_depth + 1);

  for (unsigned i = 0; i < kNumSlot; ++i) {
    uint64_t key_hash;
    if constexpr (std::is_pointer_v<T>){
      if (_[i].key != (T)INVALID)
       {
         key_hash = h(_[i].key, strlen(_[i].key));
       } 
    }else{
      key_hash = h(&_[i].key, sizeof(Key_t));
    }
    if ((_[i].key != (T)INVALID) && (key_hash >> (8*8-local_depth-1) == new_pattern))
    {
        split->Insert4split
        (_[i].key, _[i].value, (key_hash & kMask)*kNumPairPerCacheLine);
        _[i].key = (T)INVALID;
        split->count++;
        count--;
    }
  }

  //clflush((char*)split, sizeof(Segment<T>));
#ifdef PMEM
  Allocator::Persist(split, sizeof(Segment<T>));
#endif
  if constexpr (std::is_pointer_v<T>){
    //clflush((char*)this, sizeof(Segment<T>));
#ifdef PMEM
    Allocator::Persist(this,sizeof(Segment<T>));
#endif
  }
  return split;
#else
  Segment** split = new Segment*[2];
  split[0] = new Segment(local_depth+1);
  split[1] = new Segment(local_depth+1);

  for (unsigned i = 0; i < kNumSlot; ++i) {
    auto key_hash = h(&_[i].key, sizeof(Key_t));
    if (key_hash & ((size_t) 1 << ((sizeof(Key_t)*8 - local_depth - 1)))) {
      split[1]->Insert4split
        (_[i].key, _[i].value, (key_hash & kMask)*kNumPairPerCacheLine);
    } else {
      split[0]->Insert4split
        (_[i].key, _[i].value, (key_hash & kMask)*kNumPairPerCacheLine);
    }
  }

  clflush((char*)split[0], sizeof(Segment));
  clflush((char*)split[1], sizeof(Segment));

  return split;
#endif
}

template<class T>
CCEH<T>::CCEH(int initCap)
{
  //dir = new Directory<T>(initCap);
  Directory<T>::New(&dir, initCap); 
  auto dir_entry = dir->sa->_;
  for (int i = 0; i < dir->capacity; ++i)
  {
    //dir_entry[i] = new Segment<T>(dir->sa->global_depth);
    Segment<T>::New(&dir_entry[i], dir->sa->global_depth);
    dir_entry[i]->pattern = i;
  }
  seg_num = 0;
  restart = 0;
}

template<class T>
CCEH<T>::~CCEH(void)
{ }

template<class T>
void CCEH<T>::Directory_Doubling(int x, Segment<T> *s0, Segment<T> *s1){
  Seg_array<T> *sa = dir->sa;
  Segment<T>** d = sa->_;
  auto global_depth = sa->global_depth;
  //printf("Directory_Doubling towards %lld\n", global_depth+1);

  /* new segment array*/
  //auto new_sa = new Seg_array<T>(2*dir->capacity);
  Seg_array<T>::New(&dir->new_sa, 2*dir->capacity);
  auto dd = dir->new_sa->_;

  for (unsigned i = 0; i < dir->capacity; ++i) {
    if (i == x) {
      dd[2*i] = s0;
      dd[2*i+1] = s1;
    } else {
      dd[2*i] = d[i];
      dd[2*i+1] = d[i];
    }
  }

#ifdef PMEM
  Allocator::Persist(dd, sizeof(Segment<T>*)*2*dir->capacity);
  Allocator::Persist(dir->new_sa, sizeof(Seg_array<T>));
#endif
  dir->sa = dir->new_sa;
  dir->new_sa = nullptr;
#ifdef PMEM
  Allocator::Persist(&dir->sa, sizeof(dir->sa));
  Allocator::Persist(&dir->new_sa, sizeof(dir->new_sa));
#endif
  dir->capacity *= 2;
#ifdef PMEM
  Allocator::Persist(&dir->capacity, sizeof(dir->capacity));
#endif
  //delete [] d;
  //delete sa;
printf("Done!!Directory_Doubling towards %lld\n", global_depth);
}

template<class T>
void CCEH<T>::Lock_Directory(){  
  while (!dir->Acquire()) {
      asm("nop");
  }
}

template<class T>
void CCEH<T>::Unlock_Directory(){
  while (!dir->Release()) {
      asm("nop");
  }
}

template<class T>
void CCEH<T>::Directory_Update(int x, Segment<T> *s0, Segment<T> *s1){
//  printf("directory update for %d\n", x);
  Segment<T>** dir_entry = dir->sa->_;
  auto global_depth = dir->sa->global_depth;
  unsigned depth_diff = global_depth - s0->local_depth;
  if (depth_diff == 0) {
    if (x%2 == 0) {
      dir_entry[x+1] = s1;
#ifdef PMEM
      Allocator::Persist(&dir_entry[x+1], sizeof(Segment<T>*));
#endif
    } else {
      dir_entry[x] = s1;
#ifdef PMEM
      Allocator::Persist(&dir_entry[x], sizeof(Segment<T>*));
#endif
    }
  } else {
    int chunk_size = pow(2, global_depth - (s0->local_depth - 1));
    x = x - (x % chunk_size);
    for (unsigned i = 0; i < chunk_size/2; ++i) {
      dir_entry[x+chunk_size/2+i] = s1;
    }
#ifdef PMEM
      Allocator::Persist(&dir_entry[x + chunk_size / 2],
                         sizeof(Segment<T>*) * chunk_size / 2);
#endif
  }
//  printf("Done!directory update for %d\n", x);
}

template<class T>
void CCEH<T>::Insert(T key, Value_t value) {
STARTOVER:
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>){
    key_hash = h(key, strlen(key));
  }else{
    key_hash = h(&key, sizeof(key));
  }
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment<T>* target = dir_entry[x];
  if (old_sa != dir->sa)
  {
    goto RETRY;
  }

  auto ret = target->Insert(pool_addr, key, value, y, key_hash);
  //std::cout<<"I am inserting "<<key<<", the x = "<<x<<std::endl;

  if (ret == 1) {
    Segment<T>* s = target->Split(pool_addr);
    if (s == nullptr) {
      goto RETRY;
    }
    //printf("Done!!Segment split for key %lld\n", key);
    target->local_depth += 1;
    //clflush((char*)&target->local_depth, sizeof(target->local_depth));
#ifdef PMEM
    Allocator::Persist(&target->local_depth,
                         sizeof(target->local_depth));
#endif

    target->pattern = (key_hash >> (8*sizeof(key_hash)-target->local_depth+1)) << 1;
    s->pattern = ((key_hash >> (8*sizeof(key_hash)-s->local_depth+1)) << 1) + 1;

    // Directory management
    Lock_Directory();
    { // CRITICAL SECTION - directory update
      auto sa = dir->sa;
      dir_entry = sa->_;

      x = (key_hash >> (8*sizeof(key_hash)-sa->global_depth));
      target = dir_entry[x];
#ifdef INPLACE
      if (target->local_depth-1<sa->global_depth)
      {
#else
      if(target->local_depth< sa->global_depth){
#endif
        Directory_Update(x,target,s);
      } else {  // directory doubling
        Directory_Doubling(x,target,s);
      }
#ifdef INPLACE
      target->sema = 0;
      target->release_lock(pool_addr);
#endif
    }  // End of critical section
    Unlock_Directory();
    goto RETRY;
  } else if (ret == 2) {
    goto STARTOVER;
  } 
}

// TODO
template<class T>
bool CCEH<T>::Delete(T key) {
  return false;
}

template<class T>
Value_t CCEH<T>::Get(T key) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>){
    key_hash = h(key, strlen(key));
  }else{
    key_hash = h(&key, sizeof(key));
  }
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment<T>* dir_ = dir_entry[x];

#ifdef INPLACE
  auto sema = dir_->sema;
  if (sema == -1)
  {
    goto RETRY;
  }
  //dir_->mutex.lock_shared();
  dir_->get_rd_lock(pool_addr);

  if ((key_hash >> (8*sizeof(key_hash)-dir_->local_depth)) != dir_->pattern || dir_->sema == -1){
    //dir_->mutex.unlock_shared();
    dir_->release_rd_lock(pool_addr);
    goto RETRY;
  } 
#endif

  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (y+i) % Segment<T>::kNumSlot;
    if constexpr (std::is_pointer_v<T>){
      if ((dir_->_[slot].key != (T)INVALID) && (strcmp(dir_->_[slot].key, key) == 0))
      {
         auto value = dir_->_[slot].value;
  #ifdef INPLACE
        //dir_->mutex.unlock_shared();
        dir_->release_rd_lock(pool_addr);
  #endif
        return value;
      }
    }else{
      if (dir_->_[slot].key == key) {
        auto value = dir_->_[slot].value;
  #ifdef INPLACE
        //dir_->mutex.unlock_shared();
        dir_->release_rd_lock(pool_addr);
  #endif
        return value;
      }
    }      
  }

#ifdef INPLACE
  //dir_->mutex.unlock_shared();
  dir_->release_rd_lock(pool_addr);
#endif
  return NONE;
}
/*
double CCEH::Utilization(void) {
  size_t sum = 0;
  std::unordered_map<Segment*, bool> set;
  for (size_t i = 0; i < dir.capacity; ++i) {
    set[dir._[i]] = true;
  }
  for (auto& elem: set) {
    for (unsigned i = 0; i < Segment::kNumSlot; ++i) {
      if (elem.first->_[i].key != INVALID) sum++;
    }
  }
  return ((double)sum)/((double)set.size()*Segment::kNumSlot)*100.0;
}

size_t CCEH::Capacity(void) {
  std::unordered_map<Segment*, bool> set;
  for (size_t i = 0; i < dir.capacity; ++i) {
    set[dir._[i]] = true;
  }
  return set.size() * Segment::kNumSlot;
}

size_t Segment::numElem(void) {
  size_t sum = 0;
  for (unsigned i = 0; i < kNumSlot; ++i) {
    if (_[i].key != INVALID) {
      sum++;
    }
  }
  return sum;
}

bool CCEH::Recovery(void) {
  bool recovered = false;
  size_t i = 0;
  while (i < dir.capacity) {
    size_t depth_cur = dir._[i]->local_depth;
    size_t stride = pow(2, global_depth - depth_cur);
    size_t buddy = i + stride;
    if (buddy == dir.capacity) break;
    for (int j = buddy - 1; i < j; j--) {
      if (dir._[j]->local_depth != depth_cur) {
        dir._[j] = dir._[i];
      }
    }
    i = i+stride;
  }
  if (recovered) {
    clflush((char*)&dir._[0], sizeof(void*)*dir.capacity);
  }
  return recovered;
}

// for debugging
Value_t CCEH::FindAnyway(Key_t& key) {
  using namespace std;
  for (size_t i = 0; i < dir.capacity; ++i) {
     for (size_t j = 0; j < Segment::kNumSlot; ++j) {
       if (dir._[i]->_[j].key == key) {
         auto key_hash = h(&key, sizeof(key));
         auto x = (key_hash >> (8*sizeof(key_hash)-global_depth));
         auto y = (key_hash & kMask) * kNumPairPerCacheLine;
         cout << bitset<32>(i) << endl << bitset<32>((x>>1)) << endl << bitset<32>(x) << endl;
         return dir._[i]->_[j].value;
       }
     }
  }
  return NONE;
}

void Directory::SanityCheck(void* addr) {
  using namespace std;
  for (unsigned i = 0; i < capacity; ++i) {
    if (_[i] == addr) {
      cout << i << " " << _[i]->sema << endl;
      exit(1);
    }
  }
}*/
