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
#include <libpmemobj.h>

const size_t kCacheLineSize = 64;
constexpr size_t kSegmentBits = 8;
constexpr size_t kMask = (1 << kSegmentBits)-1;
constexpr size_t kShift = kSegmentBits;
constexpr size_t kSegmentSize = (1 << kSegmentBits) * 16 * 4;
constexpr size_t kNumPairPerCacheLine = 4;
constexpr size_t kNumCacheLine = kCacheLineSize/sizeof(Pair);
constexpr size_t LogMask = (1 << 10)-1;
//constexpr size_t kNumSlot = kSegmentSize/sizeof(Pair);
#define INPLACE 1
//#define PERSISTENT_LOCK 1

#define SEGMENT_TYPE 1000
#define DIRECTORY_TYPE 2000
#define ARRAY_TYPE 3000
#define CCEH_TYPE 4000
#define LOG_TYPE 5000

struct Segment {
  static const size_t kNumSlot = kSegmentSize/sizeof(Pair);
  Segment(void)
  { }

  Segment(size_t depth)
  { }

  ~Segment(void) {}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

  int Insert(PMEMobjpool *pop,Key_t&, Value_t, size_t&, size_t);
  void Insert4split(Key_t&, Value_t, size_t);
  bool Put(Key_t&, Value_t, size_t);
  Segment** Split(PMEMobjpool *pop);
  
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

  char dummy[48];//To make the segment cacheline-aligned
  Pair _[kNumSlot];
  size_t local_depth;
  int64_t sema = 0;
  size_t pattern = 0;
  int count;
  PMEMrwlock rwlock;
  //std::shared_mutex mutex;

  size_t numElem(void); 
};

struct depth_pattern{
  size_t local_depth;
  size_t pattern;
  size_t sema;
};

static int create_segment(PMEMobjpool *pop, void *ptr, void *arg){
  struct Segment *se = (struct Segment *)ptr;
  struct depth_pattern *dp = (struct depth_pattern*)arg;

  memset(&se->_[0],255,sizeof(Pair)*Segment::kNumSlot);
  se->count = 0;
  se->local_depth = dp->local_depth;
  se->pattern = dp->pattern;
  se->sema = dp->sema;
  memset(&se->rwlock, 0, sizeof(PMEMrwlock));
  pmemobj_persist(pop, se, sizeof(*se));
  return 0;
}

struct Seg_array{
  //PMEMoid _;
  uint64_t _;/* the offset of the directory entry*/
  size_t global_depth;
};

static int create_seg_array(PMEMobjpool *pop, void *ptr, void *arg){
  struct Seg_array *sa = (struct Seg_array *)ptr;
  //sa->_ = OID_NULL;
  sa->_ = 0;
  sa->global_depth = *((size_t*)arg);
  pmemobj_persist(pop, sa, sizeof(*sa));
  return 0;
}

struct Directory {
  static const size_t kDefaultDirectorySize = 1024;
  //PMEMoid sa;
  //PMEMoid new_sa;
  Seg_array* sa; // Current directory
  uint64_t new_sa; //The new directory(offset)
  size_t capacity;
  bool lock;
  int sema = 0;

  Directory(void) {
  }

  Directory(size_t size) {
  }

  ~Directory(void) {
    
  }

  /*
  void get_item_num(uint64_t base_addr){
    size_t count = 0;
    Seg_array *seg = (Seg_array*)pmemobj_direct(sa);
    uint64_t *dir_entry = (uint64_t*)(base_addr+seg->_);
    Segment *ss;
    auto global_depth = seg->global_depth;
    size_t depth_diff;
    for (int i = 0; i < capacity;)
    {
      ss = (Segment*)(base_addr+dir_entry[i]);
      depth_diff = global_depth - ss->local_depth;
      count += ss->count;
      i += pow(2, depth_diff);
    }
    printf("the item stored in this hash table is %lld\n", count);
  }
  */

  bool Acquire(void) {
    bool unlocked = false;
    return CAS(&lock, &unlocked, true);
  }
  
  bool Release(void) {
    bool locked = true;
    return CAS(&lock, &locked, false);
  }
  
  void SanityCheck(void*);
  void LSBUpdate(int, int, int, int, Segment**);
};

/*
inline void atomic_pointer_update(PMEMoid& a, PMEMoid& b){  
//  assert(a.pool_uuid_lo == b.pool_uuid_lo);
  a.off = b.off; 
}*/

/* arg is the pointer to the capacity*/
static int create_directory(PMEMobjpool *pop, void *ptr, void *arg){
  struct Directory *dd = (struct Directory *)ptr;
  dd->capacity = *((size_t *)arg);
  dd->lock = false;
  dd->sema = 0;
  dd->sa = 0; 
  dd->new_sa = 0;
  pmemobj_persist(pop, dd, sizeof(*dd));
  return 0;
}

struct log_entry{
  uint64_t lock;
  PMEMoid temp;
};

class CCEH{
  public:
    CCEH(void);
    CCEH(size_t);
    ~CCEH(void);
    void Insert(PMEMobjpool *pop, Key_t& key, Value_t value);
    bool InsertOnly(PMEMobjpool *pop, Key_t&, Value_t);
    bool Delete(Key_t&);
    Value_t Get(PMEMobjpool* pop, Key_t&);
    Value_t FindAnyway(Key_t&);
    double Utilization(void);
    size_t Capacity(void);
    bool Recovery(void);
    void Directory_Doubling(PMEMobjpool *pop, Seg_array *sa, uint64_t *dir_entry, size_t global_depth, size_t x, size_t new_pattern, PMEMoid*);
    void Directory_Update(PMEMobjpool *pop, uint64_t *dir_entry, size_t global_depth, size_t x, size_t new_pattern, PMEMoid*);
    void Lock_Directory();
    void Unlock_Directory();
    int Segment_Split(PMEMobjpool *pop, size_t key_hash, Segment *_target);
    void TX_Swap(PMEMobjpool *pop, uint64_t& entry, PMEMoid *new_seg);
    void Get_Number(){
      //dir->get_item_num(base_addr);
    }
    bool Lock_log(log_entry *le){
      uint64_t temp = 0;
      return CAS(&le->lock, &temp, 1);
    }
    bool Unlock_log(log_entry *le){
      uint64_t temp = 1;
      return CAS(&le->lock, &temp, 0);
    }

    void* operator new(size_t size) {
      void *ret;
      posix_memalign(&ret, 64, size);
      return ret;
    }

    PMEMoid _dir;// the pointer to the direcotry;
    Directory *dir;
    PMEMoid _log;
    log_entry *log;
    uint64_t base_addr;
    int seg_num;
    int restart;
    int log_num;
};

void remapping(PMEMobjpool *pop, CCEH *eh){
  TX_BEGIN(pop){
    pmemobj_tx_add_range_direct(eh, sizeof(*eh));
    auto old_base = eh->base_addr;
    eh->base_addr = (uint64_t)pop;
    eh->dir = (Directory*)pmemobj_direct(eh->_dir);
    eh->log = (log_entry*)pmemobj_direct(eh->_log);
    pmemobj_tx_add_range_direct(eh->dir, sizeof(struct Directory));
    eh->dir->sa = reinterpret_cast<Seg_array*>(reinterpret_cast<uint64_t>(eh->dir->sa) - old_base + eh->base_addr); 
  }TX_ONABORT{
    printf("remapping failure\n");
  }TX_END
}

static int create_CCEH(PMEMobjpool *pop, void *ptr, void *arg){
  CCEH *eh = (CCEH *)ptr;
  eh->_dir = OID_NULL;
  eh->dir = nullptr;
  eh->_log = OID_NULL;
  eh->log = nullptr;
  eh->base_addr = 0;
  eh->seg_num = 0;
  eh->restart = 0;
  return 0;
}

void segment_allocate(PMEMobjpool* pop, uint64_t& offset, depth_pattern *dp){
  PMEMoid temp = OID_NULL;
  TX_BEGIN(pop){
    pmemobj_tx_add_range_direct(&offset, sizeof(offset));
    pmemobj_alloc(pop,&temp, sizeof(struct Segment), SEGMENT_TYPE, create_segment, dp);
    offset = temp.off;
  }TX_ONABORT{
    pmemobj_free(&temp);
    printf("segment allocation failure\n");
  }TX_END
}

void Initialize_CCEH(PMEMobjpool *pop, CCEH *eh, size_t _capacity){
  //printf("start the initialization for CCEH\n");
  eh->base_addr = (uint64_t)pop;
  pmemobj_alloc(pop, &eh->_dir, sizeof(struct Directory), DIRECTORY_TYPE, create_directory, &_capacity);

  /* Initialize the segment_array, the pair_array and the segment pointed by the pair_array*/
    /* Segment array allocation*/
  Directory *dd = (Directory*)pmemobj_direct(eh->_dir);
  eh->dir = dd;

  PMEMoid temp = OID_NULL;
  TX_BEGIN(pop){
    temp = pmemobj_tx_zalloc(sizeof(struct Seg_array), ARRAY_TYPE);
    //dd->sa = temp.off;
    dd->sa = reinterpret_cast<Seg_array*>(temp.off + eh->base_addr);
  }TX_ONABORT{
    printf("allocation error duing seg_array initialization\n");
  }TX_END

  //pmemobj_alloc(pop, &dd->sa, sizeof(Seg_array), ARRAY_TYPE, NULL, NULL);
  //Seg_array *sa = (Seg_array*)pmemobj_direct(dd->sa);
  Seg_array *sa = dd->sa;
  sa->global_depth = static_cast<size_t>(log2(dd->capacity));
  //pmemobj_alloc(pop, &sa->_, sizeof(PMEMoid)*dd->capacity, SEGMENT_TYPE+1, NULL, NULL);
  //PMEMoid* dir_entry = (PMEMoid*)pmemobj_direct(sa->_);
  temp = OID_NULL;
  TX_BEGIN(pop){
    temp = pmemobj_tx_zalloc(sizeof(uint64_t)*dd->capacity, SEGMENT_TYPE+1);
    sa->_ = temp.off;
  }TX_ONABORT{
    printf("allocation error for direcoty entry during initialization\n");
  }TX_END

  uint64_t* dir_entry = (uint64_t*)(eh->base_addr+sa->_);

  struct depth_pattern dp;
  dp.local_depth = sa->global_depth;
  dp.sema = 0;

  /* Segment allocation*/
  for (int i = 0; i < dd->capacity; ++i)
  {
    dp.pattern = i;
    segment_allocate(pop, dir_entry[i], &dp);
  }

  /* Create and Initialize log*/
  eh->log_num = 1024;
  pmemobj_zalloc(pop, &eh->_log, sizeof(struct log_entry)*eh->log_num, LOG_TYPE);
  eh->log = (log_entry*)pmemobj_direct(eh->_log);

  /* flush the new created semgnets and the direcoty entry array*/
  pmemobj_persist(pop, sa, sizeof(struct Seg_array));
  pmemobj_persist(pop, dd, sizeof(struct Directory));
  pmemobj_persist(pop, eh, sizeof(struct CCEH));
  //printf("end of the initialization of CCEH\n");
}

#endif  // EXTENDIBLE_PTR_H_

int Segment::Insert(PMEMobjpool *pop, Key_t& key, Value_t value, size_t& loc, size_t key_hash) {
#ifdef INPLACE
  if (sema == -1) {
    //printf("retrying because the sema still equals -1\n");
    return 2;
  };
  get_lock(pop);
  if ((key_hash >> (8*sizeof(key_hash)-local_depth)) != pattern || sema == -1){
  //  std::cout<<"The computing pattern is "<<(key_hash >> (8*sizeof(key_hash)-local_depth))<<std::endl;
   // std::cout<<"the pattern is "<<pattern<<std::endl;
    release_lock(pop);
    return 2;
  }

    /*uniqueness check*/
  auto slot = loc;
  for (unsigned i = 0; i < kNumCacheLine*kNumPairPerCacheLine; ++i)
  {
    slot = (loc+i)%kNumSlot;
    if (_[slot].key == key)
    {
      release_lock(pop);
      return 0;
    }
  }

  int ret = 1;
  Key_t LOCK = INVALID;
  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    slot = (loc + i) % kNumSlot;
    if ((h(&_[slot].key,sizeof(Key_t)) >> (8*sizeof(key_hash)-local_depth)) != pattern) {
      _[slot].key = INVALID;
    }
    if (CAS(&_[slot].key, &LOCK, SENTINEL)) {
      _[slot].value = value;
      mfence();
      _[slot].key = key;
      pmemobj_persist(pop, &_[slot], sizeof(Pair));
      count++;
      ret = 0;
      break;
    } else {
      LOCK = INVALID;
    }
  }
  loc = slot;
  //mutex.unlock();
  release_lock(pop);
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
      clflush((char*)&_[slot],sizeof(Pair));
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

void Segment::Insert4split(Key_t& key, Value_t value, size_t loc) {
  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (loc+i) % kNumSlot;
    if (_[slot].key == INVALID) {
      _[slot].key = key;
      _[slot].value = value;
      return;
    }
  }
}

int CCEH::Segment_Split(PMEMobjpool *pop, size_t key_hash, Segment *_target){
  if(!_target->try_get_lock(pop))
  {
    return -1;/* it represents that this segment is not available now, it needs to retry*/
  }
  _target->sema = -1;

  /* Memory leak prevention*/
  size_t new_pattern = (_target->pattern << 1) + 1;
  size_t old_pattern = _target->pattern << 1;
  //printf("segment split from %lld to %lld and %lld\n", _target->pattern, old_pattern, new_pattern);

  /* first allocate a segment outside the direcotry*/
  struct depth_pattern dp;
  dp.sema = -1;
  dp.local_depth = _target->local_depth + 1;
  dp.pattern = new_pattern;

  uint64_t log_pos = (key_hash & LogMask);
  //printf("the log position is %lld\n", log_pos);
  while(!Lock_log(&log[log_pos])){
    asm("nop");
  }
  pmemobj_alloc(pop, &log[log_pos].temp, sizeof(struct Segment), SEGMENT_TYPE, create_segment, &dp);  

  /* grab the lock in the directory, the memory allocator may be bottleneck*/
  Lock_Directory();
  //auto sa = (Seg_array*)pmemobj_direct(dir->sa);
  //Seg_array* sa = reinterpret_cast<Seg_array*>(dir->sa + base_addr);
  Seg_array *sa = dir->sa;
  auto dir_entry = (uint64_t*)(base_addr+sa->_);

  size_t x = (key_hash >> (8*sizeof(key_hash)-sa->global_depth));
  auto target = (Segment*)(base_addr+dir_entry[x]);
  assert(target == _target);

  if (target->local_depth<sa->global_depth)
  {
    Directory_Update(pop,dir_entry, sa->global_depth, x, new_pattern, &log[log_pos].temp);  
  } else {  // directory doubling
    Directory_Doubling(pop, sa, dir_entry, sa->global_depth, x, new_pattern, &log[log_pos].temp);
  }

  sa = reinterpret_cast<Seg_array*>(dir->sa);
  dir_entry = (uint64_t*)(base_addr+sa->_);

  x = (key_hash >> (8*sizeof(key_hash)-sa->global_depth));
  /* compute the brother*/
  int chunk_size = pow(2, sa->global_depth - (target->local_depth));
  x = x - (x % chunk_size);
  //target = (Segment*)pmemobj_direct(dir_entry[x]);
  target = (Segment*)(base_addr+dir_entry[x]);
  //auto brother = (Segment*)pmemobj_direct(dir_entry[x+chunk_size/2]);/* Find the wrong brother bucket...*/
  auto brother = (Segment*)(base_addr+dir_entry[x+chunk_size/2]);
  //assert(brother->sema == -1);
  //assert(target == _target);
  Unlock_Directory();
  assert(OID_IS_NULL(log[log_pos].temp));
  while(!Unlock_log(&log[log_pos])){
    asm("nop");
  }

  auto _ = target->_;
  /* Initialize the plit segment, actually for the new segment, its pattern and local depth has been updated, but sema and k-v need to be udpdate*/
  for (unsigned i = 0; i < Segment::kNumSlot; ++i) {
    auto key_hash = h(&_[i].key, sizeof(Key_t));
    if (_[i].key != INVALID && key_hash >> (8*8-target->local_depth-1) == new_pattern)
    {
        brother->Insert4split
        (_[i].key, _[i].value, (key_hash & kMask)*kNumPairPerCacheLine);
        brother->count++;
        target->count--;
    }
  }

  pmemobj_persist(pop, brother, sizeof(*brother));

  TX_BEGIN(pop){
    pmemobj_tx_add_range_direct(&target->local_depth, sizeof(target->local_depth));
    pmemobj_tx_add_range_direct(&target->pattern, sizeof(target->pattern));
    target->local_depth = target->local_depth+1;
    target->pattern = old_pattern;
  }TX_END

  brother->sema = 0;
  target->sema = 0;
  target->release_lock(pop);
  return 0;
}

CCEH::CCEH(void)
{
}

CCEH::CCEH(size_t initCap)
{
}

CCEH::~CCEH(void)
{ }

void CCEH::Directory_Doubling(PMEMobjpool *pop, Seg_array *sa, uint64_t* dir_entry, size_t global_depth, size_t x, size_t new_pattern, PMEMoid* new_seg){
  uint64_t* d = dir_entry;
  //auto target = (Segment*)pmemobj_direct(d[x]);
  auto target = (Segment*)(base_addr+d[x]);
  printf("Directory_Doubling towards %lld\n", global_depth+1);

  PMEMoid temp = OID_NULL;
  TX_BEGIN(pop){
    temp = pmemobj_tx_zalloc(sizeof(struct Seg_array), ARRAY_TYPE);
    dir->new_sa = temp.off;
  }TX_ONABORT{
    printf("allocation error duing seg_array initialization\n");
  }TX_END

  Seg_array *new_sa = reinterpret_cast<Seg_array*>(base_addr + dir->new_sa);
  new_sa->global_depth = global_depth + 1;
  //pmemobj_alloc(pop, &dir->new_sa, sizeof(Seg_array), ARRAY_TYPE, create_seg_array, &new_global_depth);
  //Seg_array *new_sa = (Seg_array*)pmemobj_direct(dir->new_sa);
  //pmemobj_zalloc(pop, &new_sa->_, sizeof(PMEMoid)*2*dir->capacity, SEGMENT_TYPE+1);
  //PMEMoid *dd = (PMEMoid*)pmemobj_direct(new_sa->_);

  temp = OID_NULL;
  TX_BEGIN(pop){
    temp = pmemobj_tx_zalloc(sizeof(uint64_t)*dir->capacity*2, SEGMENT_TYPE+1);
    new_sa->_ = temp.off;
  }TX_ONABORT{
    printf("allocation error for direcoty entry during initialization\n");
  }TX_END

  uint64_t* dd = (uint64_t*)(base_addr+new_sa->_);

  /* attach the new segment to the new_sa*/
  for (unsigned i = 0; i < dir->capacity; ++i) {  
      dd[2*i] = d[i];
      dd[2*i+1] = d[i];
  }

  struct depth_pattern dp;
  dp.sema = -1;
  dp.local_depth = target->local_depth + 1;
  dp.pattern = new_pattern;

  /* allocate the new segment*/
  TX_Swap(pop, dd[2*x+1], new_seg);
  pmemobj_persist(pop, dd, sizeof(PMEMoid)*2*dir->capacity);
  pmemobj_persist(pop, new_sa, sizeof(struct Seg_array));

  //auto old_sa = pmemobj_oid((void*)(base_addr+dir->sa));
  auto old_sa = pmemobj_oid((void*)(dir->sa));
  auto old_dd = pmemobj_oid((void*)(base_addr+sa->_));
  // use txn to make this atomic: capacity of the directory and the segment array
  TX_BEGIN(pop){
    pmemobj_tx_add_range(_dir, 0, sizeof(struct Directory));
    dir->capacity = 2*dir->capacity;
    dir->sa = reinterpret_cast<Seg_array*>(dir->new_sa + base_addr);
    dir->new_sa = 0;
    pmemobj_free(&old_dd);
    pmemobj_free(&old_sa);
  }TX_END
 // printf("Done!!Directory_Doubling towards %lld\n", global_depth);
}

void CCEH::Lock_Directory(){  
  //auto my_dir = (Directory*)pmemobj_direct(dir);
  while (!dir->Acquire()) {
      asm("nop");
  }
}

void CCEH::Unlock_Directory(){
  //auto my_dir = (Directory*)pmemobj_direct(dir);
  while (!dir->Release()) {
      asm("nop");
  }
}

void CCEH::TX_Swap(PMEMobjpool *pop, uint64_t& entry, PMEMoid *new_seg){
  TX_BEGIN(pop){
    pmemobj_tx_add_range_direct(&entry, sizeof(uint64_t));
    pmemobj_tx_add_range_direct(new_seg,sizeof(*new_seg));
    entry = new_seg->off;
    *new_seg = OID_NULL;
  }TX_ONABORT{
    printf("error!!! on txn_swap\n");
  }TX_END
}

void CCEH::Directory_Update(PMEMobjpool *pop, uint64_t *dir_entry, size_t global_depth, size_t x, size_t new_pattern, PMEMoid *new_seg){
  auto target = (Segment*)(base_addr+dir_entry[x]);

  unsigned depth_diff = global_depth - target->local_depth;
  struct depth_pattern dp;
  dp.sema = -1;
  dp.local_depth = target->local_depth + 1;
  dp.pattern = new_pattern;
  
  if (depth_diff == 1) {
    if (x%2 == 0) {
      TX_Swap(pop, dir_entry[x+1], new_seg);
#ifdef INPLACE
#else
      TX_Swap(pop, dir_entry[x], new_seg); 
#endif
    } else {
      TX_Swap(pop, dir_entry[x], new_seg); 
#ifdef INPLACE
#else
      TX_Swap(pop, dir_entry[x-1], new_seg); 
#endif
    }
  } else {
    int chunk_size = pow(2, global_depth - (target->local_depth));
    x = x - (x % chunk_size);
    /* first we need to do the atomic attach*/
    TX_Swap(pop, dir_entry[x+chunk_size/2], new_seg); 
    auto s1 = dir_entry[x+chunk_size/2]; 
    for (unsigned i = 1; i < chunk_size/2; ++i) {
      dir_entry[x+chunk_size/2+i] = s1;
    }
    pmemobj_persist(pop, &dir_entry[x+chunk_size/2], sizeof(uint64_t)*chunk_size/2);
#ifndef INPLACE
    TX_Swap(pop, dir_entry[x], new_seg);/* this dp may be uncorrect*/
    auto s0 = dir_entry[x];
    for (unsigned i = 1; i < chunk_size/2; ++i) {
      dir_entry[x+i] = s0;
    }
    pmemobj_persist(pop, &dir_entry[x], sizeof(uint64_t)*chunk_size/2);
#endif
  }
}

void CCEH::Insert(PMEMobjpool *pop, Key_t& key, Value_t value) {
STARTOVER:
  auto key_hash = h(&key, sizeof(key));
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  Seg_array *sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-sa->global_depth));
  auto dir_entry = (uint64_t*)(base_addr+sa->_);
  Segment* target = (Segment*)(dir_entry[x]+base_addr);
  auto yy = y;
  auto ret = target->Insert(pop, key, value, yy, key_hash);
  //std::cout<<"I am inserting "<<key<<", the x = "<<x<<std::endl;

  if (ret == 1) {
    //printf("Segment split for key %lld\n", key);
    Segment_Split(pop, key_hash, target);
    goto RETRY;
  } else if (ret == 2) {
    goto STARTOVER;
  } 
}

// TODO
bool CCEH::Delete(Key_t& key) {
  return false;
}

Value_t CCEH::Get(PMEMobjpool*pop, Key_t& key) {
RETRY:
  auto key_hash = h(&key, sizeof(key));
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

  Seg_array *sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-sa->global_depth));
  auto dir_entry = (uint64_t*)(base_addr+sa->_);
  Segment* dir_ = (Segment*)(dir_entry[x]+base_addr);

#ifdef INPLACE
  auto sema = dir_->sema;
  if (sema == -1)
  {
    goto RETRY;
  }
  dir_->get_rd_lock(pop);

  if ((key_hash >> (8*sizeof(key_hash)-dir_->local_depth)) != dir_->pattern || dir_->sema == -1){
    dir_->release_rd_lock(pop);
    goto RETRY;
  } 
#endif

  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (y+i) % Segment::kNumSlot;
    if (dir_->_[slot].key == key) {
      auto value = dir_->_[slot].value;
#ifdef INPLACE
      dir_->release_rd_lock(pop);
#endif
      return value;
    }
  }

#ifdef INPLACE
  dir_->release_rd_lock(pop);
#endif
  return NONE;
}