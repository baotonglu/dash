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

#define SEGMENT_TYPE 1000
#define DIRECTORY_TYPE 2000
#define ARRAY_TYPE 3000
#define CCEH_TYPE 4000

struct Segment {
  static const size_t kNumSlot = kSegmentSize/sizeof(Pair);

  Segment(void)
  { }

  Segment(size_t depth)
  { }

  ~Segment(void) {}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

  int Insert(PMEMobjpool *pop,Key_t&, Value_t, size_t, size_t);
  void Insert4split(Key_t&, Value_t, size_t);
  bool Put(Key_t&, Value_t, size_t);
  Segment* Split(PMEMobjpool *pop);
  
  void get_lock(PMEMobjpool* pop) {
    //mutex.lock();
    pmemobj_rwlock_wrlock(pop, &rwlock);
    /*
    uint64_t temp = 0;
    while(!CAS(&seg_lock, &temp, 1)){
      temp = 0;
    }*/
  }

  void release_lock(PMEMobjpool* pop) {
    //mutex.unlock();
    pmemobj_rwlock_unlock(pop, &rwlock);
    /*
    uint64_t temp = 1;
    while(!CAS(&seg_lock, &temp, 0)){
      temp = 1;
    }*/
  }

  void get_rd_lock(PMEMobjpool* pop){
    //mutex.lock_shared();
    pmemobj_rwlock_rdlock(pop, &rwlock);
    /*
    uint64_t temp = 0;
    while(!CAS(&seg_lock, &temp, 1)){
      temp = 0;
    }*/
  }

  void release_rd_lock(PMEMobjpool* pop){
    //mutex.unlock_shared();
    pmemobj_rwlock_unlock(pop, &rwlock);
    /*
    uint64_t temp = 1;
    while(!CAS(&seg_lock, &temp, 0)){
      temp = 1;
    }*/
  }

  bool try_get_lock(PMEMobjpool* pop){
    //    return mutex.try_lock();  
    if (pmemobj_rwlock_trywrlock(pop, &rwlock) == 0)
    {
      return true;
    }

    return false;
    /*
    uint64_t temp = 0;
    return CAS(&seg_lock, &temp, 1);
    */
  }

  bool try_get_rd_lock(PMEMobjpool* pop){
   // uint64_t temp = 0;
    //return CAS(&seg_lock, &temp, 1);
    //return mutex.try_lock_shared();
    if (pmemobj_rwlock_tryrdlock(pop, &rwlock) == 0)
    {
      return true;
    }

    return false;
  }

  Pair _[kNumSlot];
  size_t local_depth;
  int64_t sema = 0;
  size_t pattern = 0;
  int count;
  PMEMrwlock rwlock;
  uint64_t seg_lock;
  //std::shared_mutex mutex;
  size_t numElem(void); 
};

struct depth_pattern{
  size_t local_depth;
  size_t pattern;
};

static int create_segment(PMEMobjpool *pop, void *ptr, void *arg){
  struct Segment *se = (struct Segment *)ptr;
  struct depth_pattern *dp = (struct depth_pattern*)arg;

  memset(&se->_[0],255,sizeof(Pair)*Segment::kNumSlot);
  memset(&se->rwlock, 0, sizeof(PMEMrwlock));
  memset(&se->seg_lock, 0, sizeof(se->seg_lock));
  //memset(&se->mutex, 0, sizeof(se->mutex));
  se->count = 0;
  se->local_depth = dp->local_depth;
  se->pattern = dp->pattern;
  se->sema = 0;
  pmemobj_persist(pop, se, sizeof(*se));
  return 0;
}

struct Seg_array{
  uint64_t *_; /*physical addr*/
  PMEMoid _arr; /* persistent pointer addr*/
  size_t global_depth;
};

// the arg is the pointer to the global depth
static int create_seg_array(PMEMobjpool *pop, void *ptr, void *arg){
  struct Seg_array *sa = (struct Seg_array*)ptr;
  sa->_ = nullptr;
  sa->global_depth = *((size_t*)arg);
  sa->_arr = OID_NULL;
  pmemobj_persist(pop, sa, sizeof(*sa));
  return 0;
}

struct Directory {
  static const size_t kDefaultDirectorySize = 1024;
  Seg_array* sa;
  PMEMoid _sa;/* this is the persistent pointer of the seg_array*/
  size_t capacity;
  bool lock;
  int sema = 0;

  Directory(void) {
  }

  Directory(size_t size) {
  }

  ~Directory(void) {
    
  }

  void get_item_num(uint64_t base_addr){
    size_t count = 0;
    Seg_array *seg = sa;
    uint64_t *dir_entry = seg->_;
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

inline void atomic_pointer_update(PMEMoid& a, PMEMoid& b){  
  a.off = b.off; 
}

/* arg is the pointer to the capacity*/
static int create_directory(PMEMobjpool *pop, void *ptr, void *arg){
  struct Directory *dd = (struct Directory *)ptr;
  dd->capacity = *((size_t *)arg);
  dd->lock = false;
  dd->sema = 0;
  dd->sa = nullptr;
  dd->_sa = OID_NULL;
  pmemobj_persist(pop, dd, sizeof(*dd));
  return 0;
}

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
    void Directory_Doubling(PMEMobjpool *pop,int x, Segment*, Segment*);
    void Directory_Update(PMEMobjpool *pop,int x, Segment*, Segment*);
    void Lock_Directory();
    void Unlock_Directory();
    void Get_Number(){
      dir->get_item_num(base_addr);
    }

    void* operator new(size_t size) {
      void *ret;
      posix_memalign(&ret, 64, size);
      return ret;
    }

    PMEMoid _dir;// the pointer to the direcotry
    Directory *dir;
    uint64_t base_addr;
    int seg_num;
    int restart;
};

/* remapp the persistent pointer to the actual address after each open*/
void remapping(PMEMobjpool *pop, CCEH *eh){
  eh->base_addr = (uint64_t)pop;
  if (!OID_IS_NULL(eh->_dir))
  {
     eh->dir = (Directory*)pmemobj_direct(eh->_dir);
    /*reset the directory lock*/
     eh->dir->lock = false;
     if (!OID_IS_NULL(eh->dir->_sa))
     {
       eh->dir->sa = (Seg_array*)pmemobj_direct(eh->dir->_sa);// remapping the address of the segment array
       auto sa = eh->dir->sa;
       if (!OID_IS_NULL(sa->_arr))
       {
         sa->_ = (uint64_t*)pmemobj_direct(sa->_arr);
       }
     }
  }
}

static int create_CCEH(PMEMobjpool *pop, void *ptr, void *arg){
  CCEH *eh = (CCEH *)ptr;
  eh->_dir = OID_NULL;
  eh->dir = nullptr;
  eh->seg_num = 0;
  pmemobj_persist(pop, eh, sizeof(*eh));
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
  eh->base_addr = (uint64_t)pop;
  pmemobj_alloc(pop, &eh->_dir, sizeof(struct Directory), DIRECTORY_TYPE, create_directory, &_capacity);

  /* Initialize the segment_array, the pair_array and the segment pointed by the pair_array*/
    /* Segment array allocation*/
  Directory *dd = (Directory*)pmemobj_direct(eh->_dir);
  eh->dir = dd; 
  pmemobj_alloc(pop, &dd->_sa, sizeof(Seg_array), ARRAY_TYPE, NULL, NULL);
  Seg_array *sa = (Seg_array*)pmemobj_direct(dd->_sa);
  dd->sa = sa;
  sa->global_depth = static_cast<size_t>(log2(dd->capacity));
  //pmemobj_alloc(pop, &sa->_, sizeof(PMEMoid)*dd->capacity, SEGMENT_TYPE+1, NULL, NULL);
  //PMEMoid* dir_entry = (PMEMoid*)pmemobj_direct(sa->_);
  TX_BEGIN(pop){
    pmemobj_tx_add_range_direct(sa,sizeof(*sa));
    sa->_arr = pmemobj_tx_zalloc(sizeof(uint64_t)*dd->capacity, SEGMENT_TYPE+1);
    sa->_ = (uint64_t*)pmemobj_direct(sa->_arr);
  }TX_ONABORT{
    printf("allocation error for direcoty entry during initialization\n");
  }TX_END

  struct depth_pattern dp;
  dp.local_depth = sa->global_depth;

  /* Segment allocation*/
  for (int i = 0; i < dd->capacity; ++i)
  {
    dp.pattern = i;
    segment_allocate(pop, sa->_[i], &dp);
  }

  /* flush the new created semgnets and the direcoty entry array*/
  pmemobj_persist(pop, sa, sizeof(struct Seg_array));
  pmemobj_persist(pop, dd, sizeof(struct Directory));
  pmemobj_persist(pop, eh, sizeof(struct CCEH));
}

#endif  // EXTENDIBLE_PTR_H_

int Segment::Insert(PMEMobjpool *pop, Key_t& key, Value_t value, size_t loc, size_t key_hash) {
#ifdef INPLACE
  if (sema == -1) {
    //printf("retrying because the sema still equals -1\n");
    return 2;
  };
  get_lock(pop);
  if ((key_hash >> (8*sizeof(key_hash)-local_depth)) != pattern || sema == -1){
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

Segment* Segment::Split(PMEMobjpool *pop){
  using namespace std;
#ifndef INPLACE
  int64_t lock = 0;
  if (!CAS(&sema, &lock, -1)) return nullptr;
#else
  if(!try_get_lock(pop))
  {
    return nullptr;
  }
  sema = -1;
#endif

#ifdef INPLACE
  /* How to get the pool handler?*/
  PMEMoid _split;
  size_t new_pattern = (pattern << 1) + 1;
  size_t old_pattern = pattern << 1;
  struct depth_pattern dp;
  dp.local_depth = local_depth + 1;
  dp.pattern = new_pattern; 

  //pmemobj_alloc(pop, &_split, sizeof(Segment*)*2, SEGMENT_TYPE+2, NULL, NULL);
  //Segment** split = (Segment**)pmemobj_direct(_split);
  //memobj_alloc(pop, &split, sizeof(struct Segment), SEGMENT_TYPE, create_segment, &_local_depth);

  //Segment** split = new Segment*[2];
  //split[0] = this;
  pattern = old_pattern;
  PMEMoid split_1;
  pmemobj_alloc(pop, &split_1, sizeof(struct Segment), SEGMENT_TYPE, create_segment, &dp);
  Segment *split = (Segment*)pmemobj_direct(split_1);

  for (unsigned i = 0; i < kNumSlot; ++i) {
    auto key_hash = h(&_[i].key, sizeof(Key_t));
    if (_[i].key != INVALID && key_hash >> (8*8-local_depth-1) == new_pattern)
    {
        split->Insert4split
        (_[i].key, _[i].value, (key_hash & kMask)*kNumPairPerCacheLine);
        split->count++;
        count--;
    }
  }

  assert(count >= 0);
  pmemobj_persist(pop, split, sizeof(struct Segment));

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


CCEH::CCEH(void)
{
}

CCEH::CCEH(size_t initCap)
{
}

CCEH::~CCEH(void)
{ }


void CCEH::Directory_Doubling(PMEMobjpool *pop, int x, Segment *_s0, Segment *_s1){
  Seg_array *sa = dir->sa;
  uint64_t* d = sa->_;
  auto global_depth = sa->global_depth;
  auto s1 = pmemobj_oid(_s1);
  auto s0 = pmemobj_oid(_s0);
//  printf("Directory_Doubling towards %lld\n", global_depth);

  /* new segment array*/
  PMEMoid _new_sa;
  pmemobj_alloc(pop, &_new_sa, sizeof(struct Seg_array), ARRAY_TYPE, NULL, NULL);
  Seg_array *new_sa = (Seg_array*)pmemobj_direct(_new_sa);
  new_sa->global_depth = global_depth + 1;
  pmemobj_alloc(pop, &new_sa->_arr, sizeof(uint64_t)*2*dir->capacity, SEGMENT_TYPE+1, NULL, NULL);
  new_sa->_ = (uint64_t*)pmemobj_direct(new_sa->_arr);
  uint64_t *dd = new_sa->_;

  for (unsigned i = 0; i < dir->capacity; ++i) {
    if (i == x) {
      dd[2*i] = s0.off;
      dd[2*i+1] = s1.off;
    } else {
      dd[2*i] = d[i];
      dd[2*i+1] = d[i];
    }
  }

  pmemobj_persist(pop, dd, sizeof(uint64_t)*2*dir->capacity);
  pmemobj_persist(pop, new_sa, sizeof(struct Seg_array));

  auto old_sa = dir->_sa;
  auto old_dd = sa->_arr;
  /* use txn to make this atomic: capacity of the directory and the segment array*/
  TX_BEGIN(pop){
    //pmemobj_tx_add_range_direct(&my_dir->capacity,sizeof(my_dir->capacity));
    pmemobj_tx_add_range(_dir, 0, sizeof(struct Directory));
    dir->capacity = 2*dir->capacity;
    dir->_sa = _new_sa;
    dir->sa = new_sa;
    pmemobj_free(&old_dd);
    pmemobj_free(&old_sa);
  }TX_END
 // printf("Done!!Directory_Doubling towards %lld\n", global_depth);
  
}

void CCEH::Lock_Directory(){  
  while (!dir->Acquire()) {
      asm("nop");
  }
}

void CCEH::Unlock_Directory(){
  while (!dir->Release()) {
      asm("nop");
  }
}

void CCEH::Directory_Update(PMEMobjpool *pop, int x, Segment *_s0, Segment *_s1){
//  printf("directory update for %d\n", x);
  uint64_t* dir_entry = dir->sa->_;
  auto global_depth = dir->sa->global_depth;
  unsigned depth_diff = global_depth - _s0->local_depth;
  auto s1 = pmemobj_oid(_s1);
  auto s0 = pmemobj_oid(_s0);
  if (depth_diff == 0) {
    if (x%2 == 0) {
      dir_entry[x+1] = s1.off;
#ifdef INPLACE
      pmemobj_persist(pop, &dir_entry[x+1],sizeof(uint64_t));
#else
      dir_entry[x] = s0.off;
      pmemobj_persist(pop, &(dir_entry[x], sizeof(uint64_t)*2));
#endif
    } else {
      dir_entry[x] = s1.off;
#ifdef INPLACE
      pmemobj_persist(pop, &dir_entry[x], sizeof(uint64_t));
#else
      dir_entry[x-1] = s0.off;
      pmemobj_persist(pop, &(dir_entry[x-1]),sizeof(uint64_t)*2);
#endif
    }
  } else {
    int chunk_size = pow(2, global_depth - (_s0->local_depth - 1));
    x = x - (x % chunk_size);
    for (unsigned i = 0; i < chunk_size/2; ++i) {
      dir_entry[x+chunk_size/2+i] = s1.off;
    }
    //clflush((char*)&dir_entry[x+chunk_size/2], sizeof(void*)*chunk_size/2);
    pmemobj_persist(pop, &dir_entry[x+chunk_size/2], sizeof(uint64_t)*chunk_size/2);
#ifndef INPLACE
    for (unsigned i = 0; i < chunk_size/2; ++i) {
      dir_entry[x+i] = s0.off;
    }
    //clflush((char*)&dir._[x], sizeof(void*)*chunk_size/2);
    pmemobj_persist(pop, &dir_entry[x], sizeof(uint64_t)*chunk_size/2);
#endif
  }
//  printf("Done!directory update for %d\n", x);
}

void CCEH::Insert(PMEMobjpool *pop, Key_t& key, Value_t value) {
STARTOVER:
  auto key_hash = h(&key, sizeof(key));
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment* target = (Segment*)(base_addr + dir_entry[x]);
  if (old_sa != dir->sa)
  {
    goto RETRY;
  }

  auto ret = target->Insert(pop, key, value, y, key_hash);
  //std::cout<<"I am inserting "<<key<<", the x = "<<x<<std::endl;

  if (ret == 1) {
    //printf("Segment split for key %lld\n", key);
    Segment* s = target->Split(pop);
    if (s == nullptr) {
      goto RETRY;
    }
    //printf("Done!!Segment split for key %lld\n", key);
    target->local_depth += 1;
    pmemobj_persist(pop, &target->local_depth, sizeof(target->local_depth));

    target->pattern = (key_hash >> (8*sizeof(key_hash)-target->local_depth+1)) << 1;
    s->pattern = ((key_hash >> (8*sizeof(key_hash)-s->local_depth+1)) << 1) + 1;

    // Directory management
    Lock_Directory();
    { // CRITICAL SECTION - directory update
      auto sa = dir->sa;
      dir_entry = sa->_;

      x = (key_hash >> (8*sizeof(key_hash)-sa->global_depth));
      assert(target == (Segment*)(dir_entry[x]+base_addr));
      target = (Segment*)(dir_entry[x]+base_addr);
#ifdef INPLACE
      if (target->local_depth-1<sa->global_depth)
      {
        //std::cout<<"Segment split for key "<<key<<"\n";
#else
      if(target->local_depth< sa->global_depth){
#endif
        Directory_Update(pop,x,target,s);
      } else {  // directory doubling
        Directory_Doubling(pop,x,target,s);
      }
#ifdef INPLACE
      target->sema = 0;
      //s[0]->mutex.unlock();
      target->release_lock(pop);
#endif
    }  // End of critical section
    Unlock_Directory();
    goto RETRY;
  } else if (ret == 2) {
    goto STARTOVER;
  } 
}

// This function does not allow resizing
bool CCEH::InsertOnly(PMEMobjpool* pop, Key_t& key, Value_t value) {
  auto key_hash = h(&key, sizeof(key));
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment* target = (Segment*)(base_addr + dir_entry[x]);
  if (old_sa != dir->sa)
  {
    goto RETRY;
  }

  auto ret = target->Insert(pop, key, value, y, key_hash);
  if (ret == 0)
  {
    return true;
  }
  return false;
}

// TODO
bool CCEH::Delete(Key_t& key) {
  return false;
}

Value_t CCEH::Get(PMEMobjpool*pop, Key_t& key) {
  auto key_hash = h(&key, sizeof(key));
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment* dir_ = (Segment*)(base_addr + dir_entry[x]);
  if (old_sa != dir->sa)
  {
    goto RETRY;
  }

#ifdef INPLACE
  auto sema = dir_->sema;

  if (sema == -1)
  {
    goto RETRY;
  }
  //dir_->mutex.lock_shared();
  dir_->get_rd_lock(pop);

  if ((key_hash >> (8*sizeof(key_hash)-dir_->local_depth)) != dir_->pattern || dir_->sema == -1){
    //dir_->mutex.unlock_shared();
    dir_->release_rd_lock(pop);
    goto RETRY;
  } 

#endif

  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (y+i) % Segment::kNumSlot;
    if (dir_->_[slot].key == key) {
      auto value = dir_->_[slot].value;
#ifdef INPLACE
      //dir_->mutex.unlock_shared();
      dir_->release_rd_lock(pop);
#endif
      return value;
    }
  }

#ifdef INPLACE
  //dir_->mutex.unlock_shared();
  dir_->release_rd_lock(pop);
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
