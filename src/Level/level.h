#ifndef LEVEL_HASHING_H_
#define LEVEL_HASHING_H_

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdint.h>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <libpmemobj.h>
#include "util/pair.h"
#include "util/hash.h"
#include "util/persist.h"
#include <sys/stat.h>
#include <inttypes.h> 
#define ASSOC_NUM 7
#define NODE_TYPE 1000
#define LEVEL_TYPE 2000
#define LOCK_TYPE 3000

#define BATCH 1
#define UNIQUE_CHECK 1

struct Entry {
  Key_t key;
  Value_t value;
  Entry() {
    key = INVALID;
    value = NONE;
  }
};

struct Node {
  uint8_t token[ASSOC_NUM];
  Entry slot[ASSOC_NUM];
  char dummy[1];
  void* operator new[] (size_t size) {
    void* ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }

  void* operator new (size_t size) {
    void* ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }
};

class LevelHashing{
  public:
    PMEMoid _buckets[2];
    PMEMoid _interim_level_buckets;
    Node* buckets[2];
    Node* interim_level_buckets; 
    //uint64_t level_item_num[2];

    uint64_t levels;
    uint64_t addr_capacity;
    uint64_t total_capacity;
//    uint32_t occupied;
    uint64_t f_seed;
    uint64_t s_seed;
    uint32_t resize_num;
    //int32_t resizing_lock = 0;
    std::atomic<int64_t> resizing_lock;
    PMEMoid _mutex;
    PMEMoid _old_mutex;
    int prev_nlocks;
    int nlocks;
    int locksize;
    bool resizing;

    //void generate_seeds(void);
    void resize(PMEMobjpool *pop);
    int b2t_movement(PMEMobjpool *pop, uint64_t );
    uint8_t try_movement(PMEMobjpool *pop, uint64_t , uint64_t , Key_t& , Value_t);
    uint64_t F_HASH(uint64_t key);
    uint64_t S_HASH(uint64_t key);
  
    LevelHashing(void);
    LevelHashing(size_t);
    ~LevelHashing(void);
     void display_size(){
      std::cout<<"The Node size is "<<sizeof(Node)<<std::endl;
      std::cout<<"The Entry size is "<<sizeof(Entry)<<std::endl;
      
      if (buckets[0])
      {
        printf("0x%" PRIXPTR "\n",(uintptr_t)buckets[0]);
      }

      if (buckets[1])
      {
        printf("0x%" PRIXPTR "\n",(uintptr_t)buckets[1]);
      }
    }

    //bool InsertOnly(Key_t&, Value_t);
    void Insert(PMEMobjpool*, Key_t&, Value_t);
    bool Delete(PMEMobjpool*, Key_t&);
    Value_t Get(PMEMobjpool*, Key_t&);
    double Utilization(void);
    size_t Capacity(void) {
      return (addr_capacity + addr_capacity/2)*ASSOC_NUM;
    }
};

#define F_IDX(hash, capacity) (hash % (capacity/2))
#define S_IDX(hash, capacity) ((hash % (capacity/2)) + (capacity/2))

uint64_t LevelHashing::F_HASH(uint64_t key){
  return h(&key, sizeof(Key_t), f_seed);
}

uint64_t LevelHashing::S_HASH(uint64_t key){
  return h(&key, sizeof(Key_t), s_seed);
}

void generate_seeds(LevelHashing *level) {
  srand(time(NULL));
  do
  {
    level->f_seed = rand();
    level->s_seed = rand();
    level->f_seed = level->f_seed << (rand() % 63);
    level->s_seed = level->s_seed << (rand() % 63);
  } while (level->f_seed == level->s_seed);
}

LevelHashing::LevelHashing(void){
}

LevelHashing::~LevelHashing(void){
 
}
/*
LevelHashing::LevelHashing(size_t _levels)
  : levels{_levels},
  resize_num{0},
  resizing_lock{0}
{
  resizing = false;
  addr_capacity = pow(2, levels);
  total_capacity = pow(2, levels) + pow(2, levels-1);
  locksize = 256;
  nlocks = (3*addr_capacity/2)/locksize+1;
  mutex = new std::shared_mutex[nlocks];

  generate_seeds();
  buckets[0] = new Node[addr_capacity];
  buckets[1] = new Node[addr_capacity/2];
  level_item_num[0] = 0;
  level_item_num[1] = 0;
  interim_level_buckets = NULL;
}*/

void* cache_align(void* ptr){
  //ptr +=  48;
  uint64_t pp = (uint64_t)ptr;
  pp += 48;
  return (void *)pp;
}

void initialize_level(PMEMobjpool *pop, LevelHashing *level, void *arg){
  TX_BEGIN(pop){
    pmemobj_tx_add_range_direct(level, sizeof(*level));
    /* modify*/
    level->levels = *((int*)arg);
    level->resize_num = 0;
    level->resizing_lock = 0;
    level->resizing = false;
    level->addr_capacity = pow(2, level->levels);
    level->total_capacity = pow(2, level->levels) + pow(2, level->levels-1);
    level->locksize = 128;
    level->nlocks = (3*level->addr_capacity/2)/level->locksize+1;

    generate_seeds(level);
    //level->level_item_num[0] = 0;
    //level->level_item_num[1] = 0;
    level->_interim_level_buckets = OID_NULL;
    level->interim_level_buckets = NULL;
    /*allocate*/
    level->_mutex = pmemobj_tx_zalloc(sizeof(PMEMrwlock)*level->nlocks, LOCK_TYPE);
    level->_buckets[0] = pmemobj_tx_zalloc(sizeof(Node)*level->addr_capacity + 1, NODE_TYPE);
    level->_buckets[1] = pmemobj_tx_zalloc(sizeof(Node)*level->addr_capacity/2 + 1, NODE_TYPE);
    level->_old_mutex = OID_NULL;
    level->prev_nlocks = 0;

    /* Intialize pointer*/
    level->buckets[0] = (Node*)cache_align(pmemobj_direct(level->_buckets[0]));
    level->buckets[1] = (Node*)cache_align(pmemobj_direct(level->_buckets[1]));
  }TX_END
}

void remapping(LevelHashing *level){
    level->buckets[0] = (Node*)cache_align(pmemobj_direct(level->_buckets[0]));
    level->buckets[1] = (Node*)cache_align(pmemobj_direct(level->_buckets[1]));
}

void LevelHashing::Insert(PMEMobjpool *pop, Key_t& key, Value_t value) {
RETRY:
//std::cout<<"Inserting key "<<key<<std::endl;
  while (resizing_lock.load() == 1) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);

  PMEMrwlock* mutex = (PMEMrwlock*)pmemobj_direct(_mutex);

#ifdef UNIQUE_CHECK
  uint32_t lock_idx = f_idx / locksize;
  while(pmemobj_rwlock_trywrlock(pop, &mutex[lock_idx]) != 0){
    if (resizing == true)
    {
      goto RETRY;
    }
  }

  if(resizing == true){
    pmemobj_rwlock_unlock(pop, &mutex[lock_idx]);
    goto RETRY;
  }

  bool inserted = false;
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < ASSOC_NUM; ++j){
      if ((buckets[i][f_idx].token[j] == 1) && (buckets[i][f_idx].slot[j].key == key)){
        inserted = true;
        goto UNIQUE;
      }
    }

    for(int j = 0; j < ASSOC_NUM; ++j){
      if ((buckets[i][s_idx].token[j] == 1) && (buckets[i][s_idx].slot[j].key == key)){
        inserted = true;
        goto UNIQUE;
      }
    }

    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }

UNIQUE:
  if(inserted){
    pmemobj_rwlock_unlock(pop, &mutex[lock_idx]);
    return;
  }

  pmemobj_rwlock_unlock(pop, &mutex[lock_idx]);
  f_idx = F_IDX(f_hash, addr_capacity);
  s_idx = S_IDX(s_hash, addr_capacity);
#endif

  int i, j;
  for(i = 0; i < 2; i ++){
    while(pmemobj_rwlock_trywrlock(pop, &mutex[f_idx/locksize]) != 0){
      if (resizing == true)
      {
        goto RETRY;
      }
    }

    if (resizing == true)
    {
      pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
      goto RETRY;
    }

    for(j = 0; j < ASSOC_NUM; j ++){
        if(buckets[i][f_idx].token[j] == 0){
          buckets[i][f_idx].slot[j].value = value;
          buckets[i][f_idx].slot[j].key = key;
          pmemobj_persist(pop, &buckets[i][f_idx].slot[j], sizeof(Entry));
          mfence();
          buckets[i][f_idx].token[j] = 1;
          //clflush((char*)&buckets[i][f_idx], sizeof(Node));
          pmemobj_persist(pop, &buckets[i][f_idx].token[j], sizeof(uint8_t));
          //level_item_num[i]++;
          //mutex[f_idx/locksize].unlock();
          pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
          return;
        }
    }
    pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);

    while(pmemobj_rwlock_trywrlock(pop, &mutex[s_idx/locksize]) != 0){
      if (resizing == true)
      {
        goto RETRY;
      }
    }

    if (resizing == true)
    {
      pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
      goto RETRY;
    }

    for(j = 0; j < ASSOC_NUM; j ++){
        if(buckets[i][s_idx].token[j] == 0){
          buckets[i][s_idx].slot[j].value = value;
          buckets[i][s_idx].slot[j].key = key;
          pmemobj_persist(pop, &buckets[i][s_idx].slot[j], sizeof(Entry));
          mfence();
          buckets[i][s_idx].token[j] = 1;
          pmemobj_persist(pop, &buckets[i][s_idx].token[j], sizeof(uint8_t));
          //level_item_num[i]++;
          pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
          return;
        }
    }
    pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);

    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }

  f_idx = F_IDX(f_hash, addr_capacity);
  s_idx = S_IDX(s_hash, addr_capacity);
  int empty_loc;
  int64_t lock = 0;
  //if (CAS(&resizing_lock, &lock, 1)) {
  if (resizing_lock.compare_exchange_strong(lock, 1))
  {
  
    for(i=0; i<2; i++){
      if(!try_movement(pop,f_idx, i, key, value)){
        //resizing_lock = 0;
        resizing_lock.store(0);
        return;
      }
      if(!try_movement(pop,s_idx, i, key, value)){
        //resizing_lock = 0;
        resizing_lock.store(0);
        return;
      }
      f_idx = F_IDX(f_hash, addr_capacity / 2);
      s_idx = S_IDX(s_hash, addr_capacity / 2);
    }

    if(resize_num>0){
      {
        //std::unique_lock<std::mutex> lock(mutex[f_idx/locksize]);
        //mutex[f_idx/locksize].lock();
        pmemobj_rwlock_wrlock(pop, &mutex[f_idx/locksize]);
#ifdef TIME
  cuck_timer.Start();
#endif
        empty_loc = b2t_movement(pop,f_idx);
#ifdef TIME
  cuck_timer.Stop();
  displacement += cuck_timer.GetSeconds();
#endif
        if(empty_loc != -1){
          buckets[1][f_idx].slot[empty_loc].value = value;
          buckets[1][f_idx].slot[empty_loc].key = key;
          pmemobj_persist(pop, &buckets[1][f_idx].slot[empty_loc], sizeof(Entry));
          mfence();
          buckets[1][f_idx].token[empty_loc] = 1;
          pmemobj_persist(pop, &buckets[1][f_idx].token[empty_loc], sizeof(uint8_t));
          //level_item_num[1]++;
          resizing_lock.store(0);
          //mutex[f_idx/locksize].unlock();
          pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
          return;
        }
       pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
      }
      {
        //std::unique_lock<std::mutex> lock(mutex[s_idx/locksize]);
        //mutex[s_idx/locksize].lock();
        pmemobj_rwlock_wrlock(pop, &mutex[s_idx/locksize]);
#ifdef TIME
  cuck_timer.Start();
#endif
        empty_loc = b2t_movement(pop, s_idx);
#ifdef TIME
  cuck_timer.Stop();
  displacement += cuck_timer.GetSeconds();
#endif
        if(empty_loc != -1){
          buckets[1][s_idx].slot[empty_loc].value = value;
          buckets[1][s_idx].slot[empty_loc].key = key;
          pmemobj_persist(pop, &buckets[1][s_idx].slot[empty_loc], sizeof(Entry));
          mfence();
          buckets[1][s_idx].token[empty_loc] = 1;
          //clflush((char*)&buckets[1][s_idx], sizeof(Node));
          pmemobj_persist(pop, &buckets[1][s_idx].token[empty_loc], sizeof(uint8_t));
          //level_item_num[1]++;
          //resizing_lock = 0;
          resizing_lock.store(0);
          //mutex[s_idx/locksize].unlock();
          pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
          return;
        }
       pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
      }
    }
    //timer.Start();
    resize(pop);
    //timer.Stop();
    //breakdown += timer.GetSeconds();
    //resizing_lock = 0;
    resizing_lock.store(0);
  }
  goto RETRY;
}
/*
bool LevelHashing::InsertOnly(Key_t& key, Value_t value) {
  return false;
}*/

void LevelHashing::resize(PMEMobjpool *pop) {
  std::cout<<"Resizing towards levels "<<levels+1<<std::endl;

  PMEMrwlock* mutex = (PMEMrwlock*)pmemobj_direct(_mutex);
  resizing = true;
  for (int i = 0; i < nlocks; ++i)
  {
    //mutex[i].lock();
    pmemobj_rwlock_wrlock(pop, &mutex[i]);
  }
  //std::shared_mutex* old_mutex = mutex;
  
  //nlocks = nlocks + 2*addr_capacity/locksize+1;
  size_t new_addr_capacity = pow(2, levels + 1);
  TX_BEGIN(pop){
      pmemobj_tx_add_range_direct(&nlocks, sizeof(nlocks));
      pmemobj_tx_add_range_direct(&prev_nlocks, sizeof(prev_nlocks));
      pmemobj_tx_add_range_direct(&_old_mutex, sizeof(_old_mutex));
      pmemobj_tx_add_range_direct(&_mutex, sizeof(_mutex));
      pmemobj_tx_add_range_direct(&interim_level_buckets, sizeof(interim_level_buckets));
      pmemobj_tx_add_range_direct(&_interim_level_buckets, sizeof(_interim_level_buckets));

      _old_mutex = _mutex;
      prev_nlocks = nlocks;
      nlocks = (3*2*addr_capacity/2)/locksize+1;
      _mutex = pmemobj_tx_zalloc(nlocks*sizeof(PMEMrwlock), LOCK_TYPE);
      _interim_level_buckets = pmemobj_tx_zalloc(new_addr_capacity*sizeof(Node)+1, NODE_TYPE);
      interim_level_buckets = (Node*)cache_align(pmemobj_direct(_interim_level_buckets));
  }TX_ONABORT{
    printf("resizing txn 1 fails\n");
  }TX_END

  PMEMrwlock* old_mutex = mutex;
  mutex = (PMEMrwlock*)pmemobj_direct(_mutex); 

  //clflush((char*)&interim_level_buckets, sizeof(Node));

  //uint64_t new_level_item_num = 0;
  uint64_t old_idx;
  for (old_idx = 0; old_idx < pow(2, levels - 1); old_idx ++) {
    uint64_t i, j;
    for(i = 0; i < ASSOC_NUM; i ++){
      if (buckets[1][old_idx].token[i] == 1)
      {
        Key_t key = buckets[1][old_idx].slot[i].key;
        Value_t value = buckets[1][old_idx].slot[i].value;

        uint32_t f_idx = F_IDX(F_HASH(key), new_addr_capacity);
        uint32_t s_idx = S_IDX(S_HASH(key), new_addr_capacity);

        uint8_t insertSuccess = 0;
        for(j = 0; j < ASSOC_NUM; j ++){
          if (interim_level_buckets[f_idx].token[j] == 0)
          {
            interim_level_buckets[f_idx].slot[j].value = value;
            interim_level_buckets[f_idx].slot[j].key = key;
#ifndef BATCH
            pmemobj_persist(pop, &interim_level_buckets[f_idx].slot[j], sizeof(Entry));
            mfence();
#endif
      interim_level_buckets[f_idx].token[j] = 1;
#ifndef BATCH
      //clflush((char*)&interim_level_buckets[f_idx], sizeof(Node));
      pmemobj_persist(pop, &interim_level_buckets[f_idx].token[j], sizeof(uint8_t));
#endif
            insertSuccess = 1;
      //new_level_item_num++;
            break;
          }
          else if (interim_level_buckets[s_idx].token[j] == 0)
          {
            interim_level_buckets[s_idx].slot[j].value = value;
            interim_level_buckets[s_idx].slot[j].key = key;
#ifndef BATCH
            pmemobj_persist(pop, &interim_level_buckets[s_idx].slot[j], sizeof(Entry));
            mfence();
#endif
      interim_level_buckets[s_idx].token[j] = 1;
#ifndef BATCH
      //clflush((char*)&interim_level_buckets[s_idx], sizeof(Node));
      pmemobj_persist(pop, &interim_level_buckets[s_idx].token[j], sizeof(uint8_t));
#endif
            insertSuccess = 1;
      //new_level_item_num++;
            break;
          }
        }

#ifndef BATCH
  buckets[1][old_idx].token[i] = 0;
  //clflush((char*)&buckets[1][old_idx].token[i], sizeof(uint8_t));
  pmemobj_persist(pop, &buckets[1][old_idx].token[i], sizeof(uint8_t));
#endif
      }
    }
  }

#ifdef BATCH
  pmemobj_persist(pop, &buckets[1][0],sizeof(Node)*pow(2,levels-1));
  pmemobj_persist(pop, &interim_level_buckets[0], sizeof(Node)*new_addr_capacity);
  //clflush((char*)&buckets[1][0],sizeof(Node)*pow(2,levels-1));
  //clflush((char*)&interim_level_buckets[0], sizeof(Node)*new_addr_capacity);
#endif

  TX_BEGIN(pop){
      pmemobj_tx_add_range_direct(this, sizeof(class LevelHashing));

      levels++;
      resize_num++;
      pmemobj_tx_free(_buckets[1]);
      _buckets[1] = _buckets[0];
      _buckets[0] = _interim_level_buckets;
      buckets[1] = (Node*)cache_align(pmemobj_direct(_buckets[1]));
      buckets[0] = (Node*)cache_align(pmemobj_direct(_buckets[0]));

      _interim_level_buckets = OID_NULL;
      interim_level_buckets = NULL;
     // level_item_num[1] = level_item_num[0];
      //level_item_num[0] = new_level_item_num;

      addr_capacity = new_addr_capacity;
      total_capacity = pow(2, levels) + pow(2, levels - 1);
      pmemobj_tx_free(_old_mutex);
      resizing = false;
  } TX_ONABORT{
    printf("resizing txn 2 fails\n");
  }TX_END
/*
  levels++;
  resize_num++;

  delete [] buckets[1];
  buckets[1] = buckets[0];
  buckets[0] = interim_level_buckets;
  interim_level_buckets = NULL;

  level_item_num[1] = level_item_num[0];
  level_item_num[0] = new_level_item_num;

  addr_capacity = new_addr_capacity;
  total_capacity = pow(2, levels) + pow(2, levels - 1);

  for(int i=0;i<prev_nlocks;i++){
    old_mutex[i].unlock();
  }
  //delete [] old_mutex;
  delete old_mutex;
  resizing = false;*/
  std::cout<<"Done! :Resizing towards levels "<<levels<<std::endl;
}

uint8_t LevelHashing::try_movement(PMEMobjpool *pop, uint64_t idx, uint64_t level_num, Key_t& key, Value_t value) {
#ifdef TIME
cuck_timer.Start();
#endif
  uint64_t i, j, jdx;
  PMEMrwlock* mutex = (PMEMrwlock*)pmemobj_direct(_mutex);
  {
    //std::unique_lock<std::mutex> *lock[2];
    //lock[0] = new std::unique_lock<std::mutex>(mutex[idx/locksize]);
    //mutex[idx/locksize].lock();
    pmemobj_rwlock_wrlock(pop, &mutex[idx/locksize]);
    for(i=0; i<ASSOC_NUM; i++){
      Key_t m_key = buckets[level_num][idx].slot[i].key;
      Value_t m_value = buckets[level_num][idx].slot[i].value;
      uint64_t f_hash = F_HASH(m_key);
      uint64_t s_hash = S_HASH(m_key);
      uint64_t f_idx = F_IDX(f_hash, addr_capacity/(1+level_num));
      uint64_t s_idx = S_IDX(s_hash, addr_capacity/(1+level_num));

      if(f_idx == idx) jdx = s_idx;
      else jdx = f_idx;

      if((jdx/locksize)!=(idx/locksize)){
        //lock[1] = new std::unique_lock<std::mutex>(mutex[jdx/locksize]);
        //mutex[jdx/locksize].lock();
         pmemobj_rwlock_wrlock(pop, &mutex[jdx/locksize]);
      }

      for(j=0; j<ASSOC_NUM; j++){
        if(buckets[level_num][jdx].token[j] == 0){
          buckets[level_num][jdx].slot[j].value = m_value;
          buckets[level_num][jdx].slot[j].key = m_key;
          pmemobj_persist(pop, &buckets[level_num][jdx].slot[j], sizeof(Entry));
          mfence();
          buckets[level_num][jdx].token[j] = 1;
          pmemobj_persist(pop, &buckets[level_num][jdx].token[j], sizeof(uint8_t));
          buckets[level_num][idx].token[i] = 0;
          //clflush((char*)&buckets[level_num][idx].token[i], sizeof(uint8_t));
          pmemobj_persist(pop, &buckets[level_num][idx].token[i], sizeof(uint8_t));

          buckets[level_num][idx].slot[i].value = value;
          buckets[level_num][idx].slot[i].key = key;
          pmemobj_persist(pop, &buckets[level_num][idx].slot[i], sizeof(Entry));
          mfence();
          buckets[level_num][idx].token[i] = 1;
          pmemobj_persist(pop, &buckets[level_num][idx].token[i], sizeof(uint8_t));
          //level_item_num[level_num]++;

          if((jdx/locksize) != (idx/locksize)) {
            //delete lock[1];
            //mutex[jdx/locksize].unlock();
            pmemobj_rwlock_unlock(pop, &mutex[jdx/locksize]);
          }
          //delete lock[0];
          //mutex[idx/locksize].unlock();
          pmemobj_rwlock_unlock(pop, &mutex[idx/locksize]);
#ifdef TIME
    cuck_timer.Stop();
    displacement += cuck_timer.GetSeconds();
#endif
          return 0;
        }
      }
      if((jdx/locksize) != (idx/locksize)) pmemobj_rwlock_unlock(pop, &mutex[jdx/locksize]);
    }
    //delete lock[0];
    pmemobj_rwlock_unlock(pop, &mutex[idx/locksize]);
  }
#ifdef TIME
  cuck_timer.Stop();
  displacement += cuck_timer.GetSeconds();
#endif
  return 1;
}

int LevelHashing::b2t_movement(PMEMobjpool *pop, uint64_t idx){
  Key_t key;
  Value_t value;
  uint64_t s_hash, f_hash;
  uint64_t s_idx, f_idx;
  uint64_t i, j;
  PMEMrwlock* mutex = (PMEMrwlock*)pmemobj_direct(_mutex);

  //std::unique_lock<std::mutex> *lock;
  for(i=0; i<ASSOC_NUM; i++){
    key = buckets[1][idx].slot[i].key;
    value = buckets[1][idx].slot[i].value;
    f_hash = F_HASH(key);
    s_hash = S_HASH(key);
    f_idx = F_IDX(f_hash, addr_capacity);
    s_idx = S_IDX(s_hash, addr_capacity);

    for(j=0; j<ASSOC_NUM; j++){
      if((idx/locksize) != (f_idx/locksize))
        //lock = new std::unique_lock<std::mutex>(mutex[f_idx/locksize]);
        //mutex[f_idx/locksize].lock();
        pmemobj_rwlock_wrlock(pop, &mutex[f_idx/locksize]);

      if(buckets[0][f_idx].token[j] == 0){
        buckets[0][f_idx].slot[j].value = value;
        buckets[0][f_idx].slot[j].key = key;
        pmemobj_persist(pop, &buckets[0][f_idx].slot[j], sizeof(Entry));
        mfence();
        buckets[0][f_idx].token[j] = 1;
        //clflush((char*)&buckets[0][f_idx], sizeof(Node));
        pmemobj_persist(pop, &buckets[0][f_idx].token[j], sizeof(uint8_t));
        buckets[1][idx].token[i] = 0;
        //clflush((char*)&buckets[1][idx].token[i], sizeof(uint8_t));
        pmemobj_persist(pop, &buckets[1][idx].token[i], sizeof(uint8_t));
        //level_item_num[0]++;
        //level_item_num[1]--;

        if((idx/locksize) != (f_idx/locksize)) pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
        return i;
      }
      //if((idx/locksize)!=(f_idx/locksize)) delete lock;
      if((idx/locksize) != (f_idx/locksize)) pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
      if((idx/locksize)!=(s_idx/locksize))
        //lock = new std::unique_lock<std::mutex>(mutex[s_idx/locksize]);
        pmemobj_rwlock_wrlock(pop, &mutex[s_idx/locksize]);

      if(buckets[0][s_idx].token[j] == 0){
        buckets[0][s_idx].slot[j].value = value;
        buckets[0][s_idx].slot[j].key = key;
        pmemobj_persist(pop, &buckets[0][s_idx].slot[j], sizeof(Entry));
        mfence();
  buckets[0][s_idx].token[j] = 1;
  //clflush((char*)&buckets[0][s_idx], sizeof(Node));
  pmemobj_persist(pop, &buckets[0][s_idx].token[j], sizeof(uint8_t));
  buckets[1][idx].token[i] = 0;
  //clflush((char*)&buckets[0][s_idx].token[j], sizeof(uint8_t));
  pmemobj_persist(pop, &buckets[0][s_idx].token[j], sizeof(uint8_t));

  //level_item_num[0]++;
  //level_item_num[1]--;

        if((idx/locksize) != (s_idx/locksize)) pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
        return i;
      }
      if((idx/locksize)!=(s_idx/locksize)) pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
    }
  }
  return -1;
}


Value_t LevelHashing::Get(PMEMobjpool *pop ,Key_t& key) {
RETRY:
  while (resizing == true) {
      asm("nop");
  }
  PMEMrwlock* mutex = (PMEMrwlock*)pmemobj_direct(_mutex);
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);
  int i = 0, j;

  for(i = 0; i < 2; i ++){
    {
      //mutex[f_idx/locksize].lock_shared();
      while(pmemobj_rwlock_tryrdlock(pop, &mutex[f_idx/locksize]) != 0){
        if (resizing == true)
        {
          goto RETRY;
        }
      }

      if (resizing == true)
      {
        //mutex[f_idx/locksize].unlock_shared();
        pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
        goto RETRY;
      }

      for(j = 0; j < ASSOC_NUM; j ++){
        if (buckets[i][f_idx].token[j] == 1 && buckets[i][f_idx].slot[j].key == key)
        {
          //mutex[f_idx/locksize].unlock_shared();
          pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
          return buckets[i][f_idx].slot[j].value;
        }
      }
      //mutex[f_idx/locksize].unlock_shared();
      pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
    }
    {
      //]mutex[s_idx/locksize].lock_shared();
      while(pmemobj_rwlock_tryrdlock(pop, &mutex[s_idx/locksize]) != 0){
        if (resizing == true)
        {
          goto RETRY;
        }
      }

      if (resizing == true)
      {
        //mutex[s_idx/locksize].unlock_shared();
        pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
        goto RETRY;
      }

      for(j = 0; j < ASSOC_NUM; j ++){
        if (buckets[i][s_idx].token[j] == 1 && buckets[i][s_idx].slot[j].key == key)
        {
          pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
          return buckets[i][s_idx].slot[j].value;
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
    }
    f_idx = F_IDX(f_hash, addr_capacity/2);
    s_idx = S_IDX(s_hash, addr_capacity/2);
  }
  return NONE;
}

bool LevelHashing::Delete(PMEMobjpool *pop ,Key_t& key) {
  RETRY:
  while (resizing == true) {
      asm("nop");
  }
  PMEMrwlock* mutex = (PMEMrwlock*)pmemobj_direct(_mutex);
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);
  int i = 0, j;

  for(i = 0; i < 2; i++){
    {
      while(pmemobj_rwlock_tryrdlock(pop, &mutex[f_idx/locksize]) != 0){
        if (resizing == true)
        {
          goto RETRY;
        }
      }

      if (resizing == true)
      {
        pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
        goto RETRY;
      }

      for(j = 0; j < ASSOC_NUM; j ++){
        if (buckets[i][f_idx].token[j] == 1 && buckets[i][f_idx].slot[j].key == key)
        {
          buckets[i][f_idx].token[j] = 0;
          pmemobj_persist(pop, &buckets[i][f_idx].token[j], sizeof(uint8_t));
          pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
          return true;
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[f_idx/locksize]);
    }
    {
      while(pmemobj_rwlock_tryrdlock(pop, &mutex[s_idx/locksize]) != 0){
        if (resizing == true)
        {
          goto RETRY;
        }
      }

      if (resizing == true)
      {
        pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
        goto RETRY;
      }

      for(j = 0; j < ASSOC_NUM; j ++){
        if (buckets[i][s_idx].token[j] == 1 && buckets[i][s_idx].slot[j].key == key)
        {
	  buckets[i][s_idx].token[j] = 0;
          pmemobj_persist(pop, &buckets[i][s_idx].token[j], sizeof(uint8_t));
          pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
          return true;
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[s_idx/locksize]);
    }
    f_idx = F_IDX(f_hash, addr_capacity/2);
    s_idx = S_IDX(s_hash, addr_capacity/2);
  }
  return false;
}

#endif  // LEVEL_HASHING_H_
