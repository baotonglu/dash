#ifndef LEVEL_HASHING_H_
#define LEVEL_HASHING_H_

/*
We do several optimization and correctness patches for level hashing, including:
(1) remove fence between storing value and storing key during insert() because
these two stores are in the same cacheline and will mot be reordered. (2) use
lock-striping rather than slot-based locking to make the lock area fitting into
the CPU cache, improving the scalability. (3) use persistent lock in PMDK
library to aovid deadlock caused by sudden system failure. (4) add uniqnuess
check during the insert opeartion to avoid inserting duplicate keys. (5) add
support for variable-length key by storing the pointer to the key object.
*/

#include <inttypes.h>
#include <libpmem.h>
#include <libpmemobj.h>
#include <stdint.h>
#include <sys/stat.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <mutex>
#include <shared_mutex>

#include "../../util/hash.h"
#include "../../util/pair.h"
#include "../../util/utils.h"
#include "../Hash.h"
#define ASSOC_NUM 7
#define NODE_TYPE 1000
#define LEVEL_TYPE 2000
#define LOCK_TYPE 3000
//#define COUNTING 1
//#define BATCH 1
#define UNIQUE_CHECK 1
namespace level {
inline bool var_compare(char *str1, char *str2, int len1, int len2) {
  if (len1 != len2) return false;
  return !memcmp(str1, str2, len1);
}

template <class T>
struct Entry {
  T key;
  Value_t value;
  Entry() {
    key = INVALID;
    value = NONE;
  }
};

/*
 * 1 cacheline setting: ASSOC_NUM = 3; dummy = 13
 * 2 cacheline setting: ASSOC_NUM = 7; dummy = 9
 * 4 cacheline setting: ASSOC_NUM = 15; dummy = 1
 */

template <class T>
struct Node {
  uint8_t token[ASSOC_NUM];
  char dummy[9];
  Entry<T> slot[ASSOC_NUM];
  void *operator new[](size_t size) {
    void *ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }

  void *operator new(size_t size) {
    void *ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }
};

template <class T>
class LevelHashing : public Hash<T> {
 public:
  PMEMobjpool *pop;
  PMEMoid _buckets[2];
  PMEMoid _interim_level_buckets;
  Node<T> *buckets[2];
  Node<T> *interim_level_buckets;
  uint64_t level_item_num[2];

  uint64_t levels;
  uint64_t addr_capacity;
  uint64_t total_capacity;
  uint64_t f_seed;
  uint64_t s_seed;
  uint32_t resize_num;
  /* lock information */
  std::atomic<int64_t> resizing_lock;
  PMEMoid _mutex;
  PMEMoid _old_mutex; /*old mutex lock*/
  PMEMrwlock *mutex;
  int nlocks;
  int locksize;
  bool resizing;

  void generate_seeds(void);
  void resize(PMEMobjpool *pop);
  int b2t_movement(PMEMobjpool *pop, uint64_t);
  uint8_t try_movement(PMEMobjpool *pop, uint64_t, uint64_t, T, Value_t);
  uint64_t F_HASH(T key);
  uint64_t S_HASH(T key);

  LevelHashing(void);
  LevelHashing(size_t);
  ~LevelHashing(void);
  void display_size() {
    std::cout << "The Node size is " << sizeof(Node<T>) << std::endl;
    std::cout << "The Entry size is " << sizeof(Entry<T>) << std::endl;

    if (buckets[0]) {
      printf("0x%" PRIXPTR "\n", (uintptr_t)buckets[0]);
    }

    if (buckets[1]) {
      printf("0x%" PRIXPTR "\n", (uintptr_t)buckets[1]);
    }
  }

  int Insert(T, Value_t);
  int Insert(T key, Value_t value, bool is_in_epoch) {
    return Insert(key, value);
  }
  bool Delete(T);
  bool Delete(T key, bool is_in_epoch) { return Delete(key); }
  Value_t Get(T);
  Value_t Get(T key, bool flag) { return Get(key); }
  void Recovery() {
    if (resizing) {
      if (!OID_IS_NULL(_old_mutex) && !OID_EQUALS(_old_mutex, _mutex)) {
        pmemobj_free(&_old_mutex);
      }
      if (!OID_IS_NULL(_interim_level_buckets)) {
        pmemobj_free(&_interim_level_buckets);
      }
      resizing = false;
      resizing_lock = 0;
    }
  }
  void getNumber() {
    std::cout << "Entry Size: " << sizeof(struct Entry<T>) << std::endl;
    std::cout << "Node Size: " << sizeof(struct Node<T>) << std::endl;
    std::cout << "First items " << level_item_num[0] << std::endl;
    std::cout << "Second items " << level_item_num[1] << std::endl;
    std::cout << "First load factor = "
              << (double)level_item_num[0] / (addr_capacity * ASSOC_NUM)
              << std::endl;
    std::cout << "Second load factor = "
              << (double)level_item_num[1] / (addr_capacity * ASSOC_NUM / 2)
              << std::endl;
  }
  size_t Capacity(void) {
    return (addr_capacity + addr_capacity / 2) * ASSOC_NUM;
  }
};

#define F_IDX(hash, capacity) (hash % (capacity / 2))
#define S_IDX(hash, capacity) ((hash % (capacity / 2)) + (capacity / 2))

template <class T>
uint64_t LevelHashing<T>::F_HASH(T key) {
  if constexpr (std::is_pointer_v<T>) {
    return h(key->key, key->length, f_seed);
  } else {
    return h(&key, sizeof(Key_t), f_seed);
  }
}

template <class T>
uint64_t LevelHashing<T>::S_HASH(T key) {
  if constexpr (std::is_pointer_v<T>) {
    return h(key->key, key->length, s_seed);
  } else {
    return h(&key, sizeof(Key_t), s_seed);
  }
}

template <class T>
void LevelHashing<T>::generate_seeds() {
  srand(time(NULL));
  do {
    f_seed = rand();
    s_seed = rand();
    f_seed = f_seed << (rand() % 63);
    s_seed = s_seed << (rand() % 63);
  } while (f_seed == s_seed);
}

template <class T>
LevelHashing<T>::LevelHashing(void) {}

template <class T>
LevelHashing<T>::~LevelHashing(void) {}

void *cache_align(void *ptr) {
  uint64_t pp = (uint64_t)ptr;
  pp += 48;
  return (void *)pp;
}

/* Initialize Function for Level Hashing*/
template <class T>
void initialize_level(PMEMobjpool *_pop, LevelHashing<T> *level, void *arg) {
  TX_BEGIN(_pop) {
    pmemobj_tx_add_range_direct(level, sizeof(*level));
    /* modify*/
    level->pop = _pop;
    level->levels = *((int *)arg);
    level->resize_num = 0;
    level->resizing_lock = 0;
    level->resizing = false;
    level->addr_capacity = pow(2, level->levels);
    level->total_capacity = pow(2, level->levels) + pow(2, level->levels - 1);
    level->locksize = 128;
    level->nlocks = (3 * level->addr_capacity / 2) / level->locksize + 1;

    level->generate_seeds();
    level->level_item_num[0] = 0;
    level->level_item_num[1] = 0;
    level->_interim_level_buckets = OID_NULL;
    level->interim_level_buckets = NULL;
    /*allocate*/
    level->_mutex = pmemobj_tx_zalloc(sizeof(PMEMrwlock) * level->nlocks,
                                      TOID_TYPE_NUM(char));
    level->_buckets[0] = pmemobj_tx_zalloc(
        sizeof(Node<T>) * level->addr_capacity + 1, TOID_TYPE_NUM(char));
    level->_buckets[1] = pmemobj_tx_zalloc(
        sizeof(Node<T>) * level->addr_capacity / 2 + 1, TOID_TYPE_NUM(char));
    level->_old_mutex = OID_NULL;

    /* Intialize pointer*/
    level->buckets[0] =
        (Node<T> *)cache_align(pmemobj_direct(level->_buckets[0]));
    level->buckets[1] =
        (Node<T> *)cache_align(pmemobj_direct(level->_buckets[1]));
    level->mutex = (PMEMrwlock *)pmemobj_direct(level->_mutex);
  }
  TX_END
}

template <class T>
void remapping(LevelHashing<T> *level) {
  level->buckets[0] =
      (Node<T> *)cache_align(pmemobj_direct(level->_buckets[0]));
  level->buckets[1] =
      (Node<T> *)cache_align(pmemobj_direct(level->_buckets[1]));
}

template <class T>
int LevelHashing<T>::Insert(T key, Value_t value) {
RETRY:
  while (resizing_lock.load() == 1) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);

#ifdef UNIQUE_CHECK
  uint32_t lock_idx = f_idx / locksize;
  while (pmemobj_rwlock_trywrlock(pop, &mutex[lock_idx]) != 0) {
    if (resizing == true) {
      goto RETRY;
    }
  }

  if (resizing == true) {
    pmemobj_rwlock_unlock(pop, &mutex[lock_idx]);
    goto RETRY;
  }

  bool inserted = false;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < ASSOC_NUM; ++j) {
      if constexpr (std::is_pointer_v<T>) {
        if ((buckets[i][f_idx].token[j] == 1) &&
            var_compare(buckets[i][f_idx].slot[j].key->key, key->key,
                        buckets[i][f_idx].slot[j].key->length, key->length)) {
          inserted = true;
          goto UNIQUE;
        }
      } else {
        if ((buckets[i][f_idx].token[j] == 1) &&
            (buckets[i][f_idx].slot[j].key == key)) {
          inserted = true;
          goto UNIQUE;
        }
      }
    }

    for (int j = 0; j < ASSOC_NUM; ++j) {
      if ((buckets[i][s_idx].token[j] == 1) &&
          (buckets[i][s_idx].slot[j].key == key)) {
        inserted = true;
        goto UNIQUE;
      }
    }

    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }

UNIQUE:
  if (inserted) {
    pmemobj_rwlock_unlock(pop, &mutex[lock_idx]);
    return -1; // unique check failure
  }

  pmemobj_rwlock_unlock(pop, &mutex[lock_idx]);
  f_idx = F_IDX(f_hash, addr_capacity);
  s_idx = S_IDX(s_hash, addr_capacity);
#endif

  int i, j;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < ASSOC_NUM; j++) {
      while (pmemobj_rwlock_trywrlock(pop, &mutex[f_idx / locksize]) != 0) {
        if (resizing == true) {
          goto RETRY;
        }
      }

      if (resizing == true) {
        pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
        goto RETRY;
      }
      if (buckets[i][f_idx].token[j] == 0) {
        buckets[i][f_idx].slot[j].value = value;
        buckets[i][f_idx].slot[j].key = key;
        pmemobj_persist(pop, &buckets[i][f_idx].slot[j], sizeof(Entry<T>));
        buckets[i][f_idx].token[j] = 1;
        pmemobj_persist(pop, &buckets[i][f_idx].token[j], sizeof(uint8_t));
#ifdef COUNTING
        level_item_num[i]++;
#endif
        pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
        return 0;
      }

      pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);

      while (pmemobj_rwlock_trywrlock(pop, &mutex[s_idx / locksize]) != 0) {
        if (resizing == true) {
          goto RETRY;
        }
      }

      if (resizing == true) {
        pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
        goto RETRY;
      }

      if (buckets[i][s_idx].token[j] == 0) {
        buckets[i][s_idx].slot[j].value = value;
        buckets[i][s_idx].slot[j].key = key;
        pmemobj_persist(pop, &buckets[i][s_idx].slot[j], sizeof(Entry<T>));
        buckets[i][s_idx].token[j] = 1;
        pmemobj_persist(pop, &buckets[i][s_idx].token[j], sizeof(uint8_t));
#ifdef COUNTING
        level_item_num[i]++;
#endif
        pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
        return 0;
      }
      pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
    }

    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }

  f_idx = F_IDX(f_hash, addr_capacity);
  s_idx = S_IDX(s_hash, addr_capacity);
  int empty_loc;
  int64_t lock = 0;
  if (resizing_lock.compare_exchange_strong(lock, 1)) {
    for (i = 0; i < 2; i++) {
      if (!try_movement(pop, f_idx, i, key, value)) {
        resizing_lock.store(0);
        return 0;
      }
      if (!try_movement(pop, s_idx, i, key, value)) {
        resizing_lock.store(0);
        return 0;
      }
      f_idx = F_IDX(f_hash, addr_capacity / 2);
      s_idx = S_IDX(s_hash, addr_capacity / 2);
    }

    if (resize_num > 0) {
      {
        pmemobj_rwlock_wrlock(pop, &mutex[f_idx / locksize]);
        empty_loc = b2t_movement(pop, f_idx);
        if (empty_loc != -1) {
          buckets[1][f_idx].slot[empty_loc].value = value;
          buckets[1][f_idx].slot[empty_loc].key = key;
          pmemobj_persist(pop, &buckets[1][f_idx].slot[empty_loc],
                          sizeof(Entry<T>));
          buckets[1][f_idx].token[empty_loc] = 1;
          pmemobj_persist(pop, &buckets[1][f_idx].token[empty_loc],
                          sizeof(uint8_t));
#ifdef COUNTING
          level_item_num[1]++;
#endif
          resizing_lock.store(0);
          pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
          return 0;
        }
        pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
      }
      {
        pmemobj_rwlock_wrlock(pop, &mutex[s_idx / locksize]);
        empty_loc = b2t_movement(pop, s_idx);
        if (empty_loc != -1) {
          buckets[1][s_idx].slot[empty_loc].value = value;
          buckets[1][s_idx].slot[empty_loc].key = key;
          pmemobj_persist(pop, &buckets[1][s_idx].slot[empty_loc],
                          sizeof(Entry<T>));
          buckets[1][s_idx].token[empty_loc] = 1;
          pmemobj_persist(pop, &buckets[1][s_idx].token[empty_loc],
                          sizeof(uint8_t));
#ifdef COUNTING
          level_item_num[1]++;
#endif
          resizing_lock.store(0);
          pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
          return 0;
        }
        pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
      }
    }
    resize(pop);
    resizing_lock.store(0);
  }
  goto RETRY;
}

template <class T>
void LevelHashing<T>::resize(PMEMobjpool *pop) {
  std::cout << "Resizing towards levels " << levels + 1 << std::endl;
  resizing = true;
  for (int i = 0; i < nlocks; ++i) {
    pmemobj_rwlock_wrlock(pop, &mutex[i]);
  }

  size_t new_addr_capacity = pow(2, levels + 1);
  _old_mutex = _mutex;
  pmemobj_persist(pop, &_old_mutex, sizeof(_old_mutex));
  nlocks = (3 * 2 * addr_capacity / 2) / locksize + 1;
  auto ret = pmemobj_zalloc(pop, &_mutex, nlocks * sizeof(PMEMrwlock),
                            TOID_TYPE_NUM(char));
  if (ret) {
    std::cout << "Ret = " << ret << std::endl;
    std::cout << "Allocation Size = " << std::hex << nlocks * sizeof(PMEMrwlock)
              << std::endl;
    LOG_FATAL("Allocation Error in New Mutex");
  }
  mutex = (PMEMrwlock *)(pmemobj_direct(_mutex));
  ret = pmemobj_zalloc(pop, &_interim_level_buckets,
                       new_addr_capacity * sizeof(Node<T>) + 1,
                       TOID_TYPE_NUM(char));
  if (ret) {
    std::cout << "Ret = " << ret << std::endl;
    std::cout << "Allocation Size = " << std::hex
              << new_addr_capacity * sizeof(Node<T>) + 1 << std::endl;
    LOG_FATAL("Allocation Error in Table");
  }
  interim_level_buckets =
      (Node<T> *)cache_align(pmemobj_direct(_interim_level_buckets));

  uint64_t new_level_item_num = 0;
  uint64_t old_idx;
  for (old_idx = 0; old_idx < pow(2, levels - 1); old_idx++) {
    uint64_t i, j;
    for (i = 0; i < ASSOC_NUM; i++) {
      if (buckets[1][old_idx].token[i] == 1) {
        T key = buckets[1][old_idx].slot[i].key;
        Value_t value = buckets[1][old_idx].slot[i].value;

        uint32_t f_idx = F_IDX(F_HASH(key), new_addr_capacity);
        uint32_t s_idx = S_IDX(S_HASH(key), new_addr_capacity);

        uint8_t insertSuccess = 0;
        for (j = 0; j < ASSOC_NUM; j++) {
          if (interim_level_buckets[f_idx].token[j] == 0) {
            interim_level_buckets[f_idx].slot[j].value = value;
            interim_level_buckets[f_idx].slot[j].key = key;
#ifndef BATCH
            pmemobj_persist(pop, &interim_level_buckets[f_idx].slot[j],
                            sizeof(Entry<T>));
#endif
            interim_level_buckets[f_idx].token[j] = 1;
#ifndef BATCH
            pmemobj_persist(pop, &interim_level_buckets[f_idx].token[j],
                            sizeof(uint8_t));
#endif
            insertSuccess = 1;
#ifdef COUNTING
            new_level_item_num++;
#endif
            break;
          } else if (interim_level_buckets[s_idx].token[j] == 0) {
            interim_level_buckets[s_idx].slot[j].value = value;
            interim_level_buckets[s_idx].slot[j].key = key;
#ifndef BATCH
            pmemobj_persist(pop, &interim_level_buckets[s_idx].slot[j],
                            sizeof(Entry<T>));
#endif
            interim_level_buckets[s_idx].token[j] = 1;
#ifndef BATCH
            pmemobj_persist(pop, &interim_level_buckets[s_idx].token[j],
                            sizeof(uint8_t));
#endif
            insertSuccess = 1;
#ifdef COUNTING
            new_level_item_num++;
#endif
            break;
          }
        }

#ifndef BATCH
        buckets[1][old_idx].token[i] = 0;
        pmemobj_persist(pop, &buckets[1][old_idx].token[i], sizeof(uint8_t));
#endif
      }
    }
  }

#ifdef BATCH
  pmemobj_persist(pop, &buckets[1][0], sizeof(Node<T>) * pow(2, levels - 1));
  pmemobj_persist(pop, &interim_level_buckets[0],
                  sizeof(Node<T>) * new_addr_capacity);
#endif

  TX_BEGIN(pop) {
    pmemobj_tx_add_range_direct(this, sizeof(class LevelHashing));

    levels++;
    resize_num++;
    pmemobj_tx_free(_buckets[1]); /*free the old bottom level*/
    _buckets[1] = _buckets[0];
    _buckets[0] = _interim_level_buckets;
    buckets[1] = (Node<T> *)cache_align(pmemobj_direct(_buckets[1]));
    buckets[0] = (Node<T> *)cache_align(pmemobj_direct(_buckets[0]));

    _interim_level_buckets = OID_NULL;
    interim_level_buckets = NULL;
#ifdef COUNTING
    level_item_num[1] = level_item_num[0];
    level_item_num[0] = new_level_item_num;
#endif
    addr_capacity = new_addr_capacity;
    total_capacity = pow(2, levels) + pow(2, levels - 1);
    pmemobj_tx_free(_old_mutex); /*free the old mutex*/
    resizing = false;
  }
  TX_ONABORT { printf("resizing txn 2 fails\n"); }
  TX_END
  std::cout << "Done! :Resizing towards levels " << levels << std::endl;
}

template <class T>
uint8_t LevelHashing<T>::try_movement(PMEMobjpool *pop, uint64_t idx,
                                      uint64_t level_num, T key,
                                      Value_t value) {
  uint64_t i, j, jdx;
  {
    pmemobj_rwlock_wrlock(pop, &mutex[idx / locksize]);
    for (i = 0; i < ASSOC_NUM; i++) {
      T m_key = buckets[level_num][idx].slot[i].key;
      Value_t m_value = buckets[level_num][idx].slot[i].value;
      uint64_t f_hash = F_HASH(m_key);
      uint64_t s_hash = S_HASH(m_key);
      uint64_t f_idx = F_IDX(f_hash, addr_capacity / (1 + level_num));
      uint64_t s_idx = S_IDX(s_hash, addr_capacity / (1 + level_num));

      if (f_idx == idx)
        jdx = s_idx;
      else
        jdx = f_idx;

      if ((jdx / locksize) != (idx / locksize)) {
        pmemobj_rwlock_wrlock(pop, &mutex[jdx / locksize]);
      }

      for (j = 0; j < ASSOC_NUM; j++) {
        if (buckets[level_num][jdx].token[j] == 0) {
          buckets[level_num][jdx].slot[j].value = m_value;
          buckets[level_num][jdx].slot[j].key = m_key;
          pmemobj_persist(pop, &buckets[level_num][jdx].slot[j],
                          sizeof(Entry<T>));
          buckets[level_num][jdx].token[j] = 1;
          pmemobj_persist(pop, &buckets[level_num][jdx].token[j],
                          sizeof(uint8_t));
          buckets[level_num][idx].token[i] = 0;
          pmemobj_persist(pop, &buckets[level_num][idx].token[i],
                          sizeof(uint8_t));

          buckets[level_num][idx].slot[i].value = value;
          buckets[level_num][idx].slot[i].key = key;
          pmemobj_persist(pop, &buckets[level_num][idx].slot[i],
                          sizeof(Entry<T>));
          buckets[level_num][idx].token[i] = 1;
          pmemobj_persist(pop, &buckets[level_num][idx].token[i],
                          sizeof(uint8_t));
#ifdef COUNTING
          level_item_num[level_num]++;
#endif
          if ((jdx / locksize) != (idx / locksize)) {
            pmemobj_rwlock_unlock(pop, &mutex[jdx / locksize]);
          }
          pmemobj_rwlock_unlock(pop, &mutex[idx / locksize]);
          return 0;
        }
      }
      if ((jdx / locksize) != (idx / locksize))
        pmemobj_rwlock_unlock(pop, &mutex[jdx / locksize]);
    }
    pmemobj_rwlock_unlock(pop, &mutex[idx / locksize]);
  }
  return 1;
}

template <class T>
int LevelHashing<T>::b2t_movement(PMEMobjpool *pop, uint64_t idx) {
  T key;
  Value_t value;
  uint64_t s_hash, f_hash;
  uint64_t s_idx, f_idx;
  uint64_t i, j;

  for (i = 0; i < ASSOC_NUM; i++) {
    key = buckets[1][idx].slot[i].key;
    value = buckets[1][idx].slot[i].value;
    f_hash = F_HASH(key);
    s_hash = S_HASH(key);
    f_idx = F_IDX(f_hash, addr_capacity);
    s_idx = S_IDX(s_hash, addr_capacity);

    for (j = 0; j < ASSOC_NUM; j++) {
      if ((idx / locksize) != (f_idx / locksize))
        pmemobj_rwlock_wrlock(pop, &mutex[f_idx / locksize]);

      if (buckets[0][f_idx].token[j] == 0) {
        buckets[0][f_idx].slot[j].value = value;
        buckets[0][f_idx].slot[j].key = key;
        pmemobj_persist(pop, &buckets[0][f_idx].slot[j], sizeof(Entry<T>));
        buckets[0][f_idx].token[j] = 1;
        pmemobj_persist(pop, &buckets[0][f_idx].token[j], sizeof(uint8_t));
        buckets[1][idx].token[i] = 0;
        pmemobj_persist(pop, &buckets[1][idx].token[i], sizeof(uint8_t));
#ifdef COUNTING
        level_item_num[0]++;
        level_item_num[1]--;
#endif
        if ((idx / locksize) != (f_idx / locksize))
          pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
        return i;
      }
      if ((idx / locksize) != (f_idx / locksize))
        pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
      if ((idx / locksize) != (s_idx / locksize))
        pmemobj_rwlock_wrlock(pop, &mutex[s_idx / locksize]);

      if (buckets[0][s_idx].token[j] == 0) {
        buckets[0][s_idx].slot[j].value = value;
        buckets[0][s_idx].slot[j].key = key;
        pmemobj_persist(pop, &buckets[0][s_idx].slot[j], sizeof(Entry<T>));
        buckets[0][s_idx].token[j] = 1;
        pmemobj_persist(pop, &buckets[0][s_idx].token[j], sizeof(uint8_t));
        buckets[1][idx].token[i] = 0;
        pmemobj_persist(pop, &buckets[0][s_idx].token[j], sizeof(uint8_t));
#ifdef COUNTING
        level_item_num[0]++;
        level_item_num[1]--;
#endif
        if ((idx / locksize) != (s_idx / locksize))
          pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
        return i;
      }
      if ((idx / locksize) != (s_idx / locksize))
        pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
    }
  }
  return -1;
}

template <class T>
Value_t LevelHashing<T>::Get(T key) {
RETRY:
  while (resizing == true) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);
  int i = 0, j;

  for (i = 0; i < 2; i++) {
    {
      while (pmemobj_rwlock_tryrdlock(pop, &mutex[f_idx / locksize]) != 0) {
        if (resizing == true) {
          goto RETRY;
        }
      }

      if (resizing == true) {
        pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
        goto RETRY;
      }

      for (j = 0; j < ASSOC_NUM; j++) {
        if constexpr (std::is_pointer_v<T>) {
          if ((buckets[i][f_idx].token[j] == 1) &&
              var_compare(buckets[i][f_idx].slot[j].key->key, key->key,
                          buckets[i][f_idx].slot[j].key->length, key->length)) {
            pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
            return buckets[i][f_idx].slot[j].value;
          }
        } else {
          if (buckets[i][f_idx].token[j] == 1 &&
              buckets[i][f_idx].slot[j].key == key) {
            pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
            return buckets[i][f_idx].slot[j].value;
          }
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
    }
    {
      while (pmemobj_rwlock_tryrdlock(pop, &mutex[s_idx / locksize]) != 0) {
        if (resizing == true) {
          goto RETRY;
        }
      }

      if (resizing == true) {
        pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
        goto RETRY;
      }

      for (j = 0; j < ASSOC_NUM; j++) {
        if constexpr (std::is_pointer_v<T>) {
          if ((buckets[i][s_idx].token[j] == 1) &&
              var_compare(buckets[i][s_idx].slot[j].key->key, key->key,
                          buckets[i][s_idx].slot[j].key->length, key->length)) {
            pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
            return buckets[i][s_idx].slot[j].value;
          }
        } else {
          if (buckets[i][s_idx].token[j] == 1 &&
              buckets[i][s_idx].slot[j].key == key) {
            pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
            return buckets[i][s_idx].slot[j].value;
          }
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
    }
    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }
  return NONE;
}

template <class T>
bool LevelHashing<T>::Delete(T key) {
RETRY:
  while (resizing == true) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);
  int i = 0, j;

  for (i = 0; i < 2; i++) {
    {
      while (pmemobj_rwlock_trywrlock(pop, &mutex[f_idx / locksize]) != 0) {
        if (resizing == true) {
          goto RETRY;
        }
      }

      if (resizing == true) {
        pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
        goto RETRY;
      }

      for (j = 0; j < ASSOC_NUM; j++) {
        if constexpr (std::is_pointer_v<T>) {
          if ((buckets[i][f_idx].token[j] == 1) &&
              var_compare(buckets[i][f_idx].slot[j].key->key, key->key,
                          buckets[i][f_idx].slot[j].key->length, key->length)) {
            buckets[i][f_idx].token[j] = 0;
            pmemobj_persist(pop, &buckets[i][f_idx].token[j], sizeof(uint8_t));
            pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
            return true;
          }
        } else {
          if (buckets[i][f_idx].token[j] == 1 &&
              buckets[i][f_idx].slot[j].key == key) {
            buckets[i][f_idx].token[j] = 0;
            pmemobj_persist(pop, &buckets[i][f_idx].token[j], sizeof(uint8_t));
            pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
            return true;
          }
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[f_idx / locksize]);
    }
    {
      while (pmemobj_rwlock_trywrlock(pop, &mutex[s_idx / locksize]) != 0) {
        if (resizing == true) {
          goto RETRY;
        }
      }

      if (resizing == true) {
        pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
        goto RETRY;
      }

      for (j = 0; j < ASSOC_NUM; j++) {
        if constexpr (std::is_pointer_v<T>) {
          if ((buckets[i][s_idx].token[j] == 1) &&
              var_compare(buckets[i][s_idx].slot[j].key->key, key->key,
                          buckets[i][s_idx].slot[j].key->length, key->length)) {
            buckets[i][s_idx].token[j] = 0;
            pmemobj_persist(pop, &buckets[i][s_idx].token[j], sizeof(uint8_t));
            pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
            return true;
          }
        } else {
          if (buckets[i][s_idx].token[j] == 1 &&
              buckets[i][s_idx].slot[j].key == key) {
            buckets[i][s_idx].token[j] = 0;
            pmemobj_persist(pop, &buckets[i][s_idx].token[j], sizeof(uint8_t));
            pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
            return true;
          }
        }
      }
      pmemobj_rwlock_unlock(pop, &mutex[s_idx / locksize]);
    }
    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }
  return false;
}
}  // namespace level
#endif  // LEVEL_HASHING_H_