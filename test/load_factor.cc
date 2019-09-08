/*
* This program is to test the load factor of the sub-hash-table
*/
#include <immintrin.h>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <shared_mutex>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "./util/hash.h"
#include "./util/pair.h"
#include "./util/persist.h"
#include "./util/random.h"

#define _INVALID 0 /* we use 0 as the invalid key*/
#define SINGLE 1
#define COUNTING 1

#define SIMD 1
#define SIMD_CMP8(src, key)                                         \
  do {                                                              \
    const __m256i key_data = _mm256_set1_epi8(key);                 \
    __m256i seg_data =                                              \
        _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src)); \
    __m256i rv_mask = _mm256_cmpeq_epi8(seg_data, key_data);        \
    mask = _mm256_movemask_epi8(rv_mask);                           \
  } while (0)

#define SSE_CMP8(src, key)                                       \
  do {                                                           \
    const __m128i key_data = _mm_set1_epi8(key);                 \
    __m128i seg_data =                                           \
        _mm_loadu_si128(reinterpret_cast<const __m128i *>(src)); \
    __m128i rv_mask = _mm_cmpeq_epi8(seg_data, key_data);        \
    mask = _mm_movemask_epi8(rv_mask);                           \
  } while (0)

#define CHECK_BIT(var, pos) ((((var) & (1 << pos)) > 0) ? (1) : (0))

const uint32_t lockSet = ((uint64_t)1 << 31);      /*locking information*/
const uint32_t lockMask = ((uint64_t)1 << 31) - 1; /*locking mask*/
const int overflowSet = 1 << 15;
const int countMask = (1 << 4) - 1;

template<class T>
struct _Pair {
  T key;
  Value_t value;
};

constexpr size_t k_PairSize = 16;  // a k-v _Pair with a bit
constexpr size_t kNumPairPerBucket = 14; /* it is determined by the usage of the fingerprint*/
constexpr size_t kFingerBits = 8;
constexpr size_t kMask = (1 << kFingerBits) - 1;
//const constexpr size_t kNumBucket = 16;
//constexpr size_t stashBucket = 2;
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
//constexpr size_t bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
//constexpr size_t stashMask = (1 << (int)log2(stashBucket)) - 1;
//constexpr uint8_t stashHighMask = ~((uint8_t)stashMask);
constexpr uint8_t preNeighborSet = 1 << 7;
constexpr uint8_t nextNeighborSet = 1 << 6;
//#define BUCKET_INDEX(hash) ((hash >> kFingerBits) % kNumBucket)
#define GET_COUNT(var) ((var)&countMask)
#define GET_BITMAP(var) (((var) >> 4) & allocMask)
#define ORG_BITMAP(var) ((~((var)&allocMask)) & allocMask)
#define PROBE_BITMAP(var) ((var)&allocMask)

template<class T>
struct Bucket {
  inline int find_empty_slot() {
    if (GET_COUNT(bitmap) == kNumPairPerBucket) {
      return -1;
    }
    auto mask = ~(GET_BITMAP(bitmap));  // Now the 1 bit should be the empty
                                        // slot
    return __builtin_ctz(mask);
  }

  /*true indicates overflow, needs extra check in the stash*/
  inline bool test_overflow() { return (overflowCount != 0) ? true : false; }

  inline bool test_stash_check() {
    int mask = *((int *)membership);
    return ((mask & overflowSet) != 0) ? true : false;
  }

  inline void clear_stash_check() {
    int mask = *((int *)membership);
    *((int *)membership) = (*((int *)membership)) & (~overflowSet);
  }

  inline void set_indicator(uint8_t meta_hash, Bucket<T> *neighbor, uint8_t pos) {
    int mask = finger_array[14];
    mask = ~mask;
    auto index = __builtin_ctz(mask);

    if (index < 4) {
      finger_array[15 + index] = meta_hash;
      finger_array[14] =
          ((uint8_t)(1 << index) | finger_array[14]); /*may be optimized*/
      finger_array[19] =
          (finger_array[19] & (~(3 << (index * 2)))) | (pos << (index * 2));
    } else {
      mask = neighbor->finger_array[14];
      mask = ~mask;
      index = __builtin_ctz(mask);
      if (index < 4) {
        neighbor->finger_array[15 + index] = meta_hash;
        neighbor->finger_array[14] =
            ((uint8_t)(1 << index) | neighbor->finger_array[14]);
        neighbor->overflowMember =
            ((uint8_t)(1 << index) | neighbor->overflowMember);
        neighbor->finger_array[19] =
            (neighbor->finger_array[19] & (~(3 << (index * 2)))) |
            (pos << (index * 2));
      } else { /*overflow, increase count*/
        overflowCount++;
      }
    }
    *((int *)membership) = (*((int *)membership)) | overflowSet;
  }

  /*both clear this bucket and its neighbor bucket*/
  /*
  inline void unset_indicator(uint8_t meta_hash, Bucket<T> *neighbor, T key,
                              uint64_t pos) {
    bool clear_success = false;
    int mask1 = finger_array[14];
    for (int i = 0; i < 4; ++i) {
      if (CHECK_BIT(mask1, i) && (finger_array[15 + i] == meta_hash) &&
          (((1 << i) & overflowMember) == 0) &&
          (((finger_array[19] >> (2 * i)) & stashMask) == pos)) {
        // printf("clear the indicator 1\n");
        finger_array[14] = finger_array[14] & ((uint8_t)(~(1 << i)));
        finger_array[19] = finger_array[19] & (~(3 << (i * 2)));
        assert(((finger_array[19] >> (i * 2)) & stashMask) == 0);
        clear_success = true;
        break;
      }
    }

    int mask2 = neighbor->finger_array[14];
    if (!clear_success) {
      for (int i = 0; i < 4; ++i) {
        if (CHECK_BIT(mask2, i) &&
            (neighbor->finger_array[15 + i] == meta_hash) &&
            (((1 << i) & neighbor->overflowMember) != 0) &&
            (((neighbor->finger_array[19] >> (2 * i)) & stashMask) == pos)) {
          // printf("clear the indicator 2\n");
          neighbor->finger_array[14] =
              neighbor->finger_array[14] & ((uint8_t)(~(1 << i)));
          neighbor->overflowMember =
              neighbor->overflowMember & ((uint8_t)(~(1 << i)));
          neighbor->finger_array[19] =
              neighbor->finger_array[19] & (~(3 << (i * 2)));
          assert(((neighbor->finger_array[19] >> (i * 2)) & stashMask) == 0);
          clear_success = true;
          break;
        }
      }
    }

    if (!clear_success) {
      overflowCount--;
    }

    mask1 = finger_array[14];
    mask2 = neighbor->finger_array[14];
    if (((mask1 & (~overflowMember)) == 0) && (overflowCount == 0) &&
        ((mask2 & neighbor->overflowMember) == 0)) {
      clear_stash_check();
    }
  }
  */

  int unique_check(uint8_t meta_hash, T key, Bucket<T> *neighbor,
                   Bucket<T> *stash, size_t stashBucket) {
    // return check_and_get(meta_hash, key) == NONE ? 0 : -1;
    if ((check_and_get(meta_hash, key, false) != NONE) ||
        (neighbor->check_and_get(meta_hash, key, true) != NONE)) {
      return -1;
    }

    if (test_stash_check()) {
      auto test_stash = false;
      if (test_overflow()) {
        test_stash = true;
      } else {
        int mask = finger_array[14];
        if (finger_array[14] != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) && (finger_array[15 + i] == meta_hash) &&
                (((1 << i) & overflowMember) == 0)) {
              test_stash = true;
              goto STASH_CHECK;
            }
          }
        }

        if (neighbor->finger_array[14] != 0) {
          mask = neighbor->finger_array[14];
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (neighbor->finger_array[15 + i] == meta_hash) &&
                (((1 << i) & neighbor->overflowMember) != 0)) {
              test_stash = true;
              break;
            }
          }
        }
      }
    STASH_CHECK:
      if (test_stash == true) {
        for (int i = 0; i < stashBucket; ++i) {
          Bucket *curr_bucket = stash + i;
          if (curr_bucket->check_and_get(meta_hash, key, false) != NONE) {
            return -1;
          }
        }
      }
    }
    return 0;
  }

  Value_t check_and_get(uint8_t meta_hash, T key, bool probe) {
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & ((~(*(int *)membership)) & allocMask);
    } else {
      mask = mask & GET_BITMAP(bitmap) & ((*(int *)membership) & allocMask);
    }

    if constexpr (std::is_pointer_v<T>){
        /*loop unrolling*/
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (strcmp(_[i].key, key) == 0)) {
            return _[i].value;
          }

          if (CHECK_BIT(mask, i + 1) && (strcmp(_[i + 1].key, key) == 0)) {
            return _[i + 1].value;
          }

          if (CHECK_BIT(mask, i + 2) && (strcmp(_[i + 2].key, key) == 0)) {
            return _[i + 2].value;
          }

          if (CHECK_BIT(mask, i + 3) && (strcmp(_[i + 3].key, key) == 0)) {
            return _[i + 3].value;
          }
        }

        if (CHECK_BIT(mask, 12) && (strcmp(_[12].key, key) == 0)) {
          return _[12].value;
        }

        if (CHECK_BIT(mask, 13) && (strcmp(_[13].key, key) == 0)) {
          return _[13].value;
        }
      }
    }else{
      /*loop unrolling*/
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (_[i].key == key)) {
            return _[i].value;
          }

          if (CHECK_BIT(mask, i + 1) && (_[i + 1].key == key)) {
            return _[i + 1].value;
          }

          if (CHECK_BIT(mask, i + 2) && (_[i + 2].key == key)) {
            return _[i + 2].value;
          }

          if (CHECK_BIT(mask, i + 3) && (_[i + 3].key == key)) {
            return _[i + 3].value;
          }
        }

        if (CHECK_BIT(mask, 12) && (_[12].key == key)) {
          return _[12].value;
        }

        if (CHECK_BIT(mask, 13) && (_[13].key == key)) {
          return _[13].value;
        }
      }
    }
    return NONE;
  }

  inline void set_hash(int index, uint8_t meta_hash,
                       bool probe) /* Do I needs the atomic instruction????*/
  {
    finger_array[index] = meta_hash;
    auto new_bitmap = bitmap | (1 << (index + 4));
    assert(GET_COUNT(bitmap) < kNumPairPerBucket);
    new_bitmap += 1;
    bitmap = new_bitmap;
    if (probe) {
      *((int *)membership) = (1 << index) | *((int *)membership);
    }
  }

  inline uint8_t get_hash(int index) { return finger_array[index]; }

  inline void unset_hash(int index, bool nt_flush = false) {
    auto new_bitmap = bitmap & (~(1 << (index + 4)));
    assert(GET_COUNT(bitmap) <= kNumPairPerBucket);
    assert(GET_COUNT(bitmap) > 0);
    new_bitmap -= 1;
    bitmap = new_bitmap;

    *((int *)membership) =
        (~(1 << index)) &
        (*((int *)membership)); /*since they are in the same cacheline,
                                   therefore no performance influence?*/
  }

  inline void get_lock() {
    auto old_value = version_lock & lockMask;
    auto new_value = version_lock | lockSet;
    while (!CAS(&version_lock, &old_value, new_value)) {
      old_value = version_lock & lockMask;
      new_value = version_lock | lockSet;
    }
  }

  inline bool try_get_lock() {
    auto old_value = version_lock & lockMask;
    auto new_value = version_lock | lockSet;
    return CAS(&version_lock, &old_value, new_value);
  }

  inline void release_lock() {
    auto old_value = version_lock;
    auto new_value = ((old_value & lockMask) + 1) & lockMask;

    while (!CAS(&version_lock, &old_value, new_value)) {
      old_value = version_lock;
      new_value = ((old_value & lockMask) + 1) & lockMask;
    }
  }

  /*if the lock is set, return true*/
  inline bool test_lock_set(uint32_t &version) {
    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    version = value & lockMask;
    return (value & lockSet) != 0;
  }

  // test whether the version has change, if change, return true
  inline bool test_lock_version_change(uint32_t old_version) {
    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    auto version = value & lockMask;
    return ((value & lockSet) != 0) || (version != old_version);
  }

  int Insert(T key, Value_t value, uint8_t meta_hash, bool probe) {
    auto slot = find_empty_slot();
    /* this branch can be optimized out*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      printf("cannot find the empty slot, for key %llu\n", key);
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
    mfence();
    set_hash(slot, meta_hash, probe);
    return 0;
  }

  /*if delete success, then return 0, else return -1*/
  int Delete(T key, uint8_t meta_hash, bool probe) {
    /*do the simd and check the key, then do the delete operation*/
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & ((~(*(int *)membership)) & allocMask);
    } else {
      mask = mask & GET_BITMAP(bitmap) & ((*(int *)membership) & allocMask);
    }
    /*loop unrolling*/
    if constexpr (std::is_pointer_v<T>){
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (strcmp(_[i].key, key) == 0)) {
            unset_hash(i, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) && (strcmp(_[i + 1].key, key) == 0)) {
            unset_hash(i + 1, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) && (strcmp(_[i + 2].key, key) == 0)) {
            unset_hash(i + 2, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) && (strcmp(_[i + 3].key, key) == 0)) {
            unset_hash(i + 3, false);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) && (strcmp(_[12].key, key) == 0)) {
          unset_hash(12, false);
          return 0;
        }

        if (CHECK_BIT(mask, 13) && (strcmp(_[13].key, key) == 0)) {
          unset_hash(13, false);
          return 0;
        }
      }
    }else{
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (_[i].key == key)) {
            unset_hash(i, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) && (_[i + 1].key == key)) {
            unset_hash(i + 1, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) && (_[i + 2].key == key)) {
            unset_hash(i + 2, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) && (_[i + 3].key == key)) {
            unset_hash(i + 3, false);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) && (_[12].key == key)) {
          unset_hash(12, false);
          return 0;
        }

        if (CHECK_BIT(mask, 13) && (_[13].key == key)) {
          unset_hash(13, false);
          return 0;
        }
      }
    }
    return -1;
}

  int Insert_with_noflush(T key, Value_t value, uint8_t meta_hash,
                          bool probe) {
    auto slot = find_empty_slot();
    /* this branch can be optimized out*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      printf("cannot find the empty slot, for key %llu\n", key);
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
    return 0;
  }

  void Insert_displace(T key, Value_t value, uint8_t meta_hash, int slot,
                       bool probe) {
    _[slot].value = value;
    _[slot].key = key;
    mfence();
    set_hash(slot, meta_hash, probe);
  }

  void Insert_displace_with_noflush(T key, Value_t value, uint8_t meta_hash,
                                    int slot, bool probe) {
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
  }

  inline int Find_org_displacement() {
    int mask = (~(*((int *)membership))) & allocMask;
    if (mask == 0) {
      return -1;
    }
    return __builtin_ctz(mask);
  }

  /*find element that it is in the probe*/
  inline int Find_probe_displacement() {
    int mask = (*((int *)membership)) & allocMask;
    if (mask == 0) {
      return -1;
    }
    return __builtin_ctz(mask);
  }

  	/*suite of function to set/unset the flush state of the node*/
	inline void setPreNonFlush(){
    overflowMember = overflowMember | preNeighborSet;
	}

	inline void unsetPreNonFlush(){
    overflowMember = overflowMember & (~preNeighborSet);
	}

	/* If the non-flush is set, return true*/
	inline bool testPreNonFlush(){
		return ((overflowMember & preNeighborSet) != 0);
	}

	inline void setNextNonFlush(){
		overflowMember = overflowMember | nextNeighborSet;
	}

	inline void unsetNextNonFlush(){
		overflowMember = overflowMember & (~nextNeighborSet);
	}

	inline bool testNextNonFlush(){
		return ((overflowMember & nextNeighborSet) != 0);
	}

  uint32_t version_lock;
  int bitmap;               // allocation bitmap + pointer bitmao + counter
  uint8_t finger_array[20]; /*only use the first 14 bytes, can be accelerated by
                               SSE instruction,0-13 for finger, 14-17 for
                               overflowed, 18 as the bitmap, 19 as the overflow bucket index indicator*/
  uint8_t membership[2];    /*Used to test whether the key originally belongs to
                               this bucket*/
  uint8_t overflowMember; /*overflowmember indicates membership of the overflow
                             fingerprint*/
  uint8_t overflowCount;

  _Pair<T> _[kNumPairPerBucket];
};

template<class T>
struct Table {
  Table(size_t _kNumBucket, size_t _stashBucket){
      kNumBucket = _kNumBucket;
      stashBucket = _stashBucket;
      bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
      stashMask = (1 << (int)log2(stashBucket)) - 1;
      stashHighMask = ~((uint8_t)stashMask);
      bucket = new Bucket<T>[kNumBucket + stashBucket];
      memset(bucket, 0, sizeof(Bucket<T>)*(kNumBucket + stashBucket));
  }

  ~Table(void) {
      delete [] bucket;
  }

  void displayLoadFactor(){
      double load_factor = (double)number / (14 * (kNumBucket + stashBucket));
      double raw_space = (double)number / (16 * (kNumBucket + stashBucket));
      std::cout << "The inserted number is " << number << std::endl;
      std::cout << "The load factor is " << load_factor << std::endl;
      std::cout << "The raw space utlization is " << raw_space << std::endl;
  }

  double getLoadFactor(){
      double load_factor = (double)number / (14 * (kNumBucket + stashBucket));
      return load_factor;
  }

  double getRawSpace(){
      double raw_space = (double)(number) / (16 * (kNumBucket + stashBucket));
      return raw_space;
  }

  uint64_t getNumber(){
    return number;
  }

  void reset(){
      memset(bucket, 0, sizeof(Bucket<T>)*(kNumBucket + stashBucket));
      number = 0;
  }

  size_t BUCKET_INDEX(size_t hash) {
      return ((hash >> kFingerBits) % kNumBucket);
  }

  int Insert(T key, Value_t value);
  int Next_displace(Bucket<T>* target, Bucket<T> *neighbor, Bucket<T> *next_neighbor, T key,
                    Value_t value, uint8_t meta_hash) {
    int displace_index = neighbor->Find_org_displacement();
    if ((GET_COUNT(next_neighbor->bitmap) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      // printf("do the displacement in next bucket, the displaced key is %lld,
      // the new key is %lld\n", neighbor->_[displace_index].key, key);
      next_neighbor->Insert(neighbor->_[displace_index].key,
                            neighbor->_[displace_index].value,
                            neighbor->finger_array[displace_index], true);
      //neighbor->setNextNonFlush();
      next_neighbor->release_lock();
      //neighbor->unsetNextNonFlush();
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace(key, value, meta_hash, displace_index, true);
      //target->setNextNonFlush();
      neighbor->release_lock();
      //target->unsetNextNonFlush();
      target->release_lock();
#ifdef COUNTING
      __sync_fetch_and_add(&number, 1);
#endif
      return 0;
    }
    return -1;
  }

  int Prev_displace(Bucket<T> *target, Bucket<T> *prev_neighbor, Bucket<T> *neighbor, T key,
                    Value_t value, uint8_t meta_hash) {
    int displace_index = target->Find_probe_displacement();
    if ((GET_COUNT(prev_neighbor->bitmap) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      // printf("do the displacement in previous bucket,the displaced key is
      // %lld, the new key is %lld\n", target->_[displace_index].key, key);
      prev_neighbor->Insert(target->_[displace_index].key,
                            target->_[displace_index].value,
                            target->finger_array[displace_index], false);
      //target->setPreNonFlush();
      prev_neighbor->release_lock();
      //target->unsetPreNonFlush();
      target->unset_hash(displace_index);
      target->Insert_displace(key, value, meta_hash, displace_index, false);
      //neighbor->setPreNonFlush();
      target->release_lock();
      //neighbor->unsetPreNonFlush();
      neighbor->release_lock();
#ifdef COUNTING
      __sync_fetch_and_add(&number, 1);
#endif
      return 0;
    }
    return -1;
  }

  int Stash_insert(Bucket<T> *target, Bucket<T> *neighbor, T key, Value_t value,
                   uint8_t meta_hash, int stash_pos) {
    for (int i = 0; i < stashBucket; ++i) {
      Bucket<T> *curr_bucket = bucket + kNumBucket + ((stash_pos + i) & stashMask);
      if (GET_COUNT(curr_bucket->bitmap) < kNumPairPerBucket) {
        // printf("insertion in the stash for key %lld\n", key);
        curr_bucket->Insert(key, value, meta_hash, false);
        target->set_indicator(meta_hash, neighbor, (stash_pos + i) & stashMask);
#ifdef COUNTING
        __sync_fetch_and_add(&number, 1);
#endif
        return 0;
      }
    }
    return -1;
  }

  Bucket<T> *bucket;
  size_t kNumBucket;
  size_t stashBucket;
  size_t bucketMask;
  size_t stashMask;
  uint8_t stashHighMask;
  size_t local_depth;
  int number;
};

/* it needs to verify whether this bucket has been deleted...*/
template<class T>
int Table<T>::Insert(T key, Value_t value) {
  auto key_hash = h(&key, sizeof(key));
  auto meta_hash = ((uint8_t)(key_hash & kMask)); 
  /*we need to first do the locking and then do the verify*/
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);
  // printf("for key %lld, target bucket is %lld, meta_hash is %d\n", key,
  // BUCKET_INDEX(key_hash), meta_hash);
  target->get_lock();
  if (!neighbor->try_get_lock()) {
    target->release_lock();
    return -2;
  }

  /*unique check, needs to check 2 hash table*/
  
  auto ret =
      target->unique_check(meta_hash, key, neighbor, bucket + kNumBucket, stashBucket);
  if (ret == -1) {
    neighbor->release_lock();
    target->release_lock();
    return 0;
  }

  if (((GET_COUNT(target->bitmap)) == kNumPairPerBucket) &&
      ((GET_COUNT(neighbor->bitmap)) == kNumPairPerBucket)) {
      Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
      // Next displacement
      if (!next_neighbor->try_get_lock()) {
        neighbor->release_lock();
        target->release_lock();
        return -2;
      }
      auto ret = Next_displace(target, neighbor, next_neighbor, key, value, meta_hash);
      if (ret == 0) {
        return 0;
      }
      next_neighbor->release_lock();

      Bucket<T> *prev_neighbor;
      int prev_index;
      if (y == 0) {
        prev_neighbor = bucket + kNumBucket - 1;
        prev_index = kNumBucket - 1;
      } else {
        prev_neighbor = bucket + y - 1;
        prev_index = y - 1;
      }
      if (!prev_neighbor->try_get_lock()) {
        target->release_lock();
        neighbor->release_lock();
        return -2;
      }

      ret = Prev_displace(target, prev_neighbor, neighbor,key, value, meta_hash);
      if (ret == 0) {
        return 0;
      }

      Bucket<T> *stash = bucket + kNumBucket;
      if (!stash->try_get_lock()) {
        neighbor->release_lock();
        target->release_lock();
        prev_neighbor->release_lock();
        return -2;
      }
      ret = Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);

      stash->release_lock();
      neighbor->release_lock();
      target->release_lock();
      prev_neighbor->release_lock();
      return ret;
  }
  
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    target->Insert(key, value, meta_hash, false);
    target->release_lock();
    neighbor->release_lock();
  } else {
    neighbor->Insert(key, value, meta_hash, true);
    neighbor->release_lock();
    target->release_lock();
  }
  /*
  if(GET_COUNT(target->bitmap) < 14){
    target->Insert(key, value, meta_hash, false);
    target->release_lock();
    neighbor->release_lock();
  }else{
    if(GET_COUNT(neighbor->bitmap) < 14){
      neighbor->Insert(key, value, meta_hash, true);
      target->release_lock();
      neighbor->release_lock();
      __sync_fetch_and_add(&number, 1);
      return 0;
    }
    target->release_lock();
    neighbor->release_lock();
    return -1;
  }*/
#ifdef COUNTING
  __sync_fetch_and_add(&number, 1);
#endif
  return 0;
}

int main(int argc, char const *argv[]){
    Table<Key_t>* table;
    int allocated_num = 20000000;
    uint64_t* workload = reinterpret_cast<uint64_t *>(malloc(sizeof(uint64_t) * allocated_num));
    int i;
    unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL};
    int length = 4;
    init_by_array64(init, length);

    int test_num = 500;
    double *load_factor = new double[test_num];
    double *raw_space = new double[test_num];
    uint64_t *insert_num = new uint64_t[test_num];
    uint64_t arr_kNumBucket[6];
    arr_kNumBucket[0] = 4;
    arr_kNumBucket[1] = 8;
    arr_kNumBucket[2] = 16;
    arr_kNumBucket[3] = 32;
    arr_kNumBucket[4] = 64;
    arr_kNumBucket[5] = 128;
    uint64_t arr_kStashBucket[6];
    for(int i = 0; i < 6; ++i){
      arr_kStashBucket[i] = 4;
    }
    int stop = 0;

    /* Generate Workload*/
    for(int i = 0; i < allocated_num; ++i){
        workload[i] = genrand64_int64();
    }

    for(int k = 0; k < 6; ++k){
        table = new Table<Key_t>(arr_kNumBucket[k], arr_kStashBucket[k]);

        for(int j = 0; j < test_num; ++j){
          for(int i = stop; i < allocated_num; ++i){
              Key_t _key = workload[i];
              Value_t _value = reinterpret_cast<Value_t>(workload[i]);
              int ret = table->Insert(_key, _value);
              if(ret == -1) {
                  stop = i + 1; 
                  break;
              }
          }
          load_factor[j] = table->getLoadFactor();
          raw_space[j] = table->getRawSpace();
          insert_num[j] = table->getNumber();
          //table->displayLoadFactor();
          table->reset();
          //std::cout << "stop = " << stop << std::endl;
        }
        double avg_load_factor = 0;
        double avg_raw_space = 0;
        uint64_t avg_insert_num = 0;
        
        for(int i = 0; i < test_num; ++i){
          avg_load_factor += load_factor[i];
          avg_raw_space += raw_space[i];
          avg_insert_num += insert_num[i];
        }
        std::cout << "---------------------For Bucket Number = " << arr_kNumBucket[k] << "; Stash Bucket Number = "<< arr_kStashBucket[k]<< std::endl;
        std::cout << "Avg insert num = " << avg_insert_num / test_num << std::endl;
        //std::cout << "Avg Load factor = " << avg_load_factor / test_num << std::endl;
        //std::cout << "Avg Raw space = " << avg_raw_space / test_num << std::endl;
        std::cout << "Avg Load factor = " << (double)avg_insert_num / test_num / (14 * (arr_kNumBucket[k] + arr_kStashBucket[k])) << std::endl;
        std::cout << "Avg Raw Sapce = " << (double)avg_insert_num / test_num / (16 * (arr_kNumBucket[k] + arr_kStashBucket[k])) << std::endl;

        for(int i = 0; i < test_num; ++i){
          load_factor[i] = 0;
          raw_space[i] = 0;
          insert_num[i] = 0;
        }
        delete table;
        stop = 0;
    }
    return 0;
}