// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
//
// Dash Linear Hashing
// Authors:
// Baotong Lu <btlu@cse.cuhk.edu.hk>
// Xiangpeng Hao <xiangpeng_hao@sfu.ca>
// Tianzheng Wang <tzwang@sfu.ca>

#pragma once

#include <immintrin.h>

#include <bitset>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <shared_mutex>
#include <thread>
#include <unordered_map>
#include <vector>

#include "../util/hash.h"
#include "../util/pair.h"
#include "Hash.h"
#include "allocator.h"
#define DOUBLE_EXPANSION 1

#ifdef PMEM
#include <libpmemobj.h>
#endif
namespace linear {
//#define PREALLOC 1
//#define COUNTING 1

template <class T>
struct _Pair {
  T key;
  Value_t value;
};

const size_t k_PairSize = 16;
const uint32_t lockSet = 1 << 31;
const uint32_t lockMask = ((uint32_t)1 << 31) - 1;
const int overflowSet = 1 << 4;
const int countMask = (1 << 4) - 1;
const uint32_t initialSet = 1 << 30;
const uint32_t versionMask = (1 << 30) - 1;
const size_t kNumPairPerBucket = 14;
const size_t kFingerBits = 8;
const uint32_t kNumBucket = 64;
const uint32_t stashBucket = 2;
const uint64_t recoverBit = 1UL << 63;
const uint64_t lockBit = 1UL << 62;
const uint8_t overflowBitmapMask = (1 << 4) - 1;
const uint8_t low2Mask = (1 << 2) - 1;
const uint32_t low4Mask = (1 << 4) - 1;
const uint32_t fixedExpandNum = 8; /*stride*/
constexpr uint32_t initialLockSet = lockSet | initialSet;
constexpr uint32_t fixedExpandBits = 31 - __builtin_clz(fixedExpandNum);
constexpr uint32_t fixedExpandMask = (1 << fixedExpandBits) - 1;
constexpr size_t kMask = (1 << kFingerBits) - 1;
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
constexpr size_t bucketMask = ((1 << (31 - __builtin_clz(kNumBucket))) - 1);
constexpr size_t stashMask = (1 << (31 - __builtin_clz(stashBucket))) - 1;
constexpr uint8_t stashHighMask = ~((uint8_t)stashMask);
constexpr uint32_t segmentSize = 64;
constexpr size_t baseShifBits =
    static_cast<uint64_t>(31 - __builtin_clz(segmentSize));
constexpr uint32_t segmentMask = (1 << baseShifBits) - 1;
constexpr size_t directorySize = 1024 * 4;
constexpr uint64_t low32Mask = ((uint64_t)1 << 32) - 1;
constexpr uint64_t high32Mask = ~low32Mask;
constexpr uint32_t shiftBits = (31 - __builtin_clz(kNumBucket)) + kFingerBits;
constexpr uint32_t expandShiftBits = fixedExpandBits + baseShifBits;
constexpr uint64_t recoverLockBit = recoverBit | lockBit;

#define BUCKET_INDEX(hash) (((hash) >> (64 - shiftBits)) & bucketMask)
#define META_HASH(hash) ((uint8_t)((hash) >> (64 - kFingerBits)))
#define GET_COUNT(var) ((var)&countMask)
#define GET_MEMBER(var) (((var) >> 4) & allocMask)
#define GET_INVERSE_MEMBER(var) ((~((var) >> 4)) & allocMask)
#define GET_BITMAP(var) ((var) >> 18)

inline bool var_compare(char *str1, char *str2, int len1, int len2) {
  if (len1 != len2) return false;
  return !memcmp(str1, str2, len1);
}

inline uint32_t pow2(int shift_index) { return (1 << shift_index); }

inline uint64_t IDX(uint64_t hashKey, uint64_t level) {
  return hashKey & ((1 << level) - 1);
}

inline void SEG_IDX_OFFSET(uint32_t bucket_idx, uint32_t &seg_idx,
                           uint32_t &seg_offset) {
#ifdef DOUBLE_EXPANSION
  int index = 31 - __builtin_clz((bucket_idx >> expandShiftBits) + 1);
  int shift =
      ((bucket_idx >> baseShifBits) - (fixedExpandNum * (pow2(index) - 1))) >>
      index;
  seg_idx = fixedExpandNum * index + shift;
  seg_offset = bucket_idx - segmentSize * ((pow2(index) - 1) * fixedExpandNum +
                                           pow2(index) * shift);
#else
  seg_idx = bucket_idx >> baseShifBits;
  seg_offset = bucket_idx & segmentMask;
#endif
}

uint64_t SUM_BUCKET(uint32_t bucket_idx) {
  int index = 31 - __builtin_clz((bucket_idx >> expandShiftBits) + 1);
  int shift =
      (((bucket_idx >> baseShifBits) - (fixedExpandNum * (pow2(index) - 1))) >>
       index);
  uint64_t sum = segmentSize * ((pow2(index) - 1) * fixedExpandNum +
                                pow2(index) * (shift + 1));
  return sum;
}

/* return the size of the segment array (#segments) that the corresponding
 * bucket resides in*/
inline uint32_t SEG_SIZE(uint32_t bucket_ix) {
  int index = 31 - __builtin_clz((bucket_ix >> expandShiftBits) + 1);
  return segmentSize * pow2(index);
}

inline uint32_t SEG_SIZE_BY_SEGARR_ID(uint32_t segarr_idx) {
  return segmentSize * pow2(segarr_idx / fixedExpandNum);
}

/* overflow Bucket*/
template <class T>
struct overflowBucket {
  overflowBucket() { memset(this, 0, sizeof(struct overflowBucket)); }

  inline int find_empty_slot() {
    if (GET_COUNT(bitmap) == kNumPairPerBucket) {
      return -1;
    }
    auto mask = ~(GET_BITMAP(bitmap));
    return __builtin_ctz(mask);
  }

  Value_t check_and_get(uint8_t meta_hash, T key) {
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    mask = mask & GET_BITMAP(bitmap);

    if constexpr (std::is_pointer_v<T>) {
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) &&
              (var_compare(_[i].key->key, key->key, _[i].key->length,
                           key->length))) {
            return _[i].value;
          }

          if (CHECK_BIT(mask, i + 1) &&
              (var_compare(_[i + 1].key->key, key->key, _[i + 1].key->length,
                           key->length))) {
            return _[i + 1].value;
          }

          if (CHECK_BIT(mask, i + 2) &&
              (var_compare(_[i + 2].key->key, key->key, _[i + 2].key->length,
                           key->length))) {
            return _[i + 2].value;
          }

          if (CHECK_BIT(mask, i + 3) &&
              (var_compare(_[i + 3].key->key, key->key, _[i + 3].key->length,
                           key->length))) {
            return _[i + 3].value;
          }
        }

        if (CHECK_BIT(mask, 12) &&
            (var_compare(_[12].key->key, key->key, _[12].key->length,
                         key->length))) {
          return _[12].value;
        }

        if (CHECK_BIT(mask, 13) &&
            (var_compare(_[13].key->key, key->key, _[13].key->length,
                         key->length))) {
          return _[13].value;
        }
      }
    } else {
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

  inline void set_hash(int index, uint8_t meta_hash) {
    finger_array[index] = meta_hash;
    uint32_t new_bitmap = bitmap | (1 << (index + 18));
    new_bitmap++;
    bitmap = new_bitmap;
  }

  inline void unset_hash(int index) {
    uint32_t new_bitmap = bitmap & (~(1 << (index + 18)));
    new_bitmap--;
    bitmap = new_bitmap;
  }

  inline void get_lock() {
    uint32_t new_value = 0;
    uint32_t old_value = 0;
    do {
      while (true) {
        old_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        if (!(old_value & lockSet)) {
          old_value &= lockMask;
          break;
        }
      }
      new_value = old_value | lockSet;
    } while (!CAS(&version_lock, &old_value, new_value));
  }

  inline void release_lock() {
    uint32_t v = version_lock;
    __atomic_store_n(&version_lock, v + 1 - lockSet, __ATOMIC_RELEASE);
  }

  int Insert(T key, Value_t value, uint8_t meta_hash) {
    auto slot = find_empty_slot();
    /* this branch can be optimized out*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
#ifdef PMEM
    Allocator::Persist(&_[slot], sizeof(_[slot]));
#endif
    set_hash(slot, meta_hash);
    return 0;
  }

  int Delete(uint8_t meta_hash, T key) {
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    mask = mask & GET_BITMAP(bitmap);
    /*loop unrolling*/
    if constexpr (std::is_pointer_v<T>) {
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) &&
              (var_compare(_[i].key->key, key->key, _[i].key->length,
                           key->length))) {
            unset_hash(i);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) &&
              (var_compare(_[i + 1].key->key, key->key, _[i + 1].key->length,
                           key->length))) {
            unset_hash(i + 1);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) &&
              (var_compare(_[i + 2].key->key, key->key, _[i + 2].key->length,
                           key->length))) {
            unset_hash(i + 2);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) &&
              (var_compare(_[i + 3].key->key, key->key, _[i + 3].key->length,
                           key->length))) {
            unset_hash(i + 3);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) &&
            (var_compare(_[12].key->key, key->key, _[12].key->length,
                         key->length))) {
          unset_hash(12);
          return 0;
        }

        if (CHECK_BIT(mask, 13) &&
            (var_compare(_[13].key->key, key->key, _[13].key->length,
                         key->length))) {
          unset_hash(13);
          return 0;
        }
      }
    } else {
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (_[i].key == key)) {
            unset_hash(i);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) && (_[i + 1].key == key)) {
            unset_hash(i + 1);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) && (_[i + 2].key == key)) {
            unset_hash(i + 2);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) && (_[i + 3].key == key)) {
            unset_hash(i + 3);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) && (_[12].key == key)) {
          unset_hash(12);
          return 0;
        }

        if (CHECK_BIT(mask, 13) && (_[13].key == key)) {
          unset_hash(13);
          return 0;
        }
      }
    }
    return -1;
  }

  int Insert_with_noflush(T key, Value_t value, uint8_t meta_hash) {
    auto slot = find_empty_slot();
    /* this branch can be optimized out*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash);
    return 0;
  }

  inline void resetLock() { version_lock = version_lock & initialSet; }

  uint32_t version_lock;
  uint32_t bitmap;
  uint8_t finger_array[kNumPairPerBucket + 2];
  overflowBucket *next;
  _Pair<T> _[kNumPairPerBucket];
};

/* normal bucket*/
template <class T>
struct Bucket {
  inline int find_empty_slot() {
    if (GET_COUNT(bitmap) == kNumPairPerBucket) {
      return -1;
    }
    auto mask = ~(GET_BITMAP(bitmap));
    return __builtin_ctz(mask);
  }

  /*true indicates overflow, needs extra check in the stash*/
  inline bool test_overflow() { return overflowCount; }

  inline bool test_stash_check() { return (overflowBitmap & overflowSet); }

  inline void clear_stash_check() {
    overflowBitmap = overflowBitmap & (~overflowSet);
  }

  inline void set_indicator(uint8_t meta_hash, Bucket<T> *neighbor,
                            uint8_t pos) {
    int mask = overflowBitmap & overflowBitmapMask;
    mask = ~mask;
    auto index = __builtin_ctz(mask);

    if (index < 4) {
      finger_array[14 + index] = meta_hash;
      overflowBitmap =
          ((uint8_t)(1 << index) | overflowBitmap); /*may be optimized*/
      overflowIndex =
          (overflowIndex & (~(3 << (index * 2)))) | (pos << (index * 2));
    } else {
      mask = neighbor->overflowBitmap & overflowBitmapMask;
      mask = ~mask;
      index = __builtin_ctz(mask);
      if (index < 4) {
        neighbor->finger_array[14 + index] = meta_hash;
        neighbor->overflowBitmap =
            ((uint8_t)(1 << index) | neighbor->overflowBitmap);
        neighbor->overflowMember =
            ((uint8_t)(1 << index) | neighbor->overflowMember);
        neighbor->overflowIndex =
            (neighbor->overflowIndex & (~(3 << (index * 2)))) |
            (pos << (index * 2));
      } else { /*overflow, increase count*/
        assert(overflowCount < 255);
        overflowCount++;
      }
    }
    overflowBitmap = overflowBitmap | overflowSet;
  }

  /*both clear this bucket and its neighbor bucket*/
  inline bool unset_indicator(uint8_t meta_hash, Bucket<T> *neighbor, T key,
                              uint64_t pos) {
    /*also needs to ensure that this meta_hash must belongs to other bucket*/
    bool clear_success = false;
    int mask1 = overflowBitmap & overflowBitmapMask;
    for (int i = 0; i < 4; ++i) {
      if (CHECK_BIT(mask1, i) && (finger_array[14 + i] == meta_hash) &&
          (((1 << i) & overflowMember) == 0) &&
          (((overflowIndex >> (2 * i)) & low2Mask) == pos)) {
        overflowBitmap = overflowBitmap & ((uint8_t)(~(1 << i)));
        overflowIndex = overflowIndex & (~(3 << (i * 2)));
        assert(((overflowIndex >> (i * 2)) & stashMask) == 0);
        clear_success = true;
        break;
      }
    }

    int mask2 = neighbor->overflowBitmap & overflowBitmapMask;
    if (!clear_success) {
      for (int i = 0; i < 4; ++i) {
        if (CHECK_BIT(mask2, i) &&
            (neighbor->finger_array[14 + i] == meta_hash) &&
            (((1 << i) & neighbor->overflowMember) != 0) &&
            (((neighbor->overflowIndex >> (2 * i)) & low2Mask) == pos)) {
          neighbor->overflowBitmap =
              neighbor->overflowBitmap & ((uint8_t)(~(1 << i)));
          neighbor->overflowMember =
              neighbor->overflowMember & ((uint8_t)(~(1 << i)));
          neighbor->overflowIndex = neighbor->overflowIndex & (~(3 << (i * 2)));
          assert(((neighbor->overflowIndex >> (i * 2)) & stashMask) == 0);
          clear_success = true;
          break;
        }
      }
    }

    if (!clear_success) {
      assert(overflowCount != 0);
      overflowCount--;
    }

    mask1 = overflowBitmap & overflowBitmapMask;
    mask2 = neighbor->overflowBitmap & overflowBitmapMask;
    if (((mask1 & (~overflowMember)) == 0) && (overflowCount == 0) &&
        ((mask2 & neighbor->overflowMember) == 0)) {
      clear_stash_check();
    }
    return true;
  }

  int unique_check(uint8_t meta_hash, T key, Bucket<T> *neighbor,
                   overflowBucket<T> *stash) {
    if ((check_and_get(meta_hash, key, false) != NONE) ||
        (neighbor->check_and_get(meta_hash, key, true) != NONE)) {
      return -1;
    }

    if (test_stash_check()) {
      auto test_stash = false;
      if (test_overflow()) {
        test_stash = true;
      } else {
        int mask = overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) && (finger_array[14 + i] == meta_hash) &&
                (((1 << i) & overflowMember) == 0)) {
              test_stash = true;
              goto STASH_CHECK;
            }
          }
        }

        mask = neighbor->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (neighbor->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & neighbor->overflowMember) != 0)) {
              test_stash = true;
              break;
            }
          }
        }
      }
    STASH_CHECK:
      if (test_stash) {
        for (int i = 0; i < stashBucket; ++i) {
          overflowBucket<T> *curr_bucket = stash + i;
          auto ret = curr_bucket->check_and_get(meta_hash, key);
          if (ret != NONE) {
            return -1;
          }
        }

        overflowBucket<T> *prev_bucket = stash;
        overflowBucket<T> *next_bucket = stash->next;
        while (next_bucket != NULL) {
          auto ret = next_bucket->check_and_get(meta_hash, key);
          if (ret != NONE) {
            return -1;
          }
          prev_bucket = next_bucket;
          next_bucket = next_bucket->next;
        }
      }
    }
    return 0;
  }

  inline int get_current_mask() {
    int mask = GET_BITMAP(bitmap) & GET_INVERSE_MEMBER(bitmap);
    return mask;
  }

  Value_t check_and_get(uint8_t meta_hash, T key, bool probe) {
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & GET_INVERSE_MEMBER(bitmap);
    } else {
      mask = mask & GET_BITMAP(bitmap) & GET_MEMBER(bitmap);
    }

    if constexpr (std::is_pointer_v<T>) {
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) &&
              (var_compare(_[i].key->key, key->key, _[i].key->length,
                           key->length))) {
            return _[i].value;
          }

          if (CHECK_BIT(mask, i + 1) &&
              (var_compare(_[i + 1].key->key, key->key, _[i + 1].key->length,
                           key->length))) {
            return _[i + 1].value;
          }

          if (CHECK_BIT(mask, i + 2) &&
              (var_compare(_[i + 2].key->key, key->key, _[i + 2].key->length,
                           key->length))) {
            return _[i + 2].value;
          }

          if (CHECK_BIT(mask, i + 3) &&
              (var_compare(_[i + 3].key->key, key->key, _[i + 3].key->length,
                           key->length))) {
            return _[i + 3].value;
          }
        }

        if (CHECK_BIT(mask, 12) &&
            (var_compare(_[12].key->key, key->key, _[12].key->length,
                         key->length))) {
          return _[12].value;
        }

        if (CHECK_BIT(mask, 13) &&
            (var_compare(_[13].key->key, key->key, _[13].key->length,
                         key->length))) {
          return _[13].value;
        }
      }
    } else {
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

  inline void set_hash(int index, uint8_t meta_hash, bool probe) {
    finger_array[index] = meta_hash;
    uint32_t new_bitmap = bitmap | (1 << (index + 18));
    if (probe) {
      new_bitmap = new_bitmap | (1 << (index + 4));
    }
    assert(GET_COUNT(bitmap) < kNumPairPerBucket);
    new_bitmap++;
    bitmap = new_bitmap;
  }

  inline uint8_t get_hash(int index) { return finger_array[index]; }

  inline void unset_hash(int index) {
    uint32_t new_bitmap =
        bitmap & (~(1 << (index + 18))) & (~(1 << (index + 4)));
    assert(GET_COUNT(bitmap) <= kNumPairPerBucket);
    assert(GET_COUNT(bitmap) > 0);
    new_bitmap--;
    bitmap = new_bitmap;
  }

  inline void get_lock() {
    uint32_t new_value = 0;
    uint32_t old_value = 0;
    do {
      while (true) {
        old_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        if (!(old_value & lockSet)) {
          old_value &= lockMask;
          break;
        }
      }
      new_value = old_value | lockSet;
    } while (!CAS(&version_lock, &old_value, new_value));
  }

  inline bool try_get_lock() {
    uint32_t v = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    if (v & lockSet) {
      return false;
    }
    auto old_value = v & lockMask;
    auto new_value = v | lockSet;
    return CAS(&version_lock, &old_value, new_value);
  }

  inline void release_lock() {
    uint32_t v = version_lock;
    __atomic_store_n(&version_lock, v + 1 - lockSet, __ATOMIC_RELEASE);
  }

  /*if the lock is set, return true*/
  inline bool test_lock_set(uint32_t &version) {
    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    version = value;
    return (value & lockSet);
  }

  // test whether the version has change, if change, return true
  inline bool test_lock_version_change(uint32_t old_version) {
    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    return (value != old_version);
  }

  /*if it has been initialized, then return true*/
  inline bool test_initialize() { return (version_lock & initialSet); }

  inline void set_initialize() { version_lock = version_lock | initialSet; }

  inline void unset_initialize() {
    version_lock = version_lock & (!initialSet);
  }

  int Insert(T key, Value_t value, uint8_t meta_hash, bool probe) {
    auto slot = find_empty_slot();
    /* this branch can be optimized out*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
#ifdef PMEM
    Allocator::Persist(&_[slot], sizeof(_[slot]));
#endif
    set_hash(slot, meta_hash, probe);
    return 0;
  }

  /*if delete success, then return 0, else return -1*/
  int Delete(uint8_t meta_hash, T key, bool probe) {
    /*do the simd and check the key, then do the delete operation*/
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & GET_INVERSE_MEMBER(bitmap);
    } else {
      mask = mask & GET_BITMAP(bitmap) & GET_MEMBER(bitmap);
    }
    /*loop unrolling*/
    if constexpr (std::is_pointer_v<T>) {
      string_key *_key = reinterpret_cast<string_key *>(key);
      /*loop unrolling*/
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) &&
              (var_compare(_[i].key->key, key->key, _[i].key->length,
                           key->length))) {
            unset_hash(i);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) &&
              (var_compare(_[i + 1].key->key, key->key, _[i + 1].key->length,
                           key->length))) {
            unset_hash(i + 1);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) &&
              (var_compare(_[i + 2].key->key, key->key, _[i + 2].key->length,
                           key->length))) {
            unset_hash(i + 2);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) &&
              (var_compare(_[i + 3].key->key, key->key, _[i + 3].key->length,
                           key->length))) {
            unset_hash(i + 3);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) &&
            (var_compare(_[12].key->key, key->key, _[12].key->length,
                         key->length))) {
          unset_hash(12);
          return 0;
        }

        if (CHECK_BIT(mask, 13) &&
            (var_compare(_[13].key->key, key->key, _[13].key->length,
                         key->length))) {
          unset_hash(13);
          return 0;
        }
      }
    } else {
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (_[i].key == key)) {
            unset_hash(i);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) && (_[i + 1].key == key)) {
            unset_hash(i + 1);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) && (_[i + 2].key == key)) {
            unset_hash(i + 2);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) && (_[i + 3].key == key)) {
            unset_hash(i + 3);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) && (_[12].key == key)) {
          unset_hash(12);
          return 0;
        }

        if (CHECK_BIT(mask, 13) && (_[13].key == key)) {
          unset_hash(13);
          return 0;
        }
      }
    }
    return -1;
  }

  int Insert_with_noflush(T key, Value_t value, uint8_t meta_hash, bool probe) {
    auto slot = find_empty_slot();

    /* this branch can be optimized out*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
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
#ifdef PMEM
    Allocator::Persist(&_[slot], sizeof(_Pair<T>));
#endif
    set_hash(slot, meta_hash, probe);
  }

  void Insert_displace_with_noflush(T key, Value_t value, uint8_t meta_hash,
                                    int slot, bool probe) {
    assert(key != 0);
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
  }

  /* Find the displacment element in this bucket*/
  inline int Find_org_displacement() {
    uint32_t mask = GET_INVERSE_MEMBER(bitmap);
    if (mask == 0) {
      return -1;
    }
    return __builtin_ctz(mask);
  }

  /*find element that it is in the probe*/
  inline int Find_probe_displacement() {
    uint32_t mask = GET_MEMBER(bitmap);
    if (mask == 0) {
      return -1;
    }
    return __builtin_ctz(mask);
  }

  inline void resetLock() { version_lock = version_lock & initialSet; }

  inline void resetOverflowFP() {
    overflowBitmap = 0;
    overflowIndex = 0;
    overflowMember = 0;
    overflowCount = 0;
    clear_stash_check();
  }

  uint32_t version_lock;
  uint32_t bitmap;          // allocation bitmap + pointer bitmao + counter
  uint8_t finger_array[18]; /*only use the first 14 bytes, can be accelerated by
                               SSE instruction,0-13 for finger, 14-17 for
                               overflowed*/
  uint8_t overflowBitmap;
  uint8_t overflowIndex;
  uint8_t overflowMember; /*overflowmember indicates membership of the overflow
                             fingerprint*/
  uint8_t overflowCount;
  uint8_t unused[2];
  _Pair<T> _[kNumPairPerBucket];
};

template <class T>
struct Table;

template <class T>
struct Directory {
  typedef Table<T> *table_p;
  uint64_t N_next;
  uint64_t recovered_index; /* Used to indicate the last segment that needs to
                               be recovered*/
  table_p _[directorySize];
  uint64_t recover_counter[directorySize];
  uint64_t crash_version; /* it does not influence the correctness*/

  static void New(PMEMoid *dir) {
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
      auto dir_ptr = reinterpret_cast<Directory<T> *>(ptr);
      dir_ptr->N_next = baseShifBits << 32;
      dir_ptr->recovered_index = 0;
      dir_ptr->crash_version = 0;
      memset(&dir_ptr->_, 0, sizeof(table_p) * directorySize);
      memset(&dir_ptr->recover_counter, 0, sizeof(uint64_t) * directorySize);
      return 0;
    };

    Allocator::Allocate(dir, kCacheLineSize, sizeof(Directory<T>), callback,
                        NULL);
  }
};

/*thread local table allcoation pool*/
template <class T>
struct TlsTablePool {
  static Table<T> *all_tables;
  static PMEMoid p_all_tables;
  static std::atomic<uint32_t> all_allocated;
  static const uint32_t kAllTables = 327680;

  static void AllocateMore() {
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg) { return 0; };
    std::pair callback_para(0, nullptr);
    Allocator::Allocate(&p_all_tables, kCacheLineSize,
                        sizeof(Table<T>) * kAllTables, callback,
                        reinterpret_cast<void *>(&callback_para));
    all_tables = reinterpret_cast<Table<T> *>(pmemobj_direct(p_all_tables));
    memset((void *)all_tables, 0, sizeof(Table<T>) * kAllTables);
    all_allocated = 0;
    printf("MORE ");
  }

  TlsTablePool() {}
  static void Initialize() { AllocateMore(); }

  Table<T> *tables = nullptr;
  static const uint32_t kTables = 128;
  uint32_t allocated = kTables;

  void TlsPrepare() {
  retry:
    uint32_t n = all_allocated.fetch_add(kTables);
    if (n == kAllTables) {
      AllocateMore();
      printf("NO MORE\n");
      abort();
      goto retry;
    }
    tables = all_tables + n;
    allocated = 0;
  }

  Table<T> *Get() {
    if (allocated == kTables) {
      TlsPrepare();
    }
    return &tables[allocated++];
  }

  /*allocate a segment from preallocated memory*/
  static Table<T> *Get(uint64_t seg_size) {
    uint32_t n = all_allocated.fetch_add(seg_size);
    if (n >= kAllTables) {
      AllocateMore();
      abort();
    }
    return &all_tables[n];
  }
};

template <class T>
std::atomic<uint32_t> TlsTablePool<T>::all_allocated(0);
template <class T>
Table<T> *TlsTablePool<T>::all_tables = nullptr;
template <class T>
PMEMoid TlsTablePool<T>::p_all_tables = OID_NULL;

/* the meta hash-table referenced by the directory*/
template <class T>
struct Table {
  Table(void) {
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = bucket + i;
      memset(curr_bucket, 0, 64);
    }

    for (int i = 0; i < stashBucket; ++i) {
      overflowBucket<T> *curr_bucket = stash + i;
      memset(curr_bucket, 0, 64);
    }
  }

  static void New(Table<T> **tbl) {
    Allocator::ZAllocate((void **)tbl, kCacheLineSize, sizeof(Table<T>));
  };

  ~Table(void) {}

  int Insert(T key, Value_t value, size_t key_hash, Directory<T> *_dir,
             uint64_t index, uint32_t old_N, uint32_t old_next);
  void Insert4split(T key, Value_t value, size_t key_hash, uint8_t meta_hash);
  void Insert4merge(T key, Value_t value, size_t key_hash, uint8_t meta_hash,
                    bool flag = false);
  void Merge(Table<T> *neighbor, bool flag = false);
  void Split(Table<T> *org_table, uint64_t base_level, int org_idx,
             Directory<T> *);
  int Insert2Org(T key, Value_t value, size_t key_hash, size_t pos);
  void PrintTableImage(Table<T> *table, uint64_t base_level);

  void getAllLocks() {
    Bucket<T> *curr_bucket;
    for (int i = 0; i < kNumBucket; ++i) {
      curr_bucket = bucket + i;
      curr_bucket->get_lock();
    }
  }

  void releaseAllLocks() {
    Bucket<T> *curr_bucket;
    for (int i = 0; i < kNumBucket; ++i) {
      curr_bucket = bucket + i;
      curr_bucket->release_lock();
    }
  }

  int Next_displace(Bucket<T> *neighbor, Bucket<T> *next_neighbor, T key,
                    Value_t value, uint8_t meta_hash) {
    int displace_index = neighbor->Find_org_displacement();
    if ((GET_COUNT(next_neighbor->bitmap) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert(neighbor->_[displace_index].key,
                            neighbor->_[displace_index].value,
                            neighbor->finger_array[displace_index], true);
      next_neighbor->release_lock();
#ifdef PMEM
      Allocator::Persist(&next_neighbor->bitmap, sizeof(next_neighbor->bitmap));
#endif
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace(key, value, meta_hash, displace_index, true);
      neighbor->release_lock();
#ifdef PMEM
      Allocator::Persist(&neighbor->bitmap, sizeof(neighbor->bitmap));
#endif
#ifdef COUNTING
      __sync_fetch_and_add(&number, 1);
#endif
      return 0;
    }
    return -1;
  }

  int Prev_displace(Bucket<T> *target, Bucket<T> *prev_neighbor, T key,
                    Value_t value, uint8_t meta_hash) {
    int displace_index = target->Find_probe_displacement();
    if ((GET_COUNT(prev_neighbor->bitmap) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert(target->_[displace_index].key,
                            target->_[displace_index].value,
                            target->finger_array[displace_index], false);
      prev_neighbor->release_lock();
#ifdef PMEM
      Allocator::Persist(&prev_neighbor->bitmap, sizeof(prev_neighbor->bitmap));
#endif
      target->unset_hash(displace_index);
      target->Insert_displace(key, value, meta_hash, displace_index, false);
      target->release_lock();
#ifdef PMEM
      Allocator::Persist(&target->bitmap, sizeof(target->bitmap));
#endif
#ifdef COUNTING
      __sync_fetch_and_add(&number, 1);
#endif
      return 0;
    }
    return -1;
  }

  /*insertion in the corresponding position of the stash*/
  int Stash_insert(Bucket<T> *target, Bucket<T> *neighbor, T key, Value_t value,
                   uint8_t meta_hash, int stash_pos) {
    for (int i = 0; i < stashBucket; ++i) {
      overflowBucket<T> *curr_bucket = stash + ((stash_pos + i) & stashMask);
      if (GET_COUNT(curr_bucket->bitmap) < kNumPairPerBucket) {
        curr_bucket->Insert(key, value, meta_hash);
#ifdef PMEM
        Allocator::Persist(&curr_bucket->bitmap, sizeof(curr_bucket->bitmap));
#endif
        target->set_indicator(meta_hash, neighbor, (stash_pos + i) & stashMask);
#ifdef COUNTING
        __sync_fetch_and_add(&number, 1);
#endif
        return 0;
      }
    }

    /*need to add the handling the overflowed chaining*/
    overflowBucket<T> *prev_bucket = stash;
    overflowBucket<T> *next_bucket = stash->next;
    while (next_bucket != NULL) {
      if (GET_COUNT(next_bucket->bitmap) < kNumPairPerBucket) {
        next_bucket->Insert(key, value, meta_hash);
#ifdef PMEM
        Allocator::Persist(&next_bucket->bitmap, sizeof(next_bucket->bitmap));
#endif
        target->set_indicator(meta_hash, neighbor, 3);
#ifdef COUNTING
        __sync_fetch_and_add(&number, 1);
#endif
        return 0;
      }
      prev_bucket = next_bucket;
      next_bucket = next_bucket->next;
    }

    /*allocate new next bucket*/
    Allocator::ZAllocate((void **)&prev_bucket->next, kCacheLineSize,
                         sizeof(overflowBucket<T>));
#ifdef PMEM
    Allocator::Persist(&prev_bucket->next, sizeof(prev_bucket->next));
#endif
    prev_bucket->next->Insert(key, value, meta_hash);
#ifdef PMEM
    Allocator::Persist(&prev_bucket->next->bitmap,
                       sizeof(prev_bucket->next->bitmap));
#endif
    target->set_indicator(meta_hash, neighbor, 3);
#ifdef COUNTING
    __sync_fetch_and_add(&number, 1);
#endif
    return -1;
  }

  inline int verify_access(Directory<T> *new_dir, uint32_t index,
                           uint32_t old_N, uint32_t old_next) {
    uint64_t new_N_next = new_dir->N_next;
    uint32_t N = new_N_next >> 32;
    uint32_t next = (uint32_t)new_N_next;

    if (((old_next <= index) && (next > index)) || (old_N != N)) {
      return -1;
    }
    return 0;
  }

  /* Get its corresponding buddy table in the right direction*/
  inline Table<T> *get_expan_table(uint64_t x, uint64_t *idx,
                                   uint64_t *base_diff, Directory<T> *dir) {
    uint64_t base_level = static_cast<uint64_t>(log2(x)) + 1;
    uint64_t diff = pow2(base_level);
    *base_diff = diff;
    auto expan_idx = x + diff;
    *idx = expan_idx;
    uint32_t dir_idx;
    uint32_t offset;
    SEG_IDX_OFFSET(static_cast<uint32_t>(expan_idx), dir_idx, offset);
    if (dir->_[dir_idx] == NULL)
      return NULL;
    else
      return (dir->_[dir_idx] + offset);
  }
  /*
   *@param idx the index of the original bucket
   */
  inline Table<T> *get_org_table(uint64_t x, uint64_t *idx, uint64_t *base_diff,
                                 Directory<T> *dir) {
    uint64_t base_level = static_cast<uint64_t>(log2(x));
    uint64_t diff = static_cast<uint64_t>(pow(2, base_level));
    *base_diff = diff;
    auto org_idx = x - diff;
    *idx = org_idx;
    uint32_t dir_idx;
    uint32_t offset;
    SEG_IDX_OFFSET(static_cast<uint32_t>(org_idx), dir_idx, offset);
    Table<T> *org_table = reinterpret_cast<Table<T> *>(
                              (uint64_t)dir->_[dir_idx] & (~recoverLockBit)) +
                          offset;
    return org_table;
  }

  void recoverMetadata() {
    Bucket<T> *curr_bucket, *neighbor_bucket;
    uint64_t knumber = 0;

    for (int i = 0; i < kNumBucket; ++i) {
      curr_bucket = bucket + i;
      curr_bucket->resetLock();
      curr_bucket->resetOverflowFP();
      neighbor_bucket = bucket + ((i + 1) & bucketMask);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        int mask = curr_bucket->get_current_mask();
        if (CHECK_BIT(mask, j) && (neighbor_bucket->check_and_get(
                                       curr_bucket->finger_array[j],
                                       curr_bucket->_[j].key, true) != NONE)) {
          curr_bucket->unset_hash(j);
        }
      }

#ifdef COUNTING
      knumber += __builtin_popcount(GET_BITMAP(curr_bucket->bitmap));
#endif
    }

    for (int i = 0; i < stashBucket; ++i) {
      auto curr_bucket = stash + i;
#ifdef COUNTING
      knumber += __builtin_popcount(GET_BITMAP(curr_bucket->bitmap));
#endif
      uint64_t key_hash;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto meta_hash = META_HASH(key_hash);
          auto org_bucket = bucket + bucket_ix;
          auto neighbor_bucket = bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->set_indicator(meta_hash, neighbor_bucket, i);
        }
      }
    }

    overflowBucket<T> *prev_bucket = stash;
    overflowBucket<T> *next_bucket = stash->next;
    while (next_bucket != NULL) {
#ifdef COUNTING
      knumber += __builtin_popcount(GET_BITMAP(next_bucket->bitmap));
#endif
      uint64_t key_hash;
      auto mask = GET_BITMAP(next_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = next_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(next_bucket->_[j].key), sizeof(Key_t));
          }
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto meta_hash = META_HASH(key_hash);
          auto org_bucket = bucket + bucket_ix;
          auto neighbor_bucket = bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->set_indicator(meta_hash, neighbor_bucket, 3);
        }
      }
      prev_bucket = next_bucket;
      next_bucket = next_bucket->next;
    }

#ifdef COUNTING
    number = knumber;
#endif
  }

  Bucket<T> bucket[kNumBucket];
  overflowBucket<T> stash[stashBucket];
  int number;
  int state; /*0: normal state; 1: split bucket; 2: expand bucket(in the
                right)*/
  uint64_t seg_version;
  char dummy[48];
  PMEMmutex lock_bit;
};

template <class T>
int Table<T>::Insert2Org(T key, Value_t value, size_t key_hash, size_t pos) {
  Bucket<T> *target_bucket = bucket + pos;
  Bucket<T> *neighbor_bucket = bucket + ((pos + 1) & bucketMask);
  uint8_t meta_hash = META_HASH(key_hash);

  int target_num = GET_COUNT(target_bucket->bitmap);
  int neighbor_num = GET_COUNT(neighbor_bucket->bitmap);

  if ((target_num == kNumPairPerBucket) &&
      (neighbor_num == kNumPairPerBucket)) {
    for (int i = 0; i < stashBucket; ++i) {
      overflowBucket<T> *curr_bucket =
          stash + ((i + (pos & stashMask)) & stashMask);
      if (GET_COUNT(curr_bucket->bitmap) < kNumPairPerBucket) {
        curr_bucket->Insert(key, value, meta_hash);
#ifdef PMEM
        Allocator::Persist(&curr_bucket->bitmap, sizeof(curr_bucket->bitmap));
#endif
        target_bucket->set_indicator(meta_hash, neighbor_bucket,
                                     ((i + (pos & stashMask)) & stashMask));
        return 0;
      }
    }
    return -1;
  }

  if (target_num <= neighbor_num) {
    target_bucket->Insert(key, value, meta_hash, false);
#ifdef PMEM
    Allocator::Persist(&target_bucket->bitmap, sizeof(target_bucket->bitmap));
#endif
  } else {
    neighbor_bucket->Insert(key, value, meta_hash, true);
#ifdef PMEM
    Allocator::Persist(&neighbor_bucket->bitmap,
                       sizeof(neighbor_bucket->bitmap));
#endif
  }

  return 0;
}

/*the base_level is used to judge the rehashed key_value should be rehashed to
 * which bucket, the org_idx the index of the original table in the hash index*/
template <class T>
void Table<T>::Split(Table<T> *org_table, uint64_t base_level, int org_idx,
                     Directory<T> *_dir) {
  Bucket<T> *curr_bucket;
  for (int i = 0; i < kNumBucket; ++i) {
    curr_bucket = org_table->bucket + i;
    curr_bucket->get_lock();
  }

  if (!org_table->bucket->test_initialize()) {
    printf("recursive initiliazation\n");
    uint64_t new_org_idx;
    uint64_t new_base_level;
    Table<T> *new_org_table =
        get_org_table(org_idx, &new_org_idx, &new_base_level, _dir);
    org_table->Split(new_org_table, new_base_level, new_org_idx, _dir);
  }

  state = 2;
  Allocator::Persist(&state, sizeof(state));
  org_table->state = 1;
  Allocator::Persist(&org_table->state, sizeof(state));

  uint32_t invalid_array[kNumBucket + stashBucket];
  /*rehashing of kv objects*/
  size_t key_hash;
  for (int i = 0; i < kNumBucket; ++i) {
    curr_bucket = org_table->bucket + i;
    auto mask = GET_BITMAP(curr_bucket->bitmap);
    uint32_t invalid_mask = 0;
    for (int j = 0; j < kNumPairPerBucket; ++j) {
      if (CHECK_BIT(mask, j)) {
        if constexpr (std::is_pointer_v<T>) {
          key_hash =
              h(curr_bucket->_[j].key->key, curr_bucket->_[j].key->length);
        } else {
          key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
        }
        auto x = key_hash % (2 * base_level);
        if (x >= base_level) {
          invalid_mask = invalid_mask | (1 << j);
          Insert4split(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]);
#ifdef COUNTING
          org_table->number--;
#endif
        }
      }
    }
    invalid_array[i] = invalid_mask;
  }

  for (int i = 0; i < stashBucket; ++i) {
    overflowBucket<T> *curr_bucket = org_table->stash + i;
    auto mask = GET_BITMAP(curr_bucket->bitmap);
    uint32_t invalid_mask = 0;
    for (int j = 0; j < kNumPairPerBucket; ++j) {
      if (CHECK_BIT(mask, j)) {
        if constexpr (std::is_pointer_v<T>) {
          key_hash =
              h(curr_bucket->_[j].key->key, curr_bucket->_[j].key->length);
        } else {
          key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
        }

        auto x = key_hash % (2 * base_level);
        if (x >= base_level) {
          invalid_mask = invalid_mask | (1 << j);
          Insert4split(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]);
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto org_bucket = org_table->bucket + bucket_ix;
          auto neighbor_bucket =
              org_table->bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->unset_indicator(curr_bucket->finger_array[j],
                                      neighbor_bucket, curr_bucket->_[j].key,
                                      i);
#ifdef COUNTING
          org_table->number--;
#endif
        }
      }
    }
    invalid_array[kNumBucket + i] = invalid_mask;
  }

  seg_version = org_table->seg_version;

  /*clear the uintialized bit in expand_table*/
  for (int i = 0; i < kNumBucket; ++i) {
    curr_bucket = bucket + i;
    curr_bucket->set_initialize();
  }

  /*flush the overflowed bucketin expand table*/
  overflowBucket<T> *prev_bucket = stash;
  overflowBucket<T> *next_bucket = prev_bucket->next;

  while (next_bucket != NULL) {
#ifdef PMEM
    Allocator::Persist(next_bucket, sizeof(struct overflowBucket<T>));
#endif
    prev_bucket = next_bucket;
    next_bucket = next_bucket->next;
  }

#ifdef PMEM
  Allocator::Persist(this, sizeof(Table));
#endif

  /*clear the bitmap in original table*/
  for (int i = 0; i < kNumBucket; ++i) {
    auto curr_bucket = org_table->bucket + i;
    curr_bucket->bitmap = curr_bucket->bitmap & (~(invalid_array[i] << 18)) &
                          (~(invalid_array[i] << 4));
    uint32_t count = __builtin_popcount(invalid_array[i]);
    curr_bucket->bitmap = curr_bucket->bitmap - count;
  }

  for (int i = 0; i < stashBucket; ++i) {
    auto curr_bucket = org_table->stash + i;
    curr_bucket->bitmap =
        curr_bucket->bitmap & (~(invalid_array[kNumBucket + i] << 18));
    uint32_t count = __builtin_popcount(invalid_array[kNumBucket + i]);
    curr_bucket->bitmap = curr_bucket->bitmap - count;
  }

  /*traverse the overflow list, also need to do the rehahsing to the original
   * bucket*/
  prev_bucket = org_table->stash;
  next_bucket = org_table->stash->next;
  while (next_bucket != NULL) {
    auto mask = GET_BITMAP(next_bucket->bitmap);
    for (int i = 0; i < kNumPairPerBucket; ++i) {
      if (CHECK_BIT(mask, i)) {
        if constexpr (std::is_pointer_v<T>) {
          key_hash =
              h(next_bucket->_[i].key->key, next_bucket->_[i].key->length);
        } else {
          key_hash = h(&(next_bucket->_[i].key), sizeof(Key_t));
        }

        auto x = key_hash % (2 * base_level);
        if (x >= base_level) {
          Insert4split(next_bucket->_[i].key, next_bucket->_[i].value, key_hash,
                       next_bucket->finger_array[i]);
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto org_bucket = org_table->bucket + bucket_ix;
          auto neighbor_bucket =
              org_table->bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->unset_indicator(next_bucket->finger_array[i],
                                      neighbor_bucket, next_bucket->_[i].key,
                                      3);
          next_bucket->unset_hash(i);
#ifdef COUNTING
          org_table->number--;
#endif
        } else {
          /*rehashing to original bucket*/
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto ret = org_table->Insert2Org(
              next_bucket->_[i].key, next_bucket->_[i].value, key_hash,
              bucket_ix); /*the unique check may avoid this restore*/
          if (ret == 0) {
            auto org_bucket = org_table->bucket + bucket_ix;
            auto neighbor_bucket =
                org_table->bucket + ((bucket_ix + 1) & bucketMask);
            org_bucket->unset_indicator(next_bucket->finger_array[i],
                                        neighbor_bucket, next_bucket->_[i].key,
                                        3);
            next_bucket->unset_hash(i);
          }
        }
      }
    }

    if (GET_COUNT(next_bucket->bitmap) == 0) {
      prev_bucket->next = next_bucket->next;
      Allocator::Free(next_bucket);
      next_bucket = prev_bucket->next;
    } else {
      prev_bucket = next_bucket;
      next_bucket = next_bucket->next;
    }
  }

#ifdef PMEM
  prev_bucket = org_table->stash;
  next_bucket = prev_bucket->next;

  while (next_bucket != NULL) {
    Allocator::Persist(next_bucket, sizeof(struct overflowBucket<T>));
    prev_bucket = next_bucket;
    next_bucket = next_bucket->next;
  }

  Allocator::Persist(org_table, sizeof(Table));
#endif

  org_table->state = 0;
  Allocator::Persist(&org_table->state, sizeof(org_table->state));
  state = 0;
  Allocator::Persist(&state, sizeof(state));
  for (int i = 0; i < kNumBucket; ++i) {
    curr_bucket = org_table->bucket + i;
    curr_bucket->release_lock();
  }
}

/*merge the neighbor table with current table*/
template <class T>
void Table<T>::Merge(Table<T> *neighbor, bool unique_check_flag) {
  if (unique_check_flag) {
    /*Restore the split/merge procedure*/
    size_t key_hash;
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = neighbor->bucket + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }

          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j],
                       true); /*this shceme may destory
                           the balanced segment*/
        }
      }
    }

    for (int i = 0; i < stashBucket; ++i) {
      auto *curr_bucket = neighbor->stash + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j],
                       true); /*this shceme may destory
                           the balanced segment*/
        }
      }
    }

    /*traverse the overflow list*/
    overflowBucket<T> *prev_bucket = neighbor->stash;
    overflowBucket<T> *curr_bucket = neighbor->stash->next;
    while (curr_bucket != NULL) {
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j],
                       true); /*this shceme may destory
                           the balanced segment*/
        }
      }
      prev_bucket = curr_bucket;
      curr_bucket = curr_bucket->next;
    }
  } else {
    /*Restore the split/merge procedure*/
    size_t key_hash;
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = neighbor->bucket + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }

          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]); /*this shceme may destory
                                                         the balanced segment*/
        }
      }
    }

    for (int i = 0; i < stashBucket; ++i) {
      auto *curr_bucket = neighbor->stash + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]); /*this shceme may destory
                                                         the balanced segment*/
        }
      }
    }

    /*traverse the overflow list*/
    overflowBucket<T> *prev_bucket = neighbor->stash;
    overflowBucket<T> *curr_bucket = neighbor->stash->next;
    while (curr_bucket != NULL) {
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]); /*this shceme may destory
                                                         the balanced segment*/
        }
      }
      prev_bucket = curr_bucket;
      curr_bucket = curr_bucket->next;
    }
  }
}

/* it needs to verify whether this bucket has been deleted...*/
template <class T>
int Table<T>::Insert(T key, Value_t value, size_t key_hash, Directory<T> *_dir,
                     uint64_t index, uint32_t old_N, uint32_t old_next) {
  /*we need to first do the locking and then do the verify*/
  uint8_t meta_hash = META_HASH(key_hash);
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);
  target->get_lock();
  if (!neighbor->try_get_lock()) {
    target->release_lock();
    return -2;
  }

  if (!target->test_initialize()) {
    neighbor->release_lock();
    target->release_lock();
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = bucket + i;
      curr_bucket->get_lock();
    }

    auto ret = verify_access(_dir, index, old_N, old_next);
    if (ret == -1) {
      for (int i = 0; i < kNumBucket; ++i) {
        Bucket<T> *curr_bucket = bucket + i;
        curr_bucket->release_lock();
      }
      return -2;
    }

    uint64_t org_idx;
    uint64_t base_level;
    /*org_idx is the index of the original table, the base_level is the index
     * diff between original table and target table*/
    Table<T> *org_table = get_org_table(index, &org_idx, &base_level, _dir);
    /*the split process splits from original table to target table*/
    Split(org_table, base_level, org_idx, _dir);
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = bucket + i;
      curr_bucket->release_lock();
    }
    return -2; /*return retry to reinsert the key-value*/
  } else {
    auto ret = verify_access(_dir, index, old_N, old_next);
    if (ret == -1) {
      neighbor->release_lock();
      target->release_lock();
      return -2;
    }

    /*the unique_check is to check whether the key has existed*/
    ret = target->unique_check(meta_hash, key, neighbor, stash);
    if (ret == -1) {
      neighbor->release_lock();
      target->release_lock();
      return -3; /* duplicate insert*/
    }

    int target_num = GET_COUNT(target->bitmap);
    int neighbor_num = GET_COUNT(neighbor->bitmap);
    if ((target_num == kNumPairPerBucket) &&
        (neighbor_num == kNumPairPerBucket)) {
      /* overflow handling */
      Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
      // Next displacement
      if (!next_neighbor->try_get_lock()) {
        neighbor->release_lock();
        target->release_lock();
        return -2;
      }
      auto ret = Next_displace(neighbor, next_neighbor, key, value, meta_hash);
      if (ret == 0) {
        target->release_lock();
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

      ret = Prev_displace(target, prev_neighbor, key, value, meta_hash);
      if (ret == 0) {
        neighbor->release_lock();
        return 0;
      }

      stash->get_lock();
      ret =
          Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);

      stash->release_lock();
      neighbor->release_lock();
      target->release_lock();
      prev_neighbor->release_lock();
      return ret;
    }

    if (target_num <= neighbor_num) {
      target->Insert(key, value, meta_hash, false);
      target->release_lock();
#ifdef PMEM
      Allocator::Persist(&target->bitmap, sizeof(target->bitmap));
#endif
      neighbor->release_lock();
    } else {
      neighbor->Insert(key, value, meta_hash, true);
      neighbor->release_lock();
#ifdef PMEM
      Allocator::Persist(&neighbor->bitmap, sizeof(neighbor->bitmap));
#endif
      target->release_lock();
    }

#ifdef COUNTING
    __sync_fetch_and_add(&number, 1);
#endif
    return 0;
  }
}

/*the insert needs to be perfectly balanced, not destory the power of balance*/
template <class T>
void Table<T>::Insert4split(T key, Value_t value, size_t key_hash,
                            uint8_t meta_hash) {
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);
  Bucket<T> *insert_target;
  bool probe = false;
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    insert_target = target;
  } else {
    insert_target = neighbor;
    probe = true;
  }

  /*some bucket may be overflowed?*/
  if (GET_COUNT(insert_target->bitmap) < kNumPairPerBucket) {
    insert_target->_[GET_COUNT(insert_target->bitmap)].key = key;
    insert_target->_[GET_COUNT(insert_target->bitmap)].value = value;
    insert_target->set_hash(GET_COUNT(insert_target->bitmap), meta_hash, probe);
#ifdef COUNTING
    ++number;
#endif
  } else {
    /*do the displacement or insertion in the stash*/
    Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
    int displace_index;
    displace_index = neighbor->Find_org_displacement();
    if (((GET_COUNT(next_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert_with_noflush(
          neighbor->_[displace_index].key, neighbor->_[displace_index].value,
          neighbor->finger_array[displace_index], true);
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace_with_noflush(key, value, meta_hash,
                                             displace_index, true);
#ifdef COUNTING
      ++number;
#endif
      return;
    }
    Bucket<T> *prev_neighbor;
    int prev_index;
    if (y == 0) {
      prev_neighbor = bucket + kNumBucket - 1;
      prev_index = kNumBucket - 1;
    } else {
      prev_neighbor = bucket + y - 1;
      prev_index = y - 1;
    }

    displace_index = target->Find_probe_displacement();
    if (((GET_COUNT(prev_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert_with_noflush(
          target->_[displace_index].key, target->_[displace_index].value,
          target->finger_array[displace_index], false);
      target->unset_hash(displace_index);
      target->Insert_displace_with_noflush(key, value, meta_hash,
                                           displace_index, false);
#ifdef COUNTING
      ++number;
#endif
      return;
    }

    Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);
  }
}

template <class T>
void Table<T>::Insert4merge(T key, Value_t value, size_t key_hash,
                            uint8_t meta_hash, bool unique_check_flag) {
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);

  if (unique_check_flag) {
    auto ret = target->unique_check(meta_hash, key, neighbor, stash);
    if (ret == -1) return;
  }
  Bucket<T> *insert_target;
  bool probe = false;
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    insert_target = target;
  } else {
    insert_target = neighbor;
    probe = true;
  }

  if (GET_COUNT(insert_target->bitmap) < kNumPairPerBucket) {
    insert_target->Insert(key, value, meta_hash, probe);
#ifdef COUNTING
    ++number;
#endif
  } else {
    /*do the displacement or insertion in the stash*/
    Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
    int displace_index;
    displace_index = neighbor->Find_org_displacement();
    if (((GET_COUNT(next_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert_with_noflush(
          neighbor->_[displace_index].key, neighbor->_[displace_index].value,
          neighbor->finger_array[displace_index], true);
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace_with_noflush(key, value, meta_hash,
                                             displace_index, true);
#ifdef COUNTING
      ++number;
#endif
      return;
    }
    Bucket<T> *prev_neighbor;
    int prev_index;
    if (y == 0) {
      prev_neighbor = bucket + kNumBucket - 1;
      prev_index = kNumBucket - 1;
    } else {
      prev_neighbor = bucket + y - 1;
      prev_index = y - 1;
    }

    displace_index = target->Find_probe_displacement();
    if (((GET_COUNT(prev_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert_with_noflush(
          target->_[displace_index].key, target->_[displace_index].value,
          target->finger_array[displace_index], false);
      target->unset_hash(displace_index);
      target->Insert_displace_with_noflush(key, value, meta_hash,
                                           displace_index, false);
#ifdef COUNTING
      ++number;
#endif
      return;
    }

    Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);
  }
}

template <class T>
class Linear : public Hash<T> {
 public:
  Linear(void);
  Linear(PMEMobjpool *_pool);
  ~Linear(void);
  int Insert(T key, Value_t value);
  bool Delete(T);
  int Insert(T key, Value_t value, bool);
  bool Delete(T, bool);
  inline Value_t Get(T);
  Value_t Get(T key, bool is_in_epoch);
  void FindAnyway(T key);
  void Recovery();
  void ShutDown() {
    clean = true;
    Allocator::Persist(&clean, sizeof(clean));
  }
  bool TryMerge(uint64_t, Table<T> *);
  void recoverSegment(Table<T> **seg_ptr, size_t, size_t, size_t);
  void getNumber() {
    uint64_t count = 0;
    uint64_t prev_length = 0;
    uint64_t after_length = 0;
    uint64_t Bucket_num = 0;
    uint64_t old_N_next = dir.N_next;
    uint32_t N = old_N_next >> 32;
    uint32_t next = (uint32_t)old_N_next;
    std::cout << "N = " << N << std::endl;
    std::cout << "next = " << next << std::endl;
    uint32_t occupied_bucket = pow2(N) + next;
    uint64_t recount_num = 0;
    uint32_t max_dir = 0;

    for (int i = 0; i < occupied_bucket; ++i) {
      uint32_t dir_idx;
      uint32_t offset;
      SEG_IDX_OFFSET(i, dir_idx, offset);
      Table<T> *curr_table = dir._[dir_idx] + offset;
      if (max_dir < dir_idx) max_dir = dir_idx;
#ifdef COUNTING
      recount_num += curr_table->number;
#endif
      for (int j = 0; j < kNumBucket; ++j) {
        Bucket<T> *curr_bucket = curr_table->bucket + j;
        count += GET_COUNT(curr_bucket->bitmap);
        int mask = GET_BITMAP(curr_bucket->bitmap);
        int micro_count = 0;
        if (mask != 0) {
          for (int k = 0; k < kNumPairPerBucket; ++k) {
            if (CHECK_BIT(mask, k)) {
              micro_count++;
            }
          }
        }
        assert(micro_count == GET_COUNT(curr_bucket->bitmap));
      }

      for (int j = 0; j < stashBucket; ++j) {
        overflowBucket<T> *curr_bucket = curr_table->stash + j;
        count += GET_COUNT(curr_bucket->bitmap);
        int mask = GET_BITMAP(curr_bucket->bitmap);
        int micro_count = 0;
        if (mask != 0) {
          for (int k = 0; k < kNumPairPerBucket; ++k) {
            if (CHECK_BIT(mask, k)) {
              micro_count++;
            }
          }
        }
        assert(micro_count == GET_COUNT(curr_bucket->bitmap));
      }

      overflowBucket<T> *prev_bucket = curr_table->stash;
      overflowBucket<T> *next_bucket = prev_bucket->next;
      while (next_bucket != NULL) {
        count += GET_COUNT(next_bucket->bitmap);
        int mask = GET_BITMAP(next_bucket->bitmap);
        int micro_count = 0;
        if (mask != 0) {
          for (int k = 0; k < kNumPairPerBucket; ++k) {
            if (CHECK_BIT(mask, k)) {
              micro_count++;
            }
          }
        }
        assert(micro_count == GET_COUNT(next_bucket->bitmap));
        prev_bucket = next_bucket;
        next_bucket = next_bucket->next;
        if (i < next) {
          prev_length++;
        } else {
          after_length++;
        }
        Bucket_num++;
      }
    }

    std::cout << "The # directory entries is " << max_dir << std::endl;
    std::cout << "Occupied table is " << occupied_bucket << std::endl;
    std::cout << "The size of overflow bucket is " << sizeof(overflowBucket<T>)
              << " ;The size of bucket is " << sizeof(Bucket<T>) << std::endl;
    Bucket_num += SUM_BUCKET(occupied_bucket - 1) * (kNumBucket + stashBucket);
    std::cout << "The size of table is that " << sizeof(Table<T>) << std::endl;
#ifdef COUNTING
    std::cout << "The recount number is " << recount_num << std::endl;
#endif
    std::cout << "the inserted num is " << count << std::endl;
    std::cout << "the bucket number is " << Bucket_num << std::endl;
    std::cout << "the local load factor = " << (double)count / (Bucket_num * kNumPairPerBucket) << std::endl;
    std::cout << "the local raw sapce utilization = " << (double)count / (Bucket_num * 16) << std::endl;
    std::cout << "the prev_length = " << prev_length << std::endl;
    std::cout << "the after_length = " << after_length << std::endl;
  }

  /**
   * @brief Expand operation
   * @param numBuckets the number of "Buckets" to expand
   * @return void
   */
  inline void Expand(uint32_t numBuckets) {
  RE_EXPAND:
    uint64_t old_N_next = dir.N_next;
    uint32_t old_N = old_N_next >> 32;
    uint32_t old_next = (uint32_t)old_N_next;
    uint32_t dir_idx, offset;
    SEG_IDX_OFFSET(
        static_cast<uint32_t>(pow(2, old_N)) + old_next + numBuckets - 1,
        dir_idx, offset);
    /*first need the reservation of the key-value*/
    Table<T> *RESERVED = reinterpret_cast<Table<T> *>(-1);
    if (dir._[dir_idx] == RESERVED) {
      goto RE_EXPAND;
    }

    Table<T> *old_value = NULL;
    if (dir._[dir_idx] == NULL) {
      /* Need to allocate the memory for new segment*/
      if (CAS(&(dir._[dir_idx]), &old_value, RESERVED)) {
#ifdef DOUBLE_EXPANSION
        uint32_t seg_size = SEG_SIZE(static_cast<uint32_t>(pow(2, old_N)) +
                                     old_next + numBuckets - 1);
#ifdef PREALLOC
        dir._[dir_idx] = TlsTablePool<T>::Get(seg_size);
#else
        Allocator::ZAllocate(&back_seg, kCacheLineSize,
                             sizeof(Table<T>) * seg_size);
        dir._[dir_idx] = reinterpret_cast<Table<T> *>(pmemobj_direct(back_seg));
        back_seg = OID_NULL;
#endif
#else
        Allocator::ZAllocate((void **)&dir._[dir_idx], kCacheLineSize,
                             sizeof(Table<T>) * segmentSize);
#endif
#ifdef PMEM
        Allocator::Persist(&dir._[dir_idx], sizeof(Table<T> *));
#endif
      } else {
        goto RE_EXPAND;
      }
    }

    uint64_t new_N_next;
    if (old_next == (static_cast<uint64_t>(pow(2, old_N)) - numBuckets)) {
      new_N_next = (((uint64_t)(old_N + 1)) << 32);
    } else {
      new_N_next = (old_N_next & high32Mask) + (old_next + numBuckets);
    }

    if (!CAS(&dir.N_next, &old_N_next, new_N_next)) {
      goto RE_EXPAND;
    }

#ifdef PMEM
    Allocator::Persist(&dir.N_next, sizeof(uint64_t));
#endif
    if ((uint32_t)new_N_next == 0) {
      printf("expand to level %lu\n", new_N_next >> 32);
    }
  }

#ifdef PMEM
  PMEMobjpool *pool_addr;
  PMEMoid back_seg;
#endif
  Directory<T> dir;
  int lock;
  bool clean;
};

template <class T>
Linear<T>::Linear(PMEMobjpool *_pool) {
  std::cout << "Start to initialize from scratch" << std::endl;
  pool_addr = _pool;
  lock = 0;
  clean = false;
  dir.N_next = baseShifBits << 32;
  std::cout << "Table size is " << sizeof(Table<T>) << std::endl;
  memset(dir._, 0, directorySize * sizeof(uint64_t));

  Allocator::ZAllocate((void **)&dir._[0], kCacheLineSize,
                       sizeof(Table<T>) * segmentSize);
  for (int j = 0; j < segmentSize; ++j) {
    Table<T> *curr_table = dir._[0] + j;
    for (int k = 0; k < kNumBucket; ++k) {
      Bucket<T> *curr_bucket = curr_table->bucket + k;
      curr_bucket->set_initialize();
    }
  }
}

template <class T>
Linear<T>::Linear(void) {
  std::cout << "Reinitialize Up for linear hashing" << std::endl;
}

template <class T>
Linear<T>::~Linear(void) {
  // TO-DO
}

template <class T>
bool Linear<T>::TryMerge(uint64_t x, Table<T> *shrunk_table) {
  /* Get all of the locks*/
  for (int i = 0; i < kNumBucket; ++i) {
    Bucket<T> *curr_bucket = shrunk_table->bucket + i;
    curr_bucket->get_lock();
  }

  if (shrunk_table->state != 0) {
    /* If its buddy is still in a unmerged state, just give up the merge
     * operation*/
    uint64_t idx, base_diff;
    Table<T> *expand_table =
        shrunk_table->get_expan_table(x, &idx, &base_diff, &dir);
    if ((expand_table != NULL) && (expand_table->state == 1)) {
      return false;
    }

    /*Get all the locks from original table*/
    Table<T> *org_table =
        shrunk_table->get_org_table(x, &idx, &base_diff, &dir);
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = org_table->bucket + i;
      curr_bucket->get_lock();
    }

    /*Do the merge operation*/
    org_table->Merge(shrunk_table);
    shrunk_table->state = 0;
  }

  for (int i = 0; i < kNumBucket; ++i) {
    Bucket<T> *curr_bucket = shrunk_table->bucket + i;
    curr_bucket->release_lock();
  }

  /*Try to reset the N_next to the normal state*/
  uint64_t old_N_next = dir->N_next;
  uint32_t right_index = pow2(old_N_next >> 32) + (uint32_t)old_N_next;

  /*Reset the N_next to the new state*/
  /*recursively merge the (org)table */
  return true;
}

template <class T>
void Linear<T>::Recovery() {
  if (clean) {
    clean = false;
    return;
  }
  Allocator::EpochRecovery();
  uint64_t old_N_next = dir.N_next;
  uint32_t N = old_N_next >> 32;
  uint32_t next = (uint32_t)old_N_next;
  uint32_t x = static_cast<uint32_t>(pow(2, N)) + next - 1;
  dir.recovered_index =
      x; /* for segments <= recovered_index, it needs to be recoverd*/
  uint32_t dir_idx;
  uint32_t offset;
  SEG_IDX_OFFSET(x, dir_idx, offset);

  for (int i = 0; i < dir_idx; ++i) {
    dir.recover_counter[i] = SEG_SIZE_BY_SEGARR_ID(i);
    dir._[i] = reinterpret_cast<Table<T> *>((uint64_t)dir._[i] | recoverBit);
  }
  std::cout << dir_idx << " segments array in the linear hashing" << std::endl;

  dir.recover_counter[dir_idx] = offset + 1;
  dir._[dir_idx] =
      reinterpret_cast<Table<T> *>((uint64_t)dir._[dir_idx] | recoverBit);

  dir.crash_version += 1;
  if (dir.crash_version == 0) {
    /* Scan all the segments to make it invalid*/
    uint32_t occupied_bucket = pow2(N) + next;
    uint64_t recount_num = 0;
    for (int i = 0; i < occupied_bucket; ++i) {
      uint32_t dir_idx;
      uint32_t offset;
      SEG_IDX_OFFSET(i, dir_idx, offset);
      Table<T> *curr_table = dir._[dir_idx] + offset;
      curr_table->seg_version = 1;
    }
  }
}

template <class T>
void Linear<T>::recoverSegment(Table<T> **seg_ptr, size_t index, size_t dir_idx,
                               size_t offset) {
RETRY:
  uint64_t snapshot = reinterpret_cast<uint64_t>(*seg_ptr);
  Table<T> *target = (Table<T> *)(snapshot & (~recoverLockBit)) + offset;

  /*No need for the recovery of this segment*/
  if ((dir.crash_version == target->seg_version) ||
      (index > dir.recovered_index)) {
    return;
  }

  /*try to get the exclusive recovery lock*/
  if (pmemobj_mutex_trylock(pool_addr, &target->lock_bit) != 0) {
    goto RETRY;
  }
  target->recoverMetadata();

  /*FIXME: handle state = 1*/
  if (target->state == 2) {
    uint64_t idx, base_diff;
    Table<T> *org_table = target->get_org_table(index, &idx, &base_diff, &dir);
    uint32_t buddy_dir_idx, buddy_offset;
    SEG_IDX_OFFSET(idx, buddy_dir_idx, buddy_offset);
    while (reinterpret_cast<uint64_t>(dir._[buddy_dir_idx]) & recoverLockBit) {
      recoverSegment(&dir._[buddy_dir_idx], idx, buddy_dir_idx, buddy_offset);
    }

    for (int i = 0; i < kNumBucket; ++i) {
      auto curr_bucket = org_table->bucket + i;
      curr_bucket->get_lock();
    }

    org_table->Merge(target, true);
    for (int i = 0; i, kNumBucket; ++i) {
      auto curr_bucket = target->bucket + i;
      curr_bucket->unset_initialize();
    }

    target->state = 0;
    Allocator::Persist(&target->state, sizeof(target->state));
    org_table->state = 0;
    Allocator::Persist(&org_table->state, sizeof(org_table->state));

    for (int i = 0; i < kNumBucket; ++i) {
      auto curr_bucket = org_table->bucket + i;
      curr_bucket->release_lock();
    }
  }

  target->seg_version = dir.crash_version;
  SUB(&dir.recover_counter[dir_idx], 1);
  if (dir.recover_counter[dir_idx] <= 0) {
    *seg_ptr = (Table<T> *)(snapshot & (~recoverLockBit));
  }
}

template <class T>
int Linear<T>::Insert(T key, Value_t value, bool is_in_epoch) {
  if (!is_in_epoch) {
    auto epoch_guard = Allocator::AquireEpochGuard();
    return Insert(key, value);
  }
  return Insert(key, value);
}

template <class T>
int Linear<T>::Insert(T key, Value_t value) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
RETRY:
  uint64_t old_N_next = dir.N_next;
  uint32_t N = old_N_next >> 32;
  uint32_t next = (uint32_t)old_N_next;
  auto x = IDX(key_hash, N);
  if (x < next) {
    x = IDX(key_hash, N + 1);
  }

  uint32_t dir_idx;
  uint32_t offset;
  SEG_IDX_OFFSET(static_cast<uint32_t>(x), dir_idx, offset);
  Table<T> *target = dir._[dir_idx] + offset;
  if (reinterpret_cast<uint64_t>(dir._[dir_idx]) & recoverLockBit) {
    recoverSegment(&dir._[dir_idx], x, dir_idx, offset);
    target =
        (Table<T> *)((uint64_t)(dir._[dir_idx]) & (~recoverLockBit)) + offset;
  }

  auto ret = target->Insert(key, value, key_hash, &dir, x, N, next);

  if (ret == -2) {
    goto RETRY;
  } else if (ret == -1) {
    Expand(2);
  } else if (ret == -3){
    return -1;
  }

  return 0;
}

template <class T>
Value_t Linear<T>::Get(T key, bool is_in_epoch) {
  if (!is_in_epoch) {
    auto epoch_guard = Allocator::AquireEpochGuard();
    return Get(key);
  }

  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = META_HASH(key_hash);
  auto y = BUCKET_INDEX(key_hash);
RETRY:
  uint64_t old_N_next = dir.N_next;
  uint32_t N = old_N_next >> 32;
  uint32_t next = (uint32_t)old_N_next;

  auto x = IDX(key_hash, N);
  if (x < next) {
    x = IDX(key_hash, N + 1);
  }

  uint32_t dir_idx;
  uint32_t offset;
  SEG_IDX_OFFSET(static_cast<uint32_t>(x), dir_idx, offset);
  Table<T> *target = dir._[dir_idx] + offset;

  if (reinterpret_cast<uint64_t>(dir._[dir_idx]) & recoverLockBit) {
    recoverSegment(&dir._[dir_idx], x, dir_idx, offset);
    target =
        (Table<T> *)((uint64_t)(dir._[dir_idx]) & (~recoverLockBit)) + offset;
  }

  Bucket<T> *target_bucket = target->bucket + y;
  Bucket<T> *neighbor_bucket = target->bucket + ((y + 1) & bucketMask);
  uint32_t old_version = target_bucket->version_lock;
  uint32_t old_neighbor_version = neighbor_bucket->version_lock;

  if ((old_version & lockSet) || (old_neighbor_version & lockSet)) {
    goto RETRY;
  }

  if (old_version & initialSet) {
    uint64_t new_N_next = dir.N_next;
    uint32_t new_N = new_N_next >> 32;
    uint32_t new_next = (uint32_t)new_N_next;
    if (((next <= x) && (new_next > x)) || (new_N != N)) {
      goto RETRY;
    }

    auto ret = target_bucket->check_and_get(meta_hash, key, false);
    if (target_bucket->test_lock_version_change(old_version)) {
      goto RETRY;
    }

    if (ret != NONE) {
      return ret;
    }

    /*no need for verification procedure, we use the version number of
     * target_bucket to test whether the bucket has ben spliteted*/
    ret = neighbor_bucket->check_and_get(meta_hash, key, true);
    if (neighbor_bucket->test_lock_version_change(old_neighbor_version)) {
      goto RETRY;
    }
    if (ret != NONE) {
      return ret;
    }

    if (target_bucket->test_stash_check()) {
      auto test_stash = false;
      if (target_bucket->test_overflow()) {
        test_stash = true;
      } else {
        // search in the original bucket
        int mask = target_bucket->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (target_bucket->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & target_bucket->overflowMember) == 0)) {
              test_stash = true;
              goto TEST_STASH;
            }
          }
        }

        mask = neighbor_bucket->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (neighbor_bucket->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & neighbor_bucket->overflowMember) != 0)) {
              test_stash = true;
              break;
            }
          }
        }
      }
    TEST_STASH:
      if (test_stash == true) {
        for (int i = 0; i < stashBucket; ++i) {
          overflowBucket<T> *curr_bucket =
              target->stash + ((i + (y & stashMask)) & stashMask);
          auto ret = curr_bucket->check_and_get(meta_hash, key);
          if (ret != NONE) {
            if (target_bucket->test_lock_version_change(old_version)) {
              goto RETRY;
            }
            return ret;
          }
        }

        overflowBucket<T> *prev_bucket = target->stash;
        overflowBucket<T> *next_bucket = target->stash->next;
        while (next_bucket != NULL) {
          auto ret = next_bucket->check_and_get(meta_hash, key);
          if (ret != NONE) {
            if (target_bucket->test_lock_version_change(old_version)) {
              goto RETRY;
            }
            return ret;
          }
          prev_bucket = next_bucket;
          next_bucket = next_bucket->next;
        }
      }
    }
  } else {
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = target->bucket + i;
      curr_bucket->get_lock();
    }

    uint64_t new_N_next = dir.N_next;
    uint32_t new_N = new_N_next >> 32;
    uint32_t new_next = (uint32_t)new_N_next;
    if (((next <= x) && (new_next > x)) || (new_N != N)) {
      for (int i = 0; i < kNumBucket; ++i) {
        Bucket<T> *curr_bucket = target->bucket + i;
        curr_bucket->release_lock();
      }
      goto RETRY;
    }

    uint64_t org_idx;
    uint64_t base_level;
    Table<T> *org_table = target->get_org_table(x, &org_idx, &base_level, &dir);
    target->Split(org_table, base_level, org_idx, &dir);

    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = target->bucket + i;
      curr_bucket->release_lock();
    }
    goto RETRY;
  }
  return NONE;
}

template <class T>
Value_t Linear<T>::Get(T key) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = META_HASH(key_hash);
  auto y = BUCKET_INDEX(key_hash);
RETRY:
  uint64_t old_N_next = dir.N_next;
  uint32_t N = old_N_next >> 32;
  uint32_t next = (uint32_t)old_N_next;

  auto x = IDX(key_hash, N);
  if (x < next) {
    x = IDX(key_hash, N + 1);
  }

  uint32_t dir_idx;
  uint32_t offset;
  SEG_IDX_OFFSET(static_cast<uint32_t>(x), dir_idx, offset);
  Table<T> *target = dir._[dir_idx] + offset;

  if (reinterpret_cast<uint64_t>(dir._[dir_idx]) & recoverLockBit) {
    recoverSegment(&dir._[dir_idx], x, dir_idx, offset);
    target =
        (Table<T> *)((uint64_t)(dir._[dir_idx]) & (~recoverLockBit)) + offset;
  }

  Bucket<T> *target_bucket = target->bucket + y;
  Bucket<T> *neighbor_bucket = target->bucket + ((y + 1) & bucketMask);
  uint32_t old_version = target_bucket->version_lock;
  uint32_t old_neighbor_version = neighbor_bucket->version_lock;

  if ((old_version & lockSet) || (old_neighbor_version & lockSet)) {
    goto RETRY;
  }

  if (old_version & initialSet) {
    uint64_t new_N_next = dir.N_next;
    uint32_t new_N = new_N_next >> 32;
    uint32_t new_next = (uint32_t)new_N_next;
    if (((next <= x) && (new_next > x)) || (new_N != N)) {
      goto RETRY;
    }

    auto ret = target_bucket->check_and_get(meta_hash, key, false);
    if (target_bucket->test_lock_version_change(old_version)) {
      goto RETRY;
    }

    if (ret != NONE) {
      return ret;
    }

    /*no need for verification procedure, we use the version number of
     * target_bucket to test whether the bucket has ben spliteted*/
    ret = neighbor_bucket->check_and_get(meta_hash, key, true);
    if (neighbor_bucket->test_lock_version_change(old_neighbor_version)) {
      goto RETRY;
    }
    if (ret != NONE) {
      return ret;
    }

    if (target_bucket->test_stash_check()) {
      auto test_stash = false;
      if (target_bucket->test_overflow()) {
        // this only occur when the bucket has more key-values than 10 that are
        // overfloed int he shared bucket area, therefore it needs to search in
        // the extra bucket
        test_stash = true;
      } else {
        // search in the original bucket
        int mask = target_bucket->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (target_bucket->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & target_bucket->overflowMember) == 0)) {
              test_stash = true;
              goto TEST_STASH;
            }
          }
        }

        mask = neighbor_bucket->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (neighbor_bucket->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & neighbor_bucket->overflowMember) != 0)) {
              test_stash = true;
              break;
            }
          }
        }
      }
    TEST_STASH:
      if (test_stash == true) {
        for (int i = 0; i < stashBucket; ++i) {
          overflowBucket<T> *curr_bucket =
              target->stash + ((i + (y & stashMask)) & stashMask);
          auto ret = curr_bucket->check_and_get(meta_hash, key);
          if (ret != NONE) {
            if (target_bucket->test_lock_version_change(old_version)) {
              goto RETRY;
            }
            return ret;
          }
        }

        overflowBucket<T> *prev_bucket = target->stash;
        overflowBucket<T> *next_bucket = target->stash->next;
        while (next_bucket != NULL) {
          auto ret = next_bucket->check_and_get(meta_hash, key);
          if (ret != NONE) {
            if (target_bucket->test_lock_version_change(old_version)) {
              goto RETRY;
            }
            return ret;
          }
          prev_bucket = next_bucket;
          next_bucket = next_bucket->next;
        }
      }
    }
  } else {
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = target->bucket + i;
      curr_bucket->get_lock();
    }

    uint64_t new_N_next = dir.N_next;
    uint32_t new_N = new_N_next >> 32;
    uint32_t new_next = (uint32_t)new_N_next;
    if (((next <= x) && (new_next > x)) || (new_N != N)) {
      for (int i = 0; i < kNumBucket; ++i) {
        Bucket<T> *curr_bucket = target->bucket + i;
        curr_bucket->release_lock();
      }
      goto RETRY;
    }

    uint64_t org_idx;
    uint64_t base_level;
    Table<T> *org_table = target->get_org_table(x, &org_idx, &base_level, &dir);
    target->Split(org_table, base_level, org_idx, &dir);

    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = target->bucket + i;
      curr_bucket->release_lock();
    }
    goto RETRY;
  }
  return NONE;
}

template <class T>
bool Linear<T>::Delete(T key, bool is_in_epoch) {
  if (!is_in_epoch) {
    auto epoch_guard = Allocator::AquireEpochGuard();
    return Delete(key);
  }
  return Delete(key);
}

/*Current version of Dash linear hashing does not support concurrent shrink
 * operation*/
template <class T>
bool Linear<T>::Delete(T key) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = META_HASH(key_hash);
  auto y = BUCKET_INDEX(key_hash);
RETRY:
  uint64_t old_N_next = dir.N_next;
  uint32_t N = old_N_next >> 32;
  uint32_t next = (uint32_t)old_N_next;

  auto x = IDX(key_hash, N);
  if (x < next) {
    x = IDX(key_hash, N + 1);
  }

  uint32_t dir_idx;
  uint32_t offset;
  SEG_IDX_OFFSET(static_cast<uint32_t>(x), dir_idx, offset);
  Table<T> *target = dir._[dir_idx] + offset;
  if (reinterpret_cast<uint64_t>(dir._[dir_idx]) & recoverLockBit) {
    recoverSegment(&dir._[dir_idx], x, dir_idx, offset);
    target =
        (Table<T> *)((uint64_t)(dir._[dir_idx]) & (~recoverLockBit)) + offset;
  }

  uint32_t old_version;
  Bucket<T> *target_bucket = target->bucket + y;
  Bucket<T> *neighbor_bucket = target->bucket + ((y + 1) & bucketMask);

  target_bucket->get_lock();
  if (!neighbor_bucket->try_get_lock()) {
    target_bucket->release_lock();
    goto RETRY;
  }

  if (!target_bucket->test_initialize()) {
    target_bucket->release_lock();
    neighbor_bucket->release_lock();
    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = target->bucket + i;
      curr_bucket->get_lock();
    }

    uint64_t new_N_next = dir.N_next;
    uint32_t new_N = new_N_next >> 32;
    uint32_t new_next = (uint32_t)new_N_next;
    if (((next <= x) && (new_next > x)) || (new_N != N)) {
      for (int i = 0; i < kNumBucket; ++i) {
        Bucket<T> *curr_bucket = target->bucket + i;
        curr_bucket->release_lock();
      }
      goto RETRY;
    }

    uint64_t org_idx;
    uint64_t base_level;
    Table<T> *org_table = target->get_org_table(x, &org_idx, &base_level, &dir);
    target->Split(org_table, base_level, org_idx, &dir);

    for (int i = 0; i < kNumBucket; ++i) {
      Bucket<T> *curr_bucket = target->bucket + i;
      curr_bucket->release_lock();
    }
    goto RETRY;
  } else {
    uint64_t new_N_next = dir.N_next;
    uint32_t new_N = new_N_next >> 32;
    uint32_t new_next = (uint32_t)new_N_next;
    if (((next <= x) && (new_next > x)) || (new_N != N)) {
      target_bucket->release_lock();
      neighbor_bucket->release_lock();
      goto RETRY;
    }

    auto ret = target_bucket->Delete(meta_hash, key, false);
    if (ret == 0) {
#ifdef COUNTING
      auto num = SUB(&target->number, 1);
#endif
      target_bucket->release_lock();
#ifdef PMEM
      Allocator::Persist(&target_bucket->bitmap, sizeof(target_bucket->bitmap));
#endif
      neighbor_bucket->release_lock();
      return true;
    }

    /*no need for verification procedure, we use the version number of
     * target_bucket to test whether the bucket has ben spliteted*/
    ret = neighbor_bucket->Delete(meta_hash, key, true);
    if (ret == 0) {
#ifdef COUNTING
      auto num = SUB(&target->number, 1);
#endif
      neighbor_bucket->release_lock();
#ifdef PMEM
      Allocator::Persist(&neighbor_bucket->bitmap,
                         sizeof(neighbor_bucket->bitmap));
#endif
      target_bucket->release_lock();
      return true;
    }

    if (target_bucket->test_stash_check()) {
      auto test_stash = false;
      if (target_bucket->test_overflow()) {
        // this only occur when the bucket has more key-values than 10 that are
        // overfloed int he shared bucket area, therefore it needs to search in
        // the extra bucket
        test_stash = true;
      } else {
        // search in the original bucket
        int mask = target_bucket->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (target_bucket->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & target_bucket->overflowMember) == 0)) {
              test_stash = true;
              goto TEST_STASH;
            }
          }
        }

        mask = neighbor_bucket->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (neighbor_bucket->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & neighbor_bucket->overflowMember) != 0)) {
              test_stash = true;
              break;
            }
          }
        }
      }
    TEST_STASH:
      if (test_stash == true) {
        overflowBucket<T> *stash = target->stash;
        stash->get_lock();
        for (int i = 0; i < stashBucket; ++i) {
          int index = ((i + (y & stashMask)) & stashMask);
          overflowBucket<T> *curr_bucket = target->stash + index;
          auto ret = curr_bucket->Delete(meta_hash, key);
          if (ret == 0) {
#ifdef COUNTING
            auto num = SUB(&target->number, 1);
#endif
            stash->release_lock();
#ifdef PMEM
            Allocator::Persist(&curr_bucket->bitmap,
                               sizeof(curr_bucket->bitmap));
#endif
            target_bucket->unset_indicator(meta_hash, neighbor_bucket, key,
                                           index);
            target_bucket->release_lock();
            neighbor_bucket->release_lock();
            return true;
          }
        }

        overflowBucket<T> *prev_bucket = target->stash;
        overflowBucket<T> *next_bucket = target->stash->next;
        while (next_bucket != NULL) {
          auto ret = next_bucket->Delete(meta_hash, key);
          if (ret == 0) {
#ifdef COUNTING
            auto num = SUB(&target->number, 1);
#endif
            stash->release_lock();
#ifdef PMEM
            Allocator::Persist(&next_bucket->bitmap,
                               sizeof(next_bucket->bitmap));
#endif
            target_bucket->unset_indicator(meta_hash, neighbor_bucket, key, 3);
            target_bucket->release_lock();
            neighbor_bucket->release_lock();
            return true;
          }
          prev_bucket = next_bucket;
          next_bucket = next_bucket->next;
        }
        stash->release_lock();
      }
    }
  }
  target_bucket->release_lock();
  neighbor_bucket->release_lock();
  return false;
}

#undef PARTITION_INDEX
#undef BUCKET_INDEX
#undef META_HASH
#undef GET_COUNT
#undef GET_BITMAP
}  // namespace linear
