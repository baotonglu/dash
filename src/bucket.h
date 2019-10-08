#pragma once
#include <immintrin.h>
#include <bitset>
#include <cassert>
#include "../util/hash.h"
#include "../util/pair.h"
#include "../util/persist.h"
#include "allocator.h"

#ifdef PMEM
#include <libpmemobj.h>
#endif

#define _INVALID 0 /* we use 0 as the invalid key*/
#define SINGLE 1
#define COUNTING 1
#define EPOCH 1

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

template <class T>
struct _Pair {
  T key;
  Value_t value;
};

constexpr size_t k_PairSize = 16;  // a k-v _Pair with a bit
constexpr size_t kNumPairPerBucket =
    14; /* it is determined by the usage of the fingerprint*/
constexpr size_t kFingerBits = 8;
constexpr size_t kMask = (1 << kFingerBits) - 1;
const constexpr size_t kNumBucket = 64;
constexpr size_t stashBucket = 2;
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
constexpr size_t bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
constexpr size_t stashMask = (1 << (int)log2(stashBucket)) - 1;
constexpr uint8_t stashHighMask = ~((uint8_t)stashMask);
constexpr uint8_t preNeighborSet = 1 << 7;
constexpr uint8_t nextNeighborSet = 1 << 6;
const uint64_t recoverBit = 1UL << 63;
const uint64_t lockBit = 1UL << 62;
#define BUCKET_INDEX(hash) ((hash >> kFingerBits) & bucketMask)
#define GET_COUNT(var) ((var)&countMask)
#define GET_BITMAP(var) (((var) >> 4) & allocMask)
#define ORG_BITMAP(var) ((~((var)&allocMask)) & allocMask)
#define PROBE_BITMAP(var) ((var)&allocMask)

inline bool var_compare(char *str1, char *str2, int len1, int len2) {
  if (len1 != len2) return false;
  return !memcmp(str1, str2, len1);
}

template <class T>
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
  inline bool test_overflow() {
    // return (overflowCount != 0) ? true : false;
    return overflowCount;
  }

  inline bool test_stash_check() {
    int mask = *((int *)membership);
    // return ((mask & overflowSet) != 0) ? true : false;
    return (mask & overflowSet);
  }

  inline void clear_stash_check() {
    int mask = *((int *)membership);
    *((int *)membership) = (*((int *)membership)) & (~overflowSet);
  }

  inline void set_indicator(uint8_t meta_hash, Bucket<T> *neighbor,
                            uint8_t pos) {
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
  inline void unset_indicator(uint8_t meta_hash, Bucket<T> *neighbor, T key,
                              uint64_t pos) {
    /*also needs to ensure that this meta_hash must belongs to other bucket*/
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
      // printf("decrease overflowCount\n");
      overflowCount--;
    }

    mask1 = finger_array[14];
    mask2 = neighbor->finger_array[14];
    if (((mask1 & (~overflowMember)) == 0) && (overflowCount == 0) &&
        ((mask2 & neighbor->overflowMember) == 0)) {
      // printf("clear the stash check\n");
      clear_stash_check();
    }
  }

  int unique_check(uint8_t meta_hash, T key, Bucket<T> *neighbor,
                   Bucket<T> *stash) {
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
    /*
    for (int i = 0; i < stashBucket; ++i) {
          Bucket *curr_bucket = stash + i;
      if (curr_bucket->check_and_get(meta_hash, key, false) != NONE) {
        return -1;
      }
    }*/
    return 0;
  }

  inline int get_current_mask() {
    int mask = GET_BITMAP(bitmap) & ((~(*(int *)membership)) & allocMask);
    return mask;
  }

  Value_t check_and_get(uint8_t meta_hash, T key, bool probe) {
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & ((~(*(int *)membership)) & allocMask);
    } else {
      mask = mask & GET_BITMAP(bitmap) & ((*(int *)membership) & allocMask);
    }

    if (mask == 0) {
      return NONE;
    }

    if constexpr (std::is_pointer_v<T>) {
      string_key *_key = reinterpret_cast<string_key *>(key);
      for (int i = 0; i < 12; i += 4) {
        // if (CHECK_BIT(mask, i) && (strcmp(_[i].key, key) == 0)) {
        //  return _[i].value;
        //}
        if (CHECK_BIT(mask, i) &&
            (var_compare((reinterpret_cast<string_key *>(_[i].key))->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[i].key))->length,
                         _key->length))) {
          return _[i].value;
        }

        if (CHECK_BIT(mask, i + 1) &&
            (var_compare((reinterpret_cast<string_key *>(_[i + 1].key))->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[i + 1].key))->length,
                         _key->length))) {
          return _[i + 1].value;
        }

        if (CHECK_BIT(mask, i + 2) &&
            (var_compare((reinterpret_cast<string_key *>(_[i + 2].key))->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[i + 2].key))->length,
                         _key->length))) {
          return _[i + 2].value;
        }

        if (CHECK_BIT(mask, i + 3) &&
            (var_compare((reinterpret_cast<string_key *>(_[i + 3].key))->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[i + 3].key))->length,
                         _key->length))) {
          return _[i + 3].value;
        }
      }

      if (CHECK_BIT(mask, 12) &&
          (var_compare((reinterpret_cast<string_key *>(_[12].key))->key,
                       _key->key,
                       (reinterpret_cast<string_key *>(_[12].key))->length,
                       _key->length))) {
        return _[12].value;
      }

      if (CHECK_BIT(mask, 13) &&
          (var_compare((reinterpret_cast<string_key *>(_[13].key))->key,
                       _key->key,
                       (reinterpret_cast<string_key *>(_[13].key))->length,
                       _key->length))) {
        return _[13].value;
      }
    } else {
      /*loop unrolling*/
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
    // #ifdef PMEM
    // Allocator::NTWrite32(reinterpret_cast<uint32_t *>(&bitmap), new_bitmap);
    // #else
    // bitmap = new_bitmap;
    // #endif
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
#ifdef PMEM
    if (nt_flush) {
      Allocator::NTWrite32(reinterpret_cast<uint32_t *>(&bitmap), new_bitmap);
    } else {
      bitmap = new_bitmap;
    }
#else
    bitmap = new_bitmap;
#endif

    *((int *)membership) =
        (~(1 << index)) &
        (*((int *)membership)); /*since they are in the same cacheline,
                                   therefore no performance influence?*/
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
    // auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    // version = value & lockMask;
    // return (value & lockSet) != 0;
    version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    return (version & lockSet) != 0;
  }

  // test whether the version has change, if change, return true
  inline bool test_lock_version_change(uint32_t old_version) {
    // auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    // return ((value & lockSet) != 0) || ((value & lockMask) != old_version);
    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    return (old_version != value);
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
#ifdef PMEM
    Allocator::Persist(&_[slot], sizeof(_[slot]));
#endif
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
    if constexpr (std::is_pointer_v<T>) {
      string_key *_key = reinterpret_cast<string_key *>(key);
      /*loop unrolling*/
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) &&
              (var_compare((reinterpret_cast<string_key *>(_[i].key))->key,
                           _key->key,
                           (reinterpret_cast<string_key *>(_[i].key))->length,
                           _key->length))) {
            unset_hash(i, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) &&
              (var_compare(
                  reinterpret_cast<string_key *>(_[i + 1].key)->key, _key->key,
                  (reinterpret_cast<string_key *>(_[i + 1].key))->length,
                  _key->length))) {
            unset_hash(i + 1, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) &&
              (var_compare(
                  reinterpret_cast<string_key *>(_[i + 2].key)->key, _key->key,
                  (reinterpret_cast<string_key *>(_[i + 2].key))->length,
                  _key->length))) {
            unset_hash(i + 2, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) &&
              (var_compare(
                  reinterpret_cast<string_key *>(_[i + 3].key)->key, _key->key,
                  (reinterpret_cast<string_key *>(_[i + 3].key))->length,
                  _key->length))) {
            unset_hash(i + 3, false);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) &&
            (var_compare(reinterpret_cast<string_key *>(_[12].key)->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[12].key))->length,
                         _key->length))) {
          unset_hash(12, false);
          return 0;
        }

        if (CHECK_BIT(mask, 13) &&
            (var_compare(reinterpret_cast<string_key *>(_[13].key)->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[13].key))->length,
                         _key->length))) {
          unset_hash(13, false);
          return 0;
        }
      }

    } else {
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

  int Insert_with_noflush(T key, Value_t value, uint8_t meta_hash, bool probe) {
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
#ifdef PMEM
    Allocator::Persist(&_[slot], sizeof(_Pair<T>));
#endif
    mfence();
    set_hash(slot, meta_hash, probe);
  }

  void Insert_displace_with_noflush(T key, Value_t value, uint8_t meta_hash,
                                    int slot, bool probe) {
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
  }

  /* Find the displacment element in this bucket*/
  /*
  int Find_displacement(int x){
          for (int i = 0; i < kNumPairPerBucket; ++i)
          {
                  auto key_hash = h(&_[i], sizeof(Key_t));
                  auto y = BUCKET_INDEX(key_hash);
                  if (x == y)
                  {
                          return i;
                  }
          }
          return -1;
  }*/

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
  inline void setPreNonFlush() {
    overflowMember = overflowMember | preNeighborSet;
  }

  inline void unsetPreNonFlush() {
    overflowMember = overflowMember & (~preNeighborSet);
  }

  /* If the non-flush is set, return true*/
  inline bool testPreNonFlush() {
    return ((overflowMember & preNeighborSet) != 0);
  }

  inline void setNextNonFlush() {
    overflowMember = overflowMember | nextNeighborSet;
  }

  inline void unsetNextNonFlush() {
    overflowMember = overflowMember & (~nextNeighborSet);
  }

  inline bool testNextNonFlush() {
    return ((overflowMember & nextNeighborSet) != 0);
  }

  inline void resetLock() { version_lock = 0; }

  inline void resetOverflowFP() {
    finger_array[18] = 0;
    finger_array[19] = 0;
    overflowMember = 0;
    overflowCount = 0;
    clear_stash_check();
  }

  uint32_t version_lock;
  int bitmap;               // allocation bitmap + pointer bitmao + counter
  uint8_t finger_array[20]; /*only use the first 14 bytes, can be accelerated by
                               SSE instruction,0-13 for finger, 14-17 for
                               overflowed, 18 as the bitmap, 19 as the overflow
                               bucket index indicator*/
  uint8_t membership[2];    /*Used to test whether the key originally belongs to
                               this bucket*/
  uint8_t overflowMember; /*overflowmember indicates membership of the overflow
                             fingerprint*/
  uint8_t overflowCount;

  _Pair<T> _[kNumPairPerBucket];
};