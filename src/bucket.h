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

/*
const int overflowSet = 1 << 15;
const int countMask = (1 << 4) - 1;
*/
const uint32_t lockSet = 1 << 31;                  /*locking information*/
const uint32_t lockMask = ((uint32_t)1 << 31) - 1; /*locking mask*/
const int overflowSet = 1 << 15;
const int countMask = (1 << 4) - 1;
const uint32_t initialSet = 1 << 30;
const uint32_t initialLockSet = lockSet | initialSet;
const uint32_t versionMask = (1 << 30) - 1;
uint64_t overflow_access;

template <class T>
struct _Pair {
  T key;
  Value_t value;
};

constexpr size_t k_PairSize = 16;
const size_t kNumPairPerBucket = 14;
const size_t kFingerBits = 8;
constexpr size_t kMask = (1 << kFingerBits) - 1;
const uint32_t kNumBucket = 64;
const uint32_t stashBucket = 2;
const uint32_t fixedExpandNum = 8;
const uint32_t fixedExpandBits = 31 - __builtin_clz(fixedExpandNum);
const uint32_t fixedExpandMask = (1 << fixedExpandBits) - 1;
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
// constexpr size_t bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
constexpr size_t bucketMask = ((1 << (31 - __builtin_clz(kNumBucket))) - 1);
// constexpr size_t stashMask = (1 << (int)log2(stashBucket)) -1;
constexpr size_t stashMask = (1 << (31 - __builtin_clz(stashBucket))) - 1;
constexpr uint8_t stashHighMask = ~((uint8_t)stashMask);
constexpr uint32_t segmentSize = 64;
// constexpr size_t baseShifBits = static_cast<uint64_t>(log2(segmentSize));
constexpr size_t baseShifBits =
    static_cast<uint64_t>(31 - __builtin_clz(segmentSize));
constexpr uint32_t segmentMask = (1 << baseShifBits) - 1;
constexpr size_t directorySize = 1024 * 2;
const uint64_t low32Mask = ((uint64_t)1 << 32) - 1;
const uint64_t high32Mask = ~low32Mask;
const uint8_t low2Mask = (1 << 2) - 1;
const uint32_t low4Mask = (1 << 4) - 1;
// constexpr size_t shiftBits = (size_t)log2(kNumBucket) + kFingerBits;
constexpr uint32_t shiftBits = (31 - __builtin_clz(kNumBucket)) + kFingerBits;
const uint32_t partitionNum = 1;
constexpr uint64_t partitionMask =
    ((1 << (31 - __builtin_clz(partitionNum))) - 1);
constexpr uint32_t partitionShifBits =
    shiftBits + (31 - __builtin_clz(partitionNum));
constexpr uint32_t expandShiftBits = fixedExpandBits + baseShifBits;
/*
constexpr size_t k_PairSize = 16;  // a k-v _Pair with a bit
constexpr size_t kNumPairPerBucket =
    14;
constexpr size_t kFingerBits = 8;
constexpr size_t kMask = (1 << kFingerBits) - 1;
const constexpr size_t kNumBucket = 64;
constexpr size_t stashBucket = 2;
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
constexpr size_t bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
constexpr size_t stashMask = (1 << (int)log2(stashBucket)) - 1;
constexpr uint8_t stashHighMask = ~((uint8_t)stashMask);
constexpr uint8_t preNeighborSet = 1 << 7;
constexpr uint8_t nextNeighborSet = 1 << 6;*/
const uint64_t recoverBit = 1UL << 63;
const uint64_t lockBit = 1UL << 62;
/*
#define BUCKET_INDEX(hash) ((hash >> kFingerBits) & bucketMask)
#define GET_COUNT(var) ((var)&countMask)
#define GET_BITMAP(var) (((var) >> 4) & allocMask)
#define ORG_BITMAP(var) ((~((var)&allocMask)) & allocMask)
#define PROBE_BITMAP(var) ((var)&allocMask)*/

#define PARTITION_INDEX(hash) \
  (((hash) >> (64 - partitionShifBits)) & partitionMask)
#define BUCKET_INDEX(hash) (((hash) >> (64 - shiftBits)) & bucketMask)
#define META_HASH(hash) ((uint8_t)((hash) >> (64 - kFingerBits)))
#define GET_COUNT(var) ((var)&countMask)
#define GET_BITMAP(var) (((var) >> 4) & allocMask)
#define ORG_BITMAP(var) ((~((var)&allocMask)) & allocMask)
#define PROBE_BITMAP(var) ((var)&allocMask)

inline bool var_compare(char *str1, char *str2, int len1, int len2) {
  if (len1 != len2) return false;
  return !memcmp(str1, str2, len1);
}