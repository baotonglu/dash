#ifndef CUCKOO
#define CUCKOO
/*
* in this version, I try to aovid the cache misses of unlocking after clflush but has no benefits 
*/
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
#include "util/hash.h"
#include "util/pair.h"
#include "util/persist.h"
#include "util/murmur3.h"
#include <immintrin.h>
#define _INVALID 0 /* we use 0 as the invalid key*/ 
#define SINGLE 1

/* we do avx2 */
#define SIMD 1
#define SIMD_CMP8(src, key)                                                        \
  do                                                                               \
  {                                                                                \
    const __m256i key_data = _mm256_set1_epi8(key);                                \
    __m256i seg_data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src)); \
    __m256i rv_mask = _mm256_cmpeq_epi8(seg_data, key_data);                       \
    mask = _mm256_movemask_epi8(rv_mask);                                          \
  } while (0)

#define SSE_CMP8(src, key)                                                         \
  do                                                                               \
  {                                                                                \
    const __m128i key_data = _mm_set1_epi8(key);                                   \
    __m128i seg_data = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src));    \
    __m128i rv_mask = _mm_cmpeq_epi8(seg_data, key_data);                          \
    mask = _mm_movemask_epi8(rv_mask);                                            \
  } while (0)

#define CHECK_BIT(var, pos) ((((var) & (1<<pos)) > 0) ? (1) : (0))

const uint8_t fingerSet = (uint8_t)(1 << 7);
const uint8_t fingerMask = (uint8_t)((1 << 7) - 1); /*LSB 8 bits as the finger hash*/
const uint8_t allocSet = (uint8_t)(1 << 7);
const uint8_t pointerSet = (uint8_t)(1 << 6);
const uint64_t lockSet = ((uint64_t)1 << 63);
const uint64_t lockMask = ((uint64_t)1 << 63) - 1;
const unsigned fpMask = ((unsigned)1 << 28) - 1; /* the least significant bits are 1s*/
const uint8_t indicatorMask = (uint8_t)((1 << 4) - 1); /*the least 4 bits are 1*/
const uint8_t indicatorSet = (uint8_t)(1 << 7); /*the least 4 bits are 1*/
const uint8_t otherMask = ((uint8_t)255 & (~indicatorMask));
const int overflowSet = 1 << 31;
const int preNeighborSet = 1 << 30;
const int nextNeighborSet = 1 << 29;
const uint8_t lowCountMask = indicatorMask; /*count of kv pair in this bucket*/
const uint8_t highCountMask = otherMask; /*count of overflowed items*/

uint64_t clflushCount;

struct _Pair{
  Key_t key;
  Value_t value;
};

constexpr size_t k_PairSize = sizeof(_Pair); //a k-v _Pair with a bit
constexpr size_t kNumPairPerBucket = 14; /* it is determined by the usage of the fingerprint*/
constexpr size_t kFingerBits = 8;
constexpr size_t kMask = (1 << kFingerBits) -1;
constexpr size_t kNumBucket = 64;
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
constexpr size_t bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
#define BUCKET_INDEX(hash) ((hash >> kFingerBits) & bucketMask)

struct Bucket {
	inline int find_empty_slot(){
		if ((count & lowCountMask) == kNumPairPerBucket)
		{
			return -1;
		}
		auto mask = ~(bitmap & allocMask); //Now the 1 bit should be the empty slot
		return __builtin_ctz(mask);
	}

	inline void set_overflow(){
		bitmap = bitmap | overflowSet;
	}

	/*true indicates overflow, needs extra check in the stash*/
	inline bool test_overflow(){
		return ((bitmap & overflowSet) != 0)?true:false;
	}

	inline void set_indicator(uint8_t meta_hash, Bucket *neighbor){
		int mask = finger_array[14] & indicatorMask;
		mask = ~mask;
		auto index = __builtin_ctz(mask);

		if (index < 4)
		{
			finger_array[15+index] = meta_hash;
			finger_array[14] = ((uint8_t)(1 << index) | finger_array[14]);/*may be optimized*/			
		}else{
			mask = neighbor->finger_array[14] & indicatorMask;
			mask = ~mask;
			index = __builtin_ctz(mask);
			if (index < 4)
			{
				neighbor->finger_array[15+index] = meta_hash;
				neighbor->finger_array[14] = ((uint8_t)(1 << index) | neighbor->finger_array[14]);
				neighbor->finger_array[14] = ((uint8_t)(1 << (4 + index)) | neighbor->finger_array[14]);
			}else{/*set the overflow indicator*/
				set_overflow();
			}
		}
		count = ((((count & highCountMask) >> 4) + 1) << 4) | (count & lowCountMask);
		assert(((count & highCountMask) >> 4) <= kNumPairPerBucket);			
	}


	/*both clear this bucket and its neighbor bucket*/
	inline void unset_indicator(uint8_t meta_hash, Bucket *neighbor, Key_t key){
		/*also needs to ensure that this meta_hash must belongs to other bucket*/
		//printf("i'm cleraing the indicato\n");
		bool clear_success = false;
		int mask = finger_array[14] & indicatorMask;
		for (int i = 0; i < 4; ++i)
		{
			if (CHECK_BIT(mask, i) && (finger_array[15 + i] == meta_hash) && (((1 << (4 + i)) & finger_array[14]) == 0))
			{
				finger_array[14] = finger_array[14] & ((uint8_t)(~(1 << i)));
				clear_success = true;
				break;
			}
		}
		
		if (!clear_success)
		{
			mask = neighbor->finger_array[14] & indicatorMask;
			for (int i = 0; i < 4; ++i)
			{
				if  (CHECK_BIT(mask, i) && (neighbor->finger_array[15 + i] == meta_hash) && (((1 << (4 + i)) & neighbor->finger_array[14]) != 0))
				{
					neighbor->finger_array[14] = neighbor->finger_array[14] & ((uint8_t)(~(1 << i)));
					neighbor->finger_array[14] = neighbor->finger_array[14] & ((uint8_t)(~(1 << (i + 4))));
					break;
				}
			}
		}

		uint8_t num = ((count & highCountMask) >> 4) - 1;
		count = (num << 4) | (count & lowCountMask);
		if (num == 0)
		{
			bitmap = bitmap & (~overflowSet);/*clear the overflow set*/
		}
	}

	int unique_check(uint8_t meta_hash, Key_t key, Bucket* neighbor, Bucket* stash){
		//return check_and_get(meta_hash, key) == NONE ? 0 : -1;
		if ((check_and_get(meta_hash,key) != NONE) || (neighbor->check_and_get(meta_hash, key) != NONE))
		{
			return -1;
		}

		if (finger_array[14] != 0)
		{
			if (test_overflow())
			{
				if (stash->check_and_get(meta_hash, key) != NONE)
				{
					return -1;
				}
			}else{
				int mask = finger_array[14] & indicatorMask;
				auto test_stash = false;
				for (int i = 0; i < 4; ++i)
				{
					if (CHECK_BIT(mask, i) && (finger_array[15+i] == meta_hash))
					{
						test_stash = true;
					}
				}

				if (neighbor->finger_array[14] != 0)
				{
					mask = neighbor->finger_array[14] & indicatorMask;
					for (int i = 0; i < 4; ++i)
						{
							if (CHECK_BIT(mask, i) && (neighbor->finger_array[15+i] == meta_hash))
							{
								test_stash = true;
							}
						}	
				}

				if (test_stash == true)
				{
					if (stash->check_and_get(meta_hash, key) != NONE)
					{
						return -1;
					}
				}
			}
		}

		return 0;
	}

	Value_t check_and_get(uint8_t meta_hash, Key_t key){
		int mask = 0;
  		SSE_CMP8(finger_array, meta_hash);
  		mask = mask & (bitmap & allocMask);
  		/*loop unrolling*/
  		if (mask != 0)
  		{
  			for (int i = 0; i < 12; i+=4)
  			{
  				if (CHECK_BIT(mask, i) && (_[i].key ==  key))
  				{
  					return _[i].value;
  				}

  				if (CHECK_BIT(mask, i+1) && (_[i+1].key ==  key))
  				{
  					return _[i+1].value;
  				}

  				if (CHECK_BIT(mask, i+2) && (_[i+2].key ==  key))
  				{
  					return _[i+2].value;
  				}

  				if (CHECK_BIT(mask, i+3) && (_[i+3].key ==  key))
  				{
  					return _[i+3].value;
  				}
  			}

  			if(CHECK_BIT(mask, 12) && (_[12].key ==  key))
			{
				return _[12].value;
			}

			if(CHECK_BIT(mask, 13) && (_[13].key ==  key))
			{
				return _[13].value;
			}
  		}
  		return NONE;
	}

	inline void set_hash(int index, uint8_t meta_hash) /* Do I needs the atomic instruction????*/
	{
	  finger_array[index] = meta_hash;
	  bitmap = bitmap | (1 << index);
	  //count++;
	  count = ((count & lowCountMask) + 1) | (count & highCountMask);
	  //assert(count <= kNumPairPerBucket);
	  //assert(count > 0);
	}

	inline uint8_t get_hash(int index){
	  return finger_array[index];
	}

	inline void unset_hash(int index){
	  bitmap = bitmap & (~(1 << index));
	  //count--;
	  count = ((count & lowCountMask) - 1) | (count & highCountMask);
	  //assert(count < kNumPairPerBucket);
	  //assert(count >= 0);
	}

	inline void get_lock(){
		auto old_value = version_lock & lockMask;
		auto new_value = version_lock | lockSet;
		while(!CAS(&version_lock, &old_value, new_value)){
			old_value = version_lock & lockMask;
			new_value = version_lock | lockSet;
		}
	}

	inline bool try_get_lock(){
		auto old_value = version_lock & lockMask;
		auto new_value = version_lock | lockSet;
		return CAS(&version_lock, &old_value, new_value);
	}

	inline void release_lock(){
		auto old_value = version_lock;
		auto new_value = ((old_value & lockMask) + 1) & lockMask;

	    while(!CAS(&version_lock, &old_value, new_value)){
	      old_value = version_lock;
	      new_value = ((old_value & lockMask) + 1) & lockMask;
	    }
	}

	/*if the lock is set, return true*/
	inline bool test_lock_set(uint64_t& version){
	    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
	    version = value & lockMask;
	    //printf("the version is %u\n", version);
	    return (value & lockSet) != 0;
	  }

  // test whether the version has change, if change, return true
	inline bool test_lock_version_change(uint64_t old_version){
	    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
	    auto version = value & lockMask;
	    //printf("the old_version is %u\n", version);
	    return ((value & lockSet) != 0) || (version != old_version); 
	  }

	int Insert(Key_t key, Value_t value, uint8_t meta_hash){
		auto slot = find_empty_slot();
		/* this branch can be optimized out*/
		assert(slot < kNumPairPerBucket);
		if (slot == -1)
		{
			printf("cannot find the empty slot, for key %llu\n", key);
			return -1;
		}
		_[slot].value = value;
		_[slot].key = key;
		mfence();
		set_hash(slot, meta_hash);
		//clflush((char*)&bitmap, sizeof(bitmap));
		return 0;
	}

	int Insert_with_noflush(Key_t key, Value_t value, uint8_t meta_hash){
		auto slot = find_empty_slot();
		/* this branch can be optimized out*/
		assert(slot < kNumPairPerBucket);
		if (slot == -1)
		{
			printf("cannot find the empty slot, for key %llu\n", key);
			return -1;
		}
		_[slot].value = value;
		_[slot].key = key;
		set_hash(slot, meta_hash);
		return 0;
	}

	void Insert_displace(Key_t key, Value_t value, uint8_t meta_hash, int slot){
		_[slot].value = value;
		_[slot].key = key;
		mfence();
		set_hash(slot, meta_hash);
		//clflush((char*)&bitmap, sizeof(bitmap));
	}

	void Insert_displace_with_noflush(Key_t key, Value_t value, uint8_t meta_hash, int slot){
		_[slot].value = value;
		_[slot].key = key;
		set_hash(slot, meta_hash);
	}

	/* Find the displacment element in this bucket*/
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
	}

	/*suite of function to set/unset the flush state of the node*/
	inline void setPreNonFlush(){
		//bitmap = bitmap | preNeighborSet;
	}

	inline void unsetPreNonFlush(){
		//bitmap = bitmap & (~preNeighborSet);
	}

	/* If the non-flush is set, return true*/
	inline bool testPreNonFlush(){
		return ((bitmap & preNeighborSet) != 0);
	}

	inline void setNextNonFlush(){
		//bitmap = bitmap | nextNeighborSet;
	}

	inline void unsetNextNonFlush(){
		//bitmap = bitmap & (~nextNeighborSet);
	}

	inline bool testNextNonFlush(){
		return ((bitmap & nextNeighborSet) != 0);
	}

	uint64_t version_lock;
	uint8_t count; /*needs to optimize this count since the comparision of it may consues many time?*/
	uint8_t finger_array[19];/*only use the first 14 bytes, can be accelerated by SSE instruction,0-13 for finger, 14-17 for overflowed, 18 as the bitmap, 19 as the btimap and overflow check indicator*/
	int bitmap; /* this bitmap use 28 bits to store the info about allocation and pointer info, the most significant bit is used to signify whether it needs to check the stash area,
	the most second significant bit is used to signify whether it pre-neighbor is flushed(if not flushed, set as 1), the most third significant bit is to signify whether the next-neighbor is flushed
	if not flushed, set as 1*/
	_Pair _[kNumPairPerBucket];
};

struct Table;

struct Directory{
	Table **_;
	size_t global_depth;
	size_t version;
	Directory(size_t capacity, size_t _version){
		version = _version;
		global_depth = static_cast<size_t>(log2(capacity));
		_ = new Table*[capacity];
	}
};

/* the meta hash-table referenced by the directory*/
struct Table {
  Table(void)
  : local_depth{0}, number{0}, next{nullptr}
  { 
    memset((void*)&bucket[0],0,sizeof(struct Bucket)*(kNumBucket+1)); 
  }

  Table(size_t depth, Table *pp)
  :local_depth{depth}, number{0}, next{pp}
  {
    memset((void*)&bucket[0],0,sizeof(struct Bucket)*(kNumBucket+1)); 
  }
  ~Table(void) {}

  void* operator new(size_t size) {
    void* ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }

  void* operator new[](size_t size) {
    void* ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

  int Insert(Key_t key, Value_t value, size_t key_hash, uint8_t meta_hash, Directory**);
  void Insert4split(Key_t key, Value_t value, size_t key_hash, uint8_t meta_hash);
  Table* Split(size_t);
  //Bucket* Split(Finger *old_f, Finger *new_f, size_t);

  struct Bucket bucket[kNumBucket+1];
  size_t local_depth;
  size_t number;
  Table *next;
};

/* it needs to verify whether this bucket has been deleted...*/
int Table::Insert(Key_t key, Value_t value, size_t key_hash, uint8_t meta_hash, Directory **_dir){
RETRY:
	/*we need to first do the locking and then do the verify*/
	auto y = BUCKET_INDEX(key_hash);
	Bucket* target = bucket + y;
	Bucket* neighbor = bucket + ((y+1) & bucketMask);
	//printf("for key %lld, target bucket is %lld, meta_hash is %d\n", key, BUCKET_INDEX(key_hash), meta_hash);
	target->get_lock();
	neighbor->get_lock();

	auto old_sa = *_dir;
	auto x = (key_hash >> (8*sizeof(key_hash) - old_sa->global_depth));
	if (old_sa->_[x] != this)/* verify process*/
	{
		neighbor->release_lock();
		target->release_lock();
		return -2;
	}

	/*unique check, needs to check 2 hash table*/
	auto ret = target->unique_check(meta_hash, key, neighbor, bucket + kNumBucket);
	if (ret == -1)
	{
		printf("unique check failure!\n");
		return -1;
	}
	
	if (((target->count & lowCountMask) == kNumPairPerBucket) && ((neighbor->count & lowCountMask) == kNumPairPerBucket))
	{
		//return -1
		// Need displacement operation, first do the displacement for the neighbor and then itself 
		//printf("start do the displacement\n");
		Bucket *next_neighbor = bucket + ((y+2) & bucketMask);
		next_neighbor->get_lock();
		int displace_index;
		displace_index = neighbor->Find_displacement((y+1) & bucketMask);
		if (((next_neighbor->count & lowCountMask) != kNumPairPerBucket) && (displace_index != -1))
		{
			//printf("do the displacement in next bucket, the displaced key is %lld, the new key is %lld\n", neighbor->_[displace_index].key, key);
			/*some set/unset may not be necessary*/
			next_neighbor->Insert(neighbor->_[displace_index].key, neighbor->_[displace_index].value, neighbor->finger_array[displace_index]);
			neighbor->setNextNonFlush();
			next_neighbor->release_lock();
			//clflush((char*)&next_neighbor->bitmap, sizeof(int));
			neighbor->unsetNextNonFlush();
			neighbor->unset_hash(displace_index);
			neighbor->Insert_displace(key, value, meta_hash, displace_index);
			target->setNextNonFlush();
			//++number;
			neighbor->release_lock();
			//clflush((char*)&neighbor->bitmap, sizeof(int));
			target->unsetNextNonFlush();
			target->release_lock();
			return 0;
		}
		/* Need to give up current lock to avoid the deadlock, and because of the behaviour of giving up lock, we need further verify whether we still access the right bucket*/
		next_neighbor->release_lock();
		neighbor->release_lock();
		target->release_lock();

		Bucket *prev_neighbor;
		int prev_index;
		/*can be optimized?*/
		if (y == 0)
		{
			prev_neighbor = bucket + kNumBucket - 1;
			prev_index = kNumBucket - 1;
		}else{
			prev_neighbor = bucket + y - 1;
			prev_index = y - 1;
		}
		//re-grab the lock
		prev_neighbor->get_lock();
		target->get_lock();

		/*REVERIFY by re-traversal*/
		old_sa = *_dir;
		x = (key_hash >> (8*sizeof(key_hash) - old_sa->global_depth));
		if (old_sa->_[x] != this)
		{
			target->release_lock();
			prev_neighbor->release_lock();
			return -2;
		}

		displace_index = target->Find_displacement(prev_index);
		if (((prev_neighbor->count & lowCountMask) != kNumPairPerBucket) && (displace_index != -1))
		{
			//printf("do the displacement in previous bucket,the displaced key is %lld, the new key is %lld\n", target->_[displace_index].key, key);
			prev_neighbor->Insert(target->_[displace_index].key, target->_[displace_index].value, target->finger_array[displace_index]);
			target->setPreNonFlush();
			prev_neighbor->release_lock();
			//clflush((char*)&prev_neighbor->bitmap, sizeof(int));
			target->unsetPreNonFlush();
			target->unset_hash(displace_index);
			target->Insert_displace(key, value, meta_hash, displace_index);
			neighbor->setPreNonFlush();
			//++number;
			target->release_lock();
			//clflush((char*)&target->bitmap, sizeof(int));
			neighbor->unsetPreNonFlush();
			return 0;
		}

		/*at last, we can get the lock of stash to insert the kv here*/
		neighbor->get_lock();
		Bucket *stash = bucket + kNumBucket;
		stash->get_lock();
		if ((stash->count & lowCountMask) < kNumPairPerBucket)
		{
			//printf("Insertion in the stash, the key is from x= %lld, y = %lld\n",x, BUCKET_INDEX(key_hash));
			stash->Insert(key, value, meta_hash);
			stash->release_lock();
			//clflush((char*)&stash->bitmap, sizeof(int));
			/*needs to udpate the indicator of current bucket*/
			target->set_indicator(meta_hash, neighbor);
			//++number;
			neighbor->release_lock();
			target->release_lock();
			prev_neighbor->release_lock();
			return 0;
		}

		stash->release_lock();
		neighbor->release_lock();
		target->release_lock();
		prev_neighbor->release_lock();

		return -1;
	}

	if ((target->count & lowCountMask)<=(neighbor->count & lowCountMask))
	{
		target->Insert(key, value, meta_hash);
		neighbor->setPreNonFlush();
		target->release_lock();
		/*do the flush on target*/
		//clflush((char*)&target->bitmap, sizeof(target->bitmap));
		//mfence();
		//flush((char*)&target->bitmap);
		//mfence();
		neighbor->unsetPreNonFlush();
		neighbor->release_lock();
	}else{
		neighbor->Insert(key, value, meta_hash);
		target->setNextNonFlush();
		neighbor->release_lock();
		//clflush((char*)&neighbor->bitmap, sizeof(neighbor->bitmap));/*do the flush on neighbor*/
		//mfence();
		//flush((char*)&neighbor->bitmap);
		//mfence();
		target->unsetNextNonFlush();
		target->release_lock();
	}
	//++number;
	return 0;
}

/*the insert needs to be perfectly balanced, not destory the power of balance*/
void Table::Insert4split(Key_t key, Value_t value, size_t key_hash, uint8_t meta_hash){
	auto y = BUCKET_INDEX(key_hash);
	Bucket* target = bucket + y;
	Bucket* neighbor = bucket + ((y+1) & bucketMask);
	
	auto insert_target = (target->count&lowCountMask)<=(neighbor->count&lowCountMask)?target:neighbor;
	//assert(insert_target->count < kNumPairPerBucket);
	/*some bucket may be overflowed?*/
	if ((insert_target->count&lowCountMask) < kNumPairPerBucket)
	{
		insert_target->_[insert_target->count & lowCountMask].key = key;
		insert_target->_[insert_target->count & lowCountMask].value = value;
		insert_target->set_hash(insert_target->count & lowCountMask, meta_hash);
		//++number;
	}else{
		/*do the displacement or insertion in the stash*/
		
		Bucket *next_neighbor = bucket + ((y+2) & bucketMask);
		int displace_index;
		displace_index = neighbor->Find_displacement((y+1) & bucketMask);
		if (((next_neighbor->count & lowCountMask) != kNumPairPerBucket) && (displace_index != -1))
		{
			//printf("do the displacement in next bucket, the displaced key is %lld, the new key is %lld\n", neighbor->_[displace_index].key, key);
			next_neighbor->Insert_with_noflush(neighbor->_[displace_index].key, neighbor->_[displace_index].value, neighbor->finger_array[displace_index]);
			neighbor->unset_hash(displace_index);
			neighbor->Insert_displace_with_noflush(key, value, meta_hash, displace_index);
			//++number;
			return;
		}
		Bucket *prev_neighbor;
		int prev_index;
		if (y == 0)
		{
			prev_neighbor = bucket + kNumBucket - 1;
			prev_index = kNumBucket - 1;
		}else{
			prev_neighbor = bucket + y - 1;
			prev_index = y - 1;
		}

		displace_index = target->Find_displacement(prev_index);
		if (((prev_neighbor->count & lowCountMask) != kNumPairPerBucket) && (displace_index != -1))
		{
			//printf("do the displacement in previous bucket,the displaced key is %lld, the new key is %lld\n", target->_[displace_index].key, key);
			prev_neighbor->Insert_with_noflush(target->_[displace_index].key, target->_[displace_index].value, target->finger_array[displace_index]);
			target->unset_hash(displace_index);
			target->Insert_displace_with_noflush(key, value, meta_hash, displace_index);
			//++number;
			return;
		}

		Bucket *stash = bucket + kNumBucket;
		if ((stash->count & lowCountMask) < kNumPairPerBucket)
		{
			//printf("Insertion in the stash, the key is from x= %lld, y = %lld\n",x, BUCKET_INDEX(key_hash));
			stash->Insert_with_noflush(key, value, meta_hash);
			target->set_indicator(meta_hash, neighbor);
			//++number;
		}
	}	
}

Table* Table::Split(size_t _key_hash){
	size_t pattern = _key_hash >> (64 - local_depth);
	size_t new_pattern = (pattern << 1) + 1;
	size_t old_pattern = pattern << 1;

	Bucket *curr_bucket;
	for (int i = 1; i < kNumBucket; ++i)
	{
		curr_bucket = bucket + i;
		curr_bucket->get_lock();
	}
	//printf("my pattern is %lld, my load factor is %f\n", pattern, ((double)number)/(kNumBucket*kNumPairPerBucket+kNumPairPerBucket));

	next = new Table(local_depth + 1, next);
	next->bucket->get_lock();/* get the first lock of the new bucket to avoid it is operated(split or merge) by other threads*/
	clflush((char*)&next, sizeof(next));
	size_t key_hash;
	for (int i = 0; i < kNumBucket; ++i)
	{
		curr_bucket = bucket + i;
		auto mask = curr_bucket->bitmap & allocMask;
		for (int j = 0; j < kNumPairPerBucket; ++j)
		{
			if (CHECK_BIT(mask, j))
			{
				key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
				if ((key_hash >> (64 - local_depth - 1)) == new_pattern)
				{
					next->Insert4split(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash, curr_bucket->finger_array[j]); /*this shceme may destory the balanced segment*/
					curr_bucket->unset_hash(j);
					//number--;
				}
			}
		}
	}

	/*split the stash bucket, the stash must be full, right?*/
	curr_bucket = bucket + kNumBucket;
	auto mask = curr_bucket->bitmap & allocMask;
	for (int j = 0; j < kNumPairPerBucket; ++j)
	{
		if (CHECK_BIT(mask, j))
		{
			key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
			if ((key_hash >> (64 - local_depth - 1)) == new_pattern)
			{
				next->Insert4split(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash, curr_bucket->finger_array[j]); /*this shceme may destory the balanced segment*/
				auto bucket_ix = BUCKET_INDEX(key_hash);
				auto org_bucket = bucket + bucket_ix;
				auto neighbor_bucket = bucket + ((bucket_ix + 1) & bucketMask);
				org_bucket->unset_indicator(curr_bucket->finger_array[j], neighbor_bucket, curr_bucket->_[j].key);
				curr_bucket->unset_hash(j);
				//number--;
			}
		}
	}

	clflush((char*)next, sizeof(struct Table));
	return next;
}

class Finger_EH{
	public:
    Finger_EH(void);
    Finger_EH(size_t);
    ~Finger_EH(void);
    int Insert(Key_t key, Value_t value);
    bool InsertOnly(Key_t, Value_t);
    bool Delete(Key_t);
    Value_t Get(Key_t);
    void Directory_Doubling(int x, Table *new_b);
    void Directory_Update(Directory *_sa ,int x, Table *new_b);
    void Lock_Directory();
    void Unlock_Directory();
    int FindAnyway(Key_t key);
    void getNumber(){
		printf("the size of the bucket is %lld\n", sizeof(struct Bucket));
		printf("the size of the Table is %lld\n", sizeof(struct Table));

		size_t _count = 0;
		size_t seg_count = 0;
		Directory *seg = dir;
		Table** dir_entry = seg->_;
		Table *ss;
		auto global_depth = seg->global_depth;
		size_t depth_diff;
		int capacity = pow(2, global_depth);
		for (int i = 0; i < capacity;)
		{
			ss = dir_entry[i];
			depth_diff = global_depth - ss->local_depth;
			_count += ss->number;
			seg_count++;
			i += pow(2, depth_diff);
		}
		printf("the item stored in this hash table is %lld\n", _count);
		printf("the segment number in this hash table is %lld\n", seg_count);
		printf("the space Utilization is %f\n", (double)_count/(seg_count*kNumPairPerBucket*(kNumBucket+1)));
	}

    inline bool Acquire(void) {
      int unlocked = 0;
      return CAS(&lock, &unlocked, 1);
    }
  
    inline bool Release(void) {
      int locked = 1;
      return CAS(&lock, &locked, 0);
    }

    inline int Test_Directory_Lock_Set(void){
      return __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    }

    Directory *dir;
    int lock;
};

Finger_EH::Finger_EH(void){

}

Finger_EH::Finger_EH(size_t initCap){
	dir = new Directory(initCap, 0);
	lock = 0;

	dir->_[initCap - 1] = new Table(dir->global_depth, nullptr);
	/* Initilize the Directory*/
	for (int i = initCap - 2; i >= 0; --i)
	{
		dir->_[i] = new Table(dir->global_depth, dir->_[i+1]);
	}
}

Finger_EH::~Finger_EH(void){
	//TO-DO
}

void Finger_EH::Lock_Directory(){
	while(!Acquire()){
		asm("nop");
	}
}

void Finger_EH::Unlock_Directory(){
  while (!Release()) {
      asm("nop");
  }
}

void Finger_EH::Directory_Doubling(int x, Table *new_b){
  Table** d = dir->_;
  auto global_depth = dir->global_depth;
  printf("Directory_Doubling towards %lld\n", global_depth+1);

  auto new_sa = new Directory(2*pow(2, global_depth), dir->version + 1);
  auto dd = new_sa->_;

  auto capacity = pow(2, global_depth);
  for (unsigned i = 0; i < capacity; ++i) {
      dd[2*i] = d[i];
      dd[2*i+1] = d[i];
  }
  dd[2*x+1] = new_b;

  clflush((char*)new_sa, sizeof(struct Directory));
  dir = new_sa;
  clflush((char*)&dir, sizeof(dir));
  /*need to delete the old directory...*/
  printf("Done!!Directory_Doubling towards %lld\n", dir->global_depth);
}

void Finger_EH::Directory_Update(Directory *_sa, int x, Table *new_b){
  //printf("directory update for %d\n", x);
  Table** dir_entry = _sa->_;
  auto global_depth = _sa->global_depth;
  unsigned depth_diff = global_depth - new_b->local_depth;
  if (depth_diff == 0) {
    if (x%2 == 0) {
      dir_entry[x+1] = new_b;
    } else {
      dir_entry[x] = new_b;
    }
  } else {
    int chunk_size = pow(2, global_depth - (new_b->local_depth - 1));
    x = x - (x % chunk_size);
    for (unsigned i = 0; i < chunk_size/2; ++i) {
      dir_entry[x+chunk_size/2+i] = new_b;
    }
  }
  //printf("Done!directory update for %d\n", x);
}

int Finger_EH::Insert(Key_t key, Value_t value) {
   auto key_hash = h(&key, sizeof(key));
   auto meta_hash = ((uint8_t)(key_hash & kMask));//the last 8 bits
RETRY:
  auto old_sa = dir;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Table* target = dir_entry[x];

  auto ret = target->Insert(key, value, key_hash, meta_hash, &dir);
  //printf("insert key %lld, x = %d, y = %d, meta_hash = %d\n", key, x, BUCKET_INDEX(key_hash), meta_hash);
  if (ret == -1)
  {
  	/*Insertion failure ,split process needs to get all the lock in this segment? it has much overfead? I do not know...need to test whether this has much overhead!*/
  	if (!target->bucket->try_get_lock())
	{
		goto RETRY;
	}

	/*verify procedure*/
	auto old_sa = dir;
	auto x = (key_hash >> (8*sizeof(key_hash) - old_sa->global_depth));
	if (old_sa->_[x] != target)/* verify process*/
	{
		target->bucket->release_lock();
		goto RETRY;
	}

  	auto new_b = target->Split(key_hash);/* also needs the verify..., and we use try lock for this rather than the spin lock*/
  	target->local_depth += 1;
  	clflush((char*)&target->local_depth, sizeof(target->local_depth));
  	/* update directory*/
REINSERT:    
    // the following three statements may be unnecessary...
    old_sa = dir;
    dir_entry = old_sa->_;
    x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
    //assert(target == dir_entry[x].bucket_p);
    if (target->local_depth-1 < old_sa->global_depth)
    {
      Directory_Update(old_sa,x,new_b);
      /*after the update, I need to recheck*/
      if (Test_Directory_Lock_Set() || old_sa->version != dir->version)
      {
        goto REINSERT;
      }
    }else{
      Lock_Directory();
      if (old_sa->version != dir->version)
      {
        Unlock_Directory();
        goto REINSERT;
      }
      Directory_Doubling(x,new_b);
      Unlock_Directory();
    }
    /*release the lock for the target bucket and the new bucket*/
    Bucket *curr_bucket;
    for (int i = 0; i < kNumBucket; ++i)
    {
      curr_bucket = target->bucket + i;
      curr_bucket->release_lock();
    }
    curr_bucket = new_b->bucket;
    curr_bucket->release_lock();
    goto RETRY;
  }else if(ret == -2){
  	goto RETRY;
  }

  return 0;
}

Value_t Finger_EH::Get(Key_t key){
  auto key_hash = h(&key, sizeof(key));
  auto meta_hash = ((uint8_t)(key_hash & kMask));//the last 8 bits
RETRY:
  auto old_sa = dir;
  auto x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  auto y = BUCKET_INDEX(key_hash);
  auto dir_entry = old_sa->_;
  Table* target = dir_entry[x];

  /* after accessing the right bucket, still need the verification process to verify that I access the right bucket?* or no need to verify?
	Need the re-traversal for the verification since assuming it go to sleep after it get the pointer to this segment, then the segment has been splitted, the element it wanna search has been re-hashed to other places
	then it is wrong, so after get the version nnumber of the lock, it needs to re-traversal to do the verifciation!
   */
  uint64_t old_version;
  Bucket *target_bucket = target->bucket + y;
  Bucket *neighbor_bucket = target->bucket + ((y+1) & bucketMask);
  //printf("Get key %lld, x = %d, y = %d, meta_hash = %d\n", key, x, BUCKET_INDEX(key_hash), meta_hash);

  if (target_bucket->test_lock_set(old_version) || neighbor_bucket->testPreNonFlush())
  {
  	goto RETRY;
  }

  /*verification procedure*/
  old_sa = dir;
  x = (key_hash >> (8*sizeof(key_hash)-old_sa->global_depth));
  if (old_sa->_[x] != target)
  {
  	goto RETRY;
  }

  //printf("for key %lld, the target table is %lld, target bucket is %lld, meta_hash is %d\n", key, x, y, meta_hash);

  auto ret = target_bucket->check_and_get(meta_hash, key);
  if (ret != NONE && !(target_bucket->test_lock_version_change(old_version)))
  {
  	return ret;
  }
  
  uint64_t _version;
  if (neighbor_bucket->test_lock_set(_version) || target_bucket->testNextNonFlush())
  {
  	goto RETRY;
  }
  /*no need for verification procedure, we use the version number of target_bucket to test whether the bucket has ben spliteted*/
  ret = neighbor_bucket->check_and_get(meta_hash, key);
  if (neighbor_bucket->test_lock_version_change(_version) ||  target_bucket->test_lock_version_change(old_version))
  {
  	goto RETRY;
  }

  if (ret != NONE)
  {
  	return ret;
  }

  /*next step is to search in the stash*/
  if ((target_bucket->count & highCountMask) != 0)
	{
		if (target_bucket->test_overflow())
		{
			Bucket *stash = target->bucket + kNumBucket;
			ret = stash->check_and_get(meta_hash, key);
		}else{
			auto test_stash = false;
			/*search in the original bucket*/
			if (target_bucket->finger_array[14] != 0)
			{
				int mask = target_bucket->finger_array[14] & indicatorMask;
				for (int i = 0; i < 4; ++i)
				{
					if (CHECK_BIT(mask, i) && (target_bucket->finger_array[15+i] == meta_hash))
					{
						test_stash = true;
						break;
					}
				}
			}
			
			/*search in the neighbor bucket*/
			if (neighbor_bucket->finger_array[14] != 0)
			{
				int mask = neighbor_bucket->finger_array[14] & indicatorMask;
				for (int i = 0; i < 4; ++i)
				{
					if (CHECK_BIT(mask, i) && (neighbor_bucket->finger_array[15+i] == meta_hash))
					{
						test_stash = true;
						break;
					}
				}	
			}

			if (test_stash == true)
			{
				Bucket *stash = target->bucket + kNumBucket;
				ret = stash->check_and_get(meta_hash, key);
			}
		}
	}

	/*at last, verify whether the lock in primary bucket change*/
	if(target_bucket->test_lock_version_change(old_version)){
		goto RETRY;
	}
	/*need to verify whether the version_lock in the primary bucket change*/
	//printf("the x = %lld, the y = %lld, the meta_hash is %d\n", x, y, meta_hash);
	return ret;
}

/*DEBUG FUNCTION: search the position of the key in this table and print correspongdign informantion in this table, to test whether it is correct*/
int Finger_EH::FindAnyway(Key_t key){
	auto key_hash = h(&key, sizeof(key));
  	auto meta_hash = ((uint8_t)(key_hash & kMask));
  	auto x = (key_hash >> (8*sizeof(key_hash)-dir->global_depth));

	size_t _count = 0;
	size_t seg_count = 0;
	Directory *seg = dir;
	Table** dir_entry = seg->_;
	Table *ss;
	auto global_depth = seg->global_depth;
	size_t depth_diff;
	int capacity = pow(2, global_depth);
	for (int i = 0; i < capacity;)
	{
		ss = dir_entry[i];
		/*search in this table*/
		Bucket *curr_bucket;
		for (int i = 0; i < kNumBucket; ++i)
		{
			curr_bucket = ss->bucket + i;
			auto ret = curr_bucket->check_and_get(meta_hash, key);
			if (ret != NONE)
			{
				printf("successfully find in the normal bucket\n");
				return 0;
			}
		}

		/*search in the stash*/
		curr_bucket = ss->bucket + kNumBucket;
		auto ret = curr_bucket->check_and_get(meta_hash, key);
		if (ret != NONE)
		{
			printf("successfully find in the stash bucket\n");
			/*Print the image in the original bucket and neighbor bucket*/
			auto bucket_ix = BUCKET_INDEX(key_hash);
			auto org_bucket = ss->bucket + bucket_ix;
			auto neighbor_bucket = ss->bucket + ((bucket_ix + 1) & bucketMask);
			printf("the segment number is %d, the bucket_ix is %d\n", x, bucket_ix);

			printf("the image of org_bucket\n");
			int mask = org_bucket->finger_array[14] & indicatorMask;
			for (int i = 0; i < 4; ++i)
			{
				printf("the hash is %d, the pos bit is %d, the alloc bit is %d\n", org_bucket->finger_array[15 + i], (org_bucket->finger_array[14] >> (4 + i)) & 1, (org_bucket->finger_array[14] >> i) & 1);
			}

			printf("the image of the neighbor bucket\n");
			mask = neighbor_bucket->finger_array[14] & indicatorMask;
			for (int i = 0; i < 4; ++i)
			{
				printf("the hash is %d, the pos bit is %d, the alloc bit is %d\n", neighbor_bucket->finger_array[15 + i], (neighbor_bucket->finger_array[14] >> (4 + i)) & 1, (neighbor_bucket->finger_array[14] >> i) & 1);
			}

			if (org_bucket->test_overflow())
			{
				printf("the org bucket has overflowed\n");
			}
			return 0;
		}

		depth_diff = global_depth - ss->local_depth;
		_count += ss->number;
		seg_count++;
		i += pow(2, depth_diff);
	}
	return -1;
}

#endif
