/*   
 *   File: clht.h
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: CLHT common interface
 *   clht.h is part of ASCYLIB
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *	      	      Distributed Programming Lab (LPD), EPFL
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef _CLHT_H_
#define _CLHT_H_

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#define true 1
#define false 0

/* #define DEBUG */

#if defined(DEBUG)
#  define DPP(x)	x++				
#else
#  define DPP(x)
#endif

#define CACHE_LINE_SIZE    64

#define MAP_INVLD 0
#define MAP_VALID 1
#define MAP_INSRT 2

#define KEY_BUCKT 3
#define ENTRIES_PER_BUCKET KEY_BUCKT

#define CLHT_DO_GC                  1
#define CLHT_PERC_FULL_HALVE        2
#define CLHT_PERC_FULL_DOUBLE       15
#define CLHT_OCCUP_AFTER_RES        40
#define CLHT_INC_EMERGENCY          2
#define CLHT_NO_EMPTY_SLOT_TRIES    16
#define CLHT_GC_HT_VERSION_USED(ht) clht_gc_thread_version(ht)
#define LOAD_FACTOR                 2

#ifndef ALIGNED
#  if __GNUC__ && !SCC
#    define ALIGNED(N) __attribute__ ((aligned (N)))
#  else
#    define ALIGNED(N)
#  endif
#endif

#define likely(x)       __builtin_expect((x), 1)
#define unlikely(x)     __builtin_expect((x), 0)

#if defined(__sparc__)
#  define PREFETCHW(x) 
#  define PREFETCH(x) 
#  define PREFETCHNTA(x) 
#  define PREFETCHT0(x) 
#  define PREFETCHT1(x) 
#  define PREFETCHT2(x) 

#  define PAUSE    asm volatile("rd    %%ccr, %%g0\n\t" \
				::: "memory")
#  define _mm_pause() PAUSE
#  define _mm_mfence() __asm__ __volatile__("membar #LoadLoad | #LoadStore | #StoreLoad | #StoreStore");
#  define _mm_lfence() __asm__ __volatile__("membar #LoadLoad | #LoadStore");
#  define _mm_sfence() __asm__ __volatile__("membar #StoreLoad | #StoreStore");


#elif defined(__tile__)
#  define _mm_lfence() arch_atomic_read_barrier()
#  define _mm_sfence() arch_atomic_write_barrier()
#  define _mm_mfence() arch_atomic_full_barrier()
#  define _mm_pause() cycle_relax()
#endif

#define CAS_U64_BOOL(a, b, c) (CAS_U64(a, b, c) == b)
static inline int
is_power_of_two (unsigned int x) 
{
  return ((x != 0) && !(x & (x - 1)));
}


typedef uintptr_t clht_addr_t;
typedef volatile uintptr_t clht_val_t;
typedef uint64_t clht_snapshot_all_t;

typedef union
{
  volatile uint64_t snapshot;
  struct
  {
#if KEY_BUCKT == 4
    uint32_t version;
#elif KEY_BUCKT == 6
    uint16_t version;
#else
    uint32_t version;
#endif
    uint8_t map[KEY_BUCKT];
  };
} clht_snapshot_t;

#if __GNUC__ > 4 && __GNUC_MINOR__ > 4
_Static_assert (sizeof(clht_snapshot_t) == 8, "sizeof(clht_snapshot_t) == 8");
#endif

typedef volatile struct ALIGNED(CACHE_LINE_SIZE) bucket_s
{
  union
  {
    volatile uint64_t snapshot;
    struct
    {
#if KEY_BUCKT == 4
      uint32_t version;
#elif KEY_BUCKT == 6
      uint16_t version;
#else
      uint32_t version;
/* #  error "KEY_BUCKT should be either 4 or 6" */
#endif
      uint8_t map[KEY_BUCKT];
    };
  };
  clht_addr_t key[KEY_BUCKT];
  clht_val_t  val[KEY_BUCKT];
  volatile struct bucket_s* padding;
} bucket_t;

#if __GNUC__ > 4 && __GNUC_MINOR__ > 4
_Static_assert (sizeof(bucket_t) % 64 == 0, "sizeof(bucket_t) == 64");
#endif

#if defined(__tile__)
typedef volatile uint32_t clht_lock_t;
#else
typedef volatile uint8_t clht_lock_t;
#endif
/* typedef volatile uint64_t clht_lock_t; */
#define CLHT_LOCK_FREE 0
#define CLHT_LOCK_ACQR 1

#define CLHT_CHECK_RESIZE(w)				\
  while (unlikely(w->resize_lock == CLHT_LOCK_ACQR))	\
    {							\
      _mm_pause();					\
      CLHT_GC_HT_VERSION_USED(w->ht);			\
    }

#define CLHT_LOCK_RESIZE(w)						\
  (CAS_U8(&w->resize_lock, CLHT_LOCK_FREE, CLHT_LOCK_ACQR) == CLHT_LOCK_FREE)

#define CLHT_RLS_RESIZE(w)			\
  w->resize_lock = CLHT_LOCK_FREE

#define TRYLOCK_ACQ(lock)			\
  TAS_U8(lock)

#define TRYLOCK_RLS(lock)			\
  lock = CLHT_LOCK_FREE


typedef struct ALIGNED(CACHE_LINE_SIZE) clht
{
  union
  {
    struct
    {
      struct clht_hashtable_s* ht;
      uint8_t next_cache_line[CACHE_LINE_SIZE - (sizeof(void*))];
      struct clht_hashtable_s* ht_oldest;
      struct ht_ts* version_list;
      size_t version_min;
      volatile clht_lock_t resize_lock;
      volatile clht_lock_t gc_lock;
      volatile clht_lock_t status_lock;
    };
    uint8_t padding[2 * CACHE_LINE_SIZE];
  };
} clht_t;

typedef struct ALIGNED(CACHE_LINE_SIZE) clht_hashtable_s
{
  union
  {
    struct
    {
      size_t num_buckets;
      bucket_t* table;
      size_t hash;
      size_t version;
      uint8_t next_cache_line[CACHE_LINE_SIZE - (3 * sizeof(size_t)) - (sizeof(void*))];
      struct clht_hashtable_s* table_tmp;
      struct clht_hashtable_s* table_prev;
      struct clht_hashtable_s* table_new;
      volatile uint32_t num_expands;
      union
      {
	volatile uint32_t num_expands_threshold;
	uint32_t num_buckets_prev;
      };
      volatile int32_t is_helper;
      volatile int32_t helper_done;
      size_t version_min;
    };
    uint8_t padding[2*CACHE_LINE_SIZE];
  };
} clht_hashtable_t;

extern uint64_t __ac_Jenkins_hash_64(uint64_t key);

/* Hash a key for a particular hashtable. */
uint64_t clht_hash(clht_hashtable_t* hashtable, clht_addr_t key);


static inline int
snap_get_empty_index(uint64_t snap)
{
  clht_snapshot_t s = { .snapshot = snap };
  int i;
  for (i = 0; i < KEY_BUCKT; i++)
    {
      if (s.map[i] == MAP_INVLD)
	{
	  return i;
	}
    }
  return -1;
}

static inline int
keys_get_empty_index(clht_addr_t* keys)
{
  int i;
  for (i = 0; i < KEY_BUCKT; i++)
    {
      if (keys[i] == 0)
	{
	  return i;
	}
    }
  return -1;
}

static inline int
buck_get_empty_index(bucket_t* b, uint64_t snap)
{
  clht_snapshot_t s = { .snapshot = snap };

  int i;
  for (i = 0; i < KEY_BUCKT; i++)
    {
      if (b->key[i] == 0 && s.map[i] != MAP_INSRT)
	{
	  return i;
	}
    }
  return -1;
}


static inline int
vals_get_empty_index(clht_val_t* vals, clht_snapshot_all_t snap)
{
  clht_snapshot_t s = { .snapshot = snap };

  int i;
  for (i = 0; i < KEY_BUCKT; i++)
    {
      if (vals[i] == 0 && s.map[i] != MAP_INSRT)
	{
	  return i;
	}
    }
  return -1;
}


static inline uint64_t
snap_set_map(uint64_t s, int index, int val)
{
  clht_snapshot_t s1 = { .snapshot = s };
  s1.map[index] = val;
  return s1.snapshot;
}

static inline uint64_t
snap_set_map_and_inc_version(uint64_t s, int index, int val)
{
  clht_snapshot_t s1 = { .snapshot =  s};
  s1.map[index] = val;
  s1.version++;
  return s1.snapshot;
}

static inline void
_mm_pause_rep(uint64_t w)
{
  while (w--)
    {
      _mm_pause();
    }
}



/* ******************************************************************************** */
/* inteface */
/* ******************************************************************************** */

/* Create a new hashtable. */
clht_t* clht_create(uint64_t num_buckets);
/* initializes the necessary per-thread structures for the hash table */
void clht_gc_thread_init(clht_t* hashtable, int id);


/* Insert a key-value pair into a hashtable. */
int clht_put(clht_t* hashtable, clht_addr_t key, clht_val_t val);
/* Retrieve a key-value pair from a hashtable. */
clht_val_t clht_get(clht_hashtable_t* hashtable, clht_addr_t key);
/* Remove a key-value pair from a hashtable. */
clht_val_t clht_remove(clht_t* hashtable, clht_addr_t key);

/* returns the size of the hash table */
size_t clht_size(clht_hashtable_t* hashtable);

/* frees the memory used by the hashtable */
void clht_gc_destroy(clht_t* hashtable);

/* prints the contents of the hash table */
void clht_print(clht_hashtable_t* hashtable);

/* string description of the type of the hash table
 For example, CLHT-LB-RESIZE */
const char* clht_type_desc();

void clht_lock_initialization(clht_t *h);

/* internal */
extern void clht_gc_thread_version(clht_hashtable_t* h);
extern void clht_gc_thread_version_max();
extern int clht_gc_get_id();

#ifdef __cplusplus
}
#endif

#endif /* _CLHT_H_ */

