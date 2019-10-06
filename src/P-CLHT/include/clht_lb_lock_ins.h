/*   
 *   File: clht_lb_lock_ins.h
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: lock-based cache-line hash table with no resizing. Only
 *    inserts use real locks. Removals work with a CAS. 
 *   clht_lb_lock_ins.h is part of ASCYLIB
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

#ifndef _CLHT_LOCK_INS_H_
#define _CLHT_LOCK_INS_H_

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "atomic_ops.h"
#include "utils.h"

#include "ssmem.h"

#define true 1
#define false 0

/* #define DEBUG */

#define CLHT_READ_ONLY_FAIL   1
#define CLHT_HELP_RESIZE      1
#define CLHT_PERC_EXPANSIONS  0.05
#define CLHT_MAX_EXPANSIONS   2
#define CLHT_PERC_FULL_DOUBLE 80	   /* % */
#define CLHT_RATIO_DOUBLE     2		  
#define CLHT_PERC_FULL_HALVE  5		   /* % */
#define CLHT_RATIO_HALVE      8		  
#define CLHT_MIN_CLHT_SIZE      8
#define CLHT_DO_CHECK_STATUS  0
#define CLHT_DO_GC            0
#define CLHT_STATUS_INVOK     500000
#define CLHT_STATUS_INVOK_IN  500000
#if defined(RTM)	       /* only for processors that have RTM */
#define CLHT_USE_RTM          1
#else
#define CLHT_USE_RTM          0
#endif

#if CLHT_DO_CHECK_STATUS == 1
#  define CLHT_CHECK_STATUS(h)				\
  if (unlikely((--check_ht_status_steps) == 0))		\
    {							\
      ht_status(h, 0, 0);				\
      check_ht_status_steps = CLHT_STATUS_INVOK;	\
    }

#else 
#  define CLHT_CHECK_STATUS(h)
#endif

#if CLHT_DO_GC == 1
#  define CLHT_GC_HT_VERSION_USED(ht) clht_gc_thread_version(ht)
#else
#  define CLHT_GC_HT_VERSION_USED(ht)
#endif

#if defined(DEBUG)
#  define DPP(x)	x++				
#else
#  define DPP(x)
#endif

#define CACHE_LINE_SIZE    64
#define ENTRIES_PER_BUCKET 3

#define KEY_BLCK -1

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
inline int is_power_of_two(unsigned int x);

typedef uintptr_t clht_addr_t;
typedef volatile uintptr_t clht_val_t;

#if defined(__tile__)
typedef volatile uint32_t clht_lock_t;
#else
typedef volatile uint8_t clht_lock_t;
#endif

typedef struct ALIGNED(CACHE_LINE_SIZE) bucket_s
{
  clht_lock_t lock;
  clht_addr_t key[ENTRIES_PER_BUCKET];
  clht_val_t  val[ENTRIES_PER_BUCKET];
  volatile struct bucket_s* next;
} bucket_t;


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
      volatile uint32_t num_expands_threshold;
      volatile int32_t is_helper;
      volatile int32_t helper_done;
      size_t version_min;
    };
    uint8_t padding[2*CACHE_LINE_SIZE];
  };
} clht_hashtable_t;

typedef struct ALIGNED(CACHE_LINE_SIZE) ht_ts
{
  union
  {
    struct
    {
      size_t version;
      clht_hashtable_t* versionp;
      int id;
      volatile struct ht_ts* next;
    };
    uint8_t padding[CACHE_LINE_SIZE];
  };
} ht_ts_t;


inline uint64_t __ac_Jenkins_hash_64(uint64_t key);

/* Hash a key for a particular hashtable. */
uint64_t clht_hash(clht_hashtable_t* hashtable, clht_addr_t key );

static inline void
_mm_pause_rep(uint64_t w)
{
  while (w--)
    {
      _mm_pause();
    }
}

#if defined(XEON) | defined(COREi7) | defined(__tile__)
#  define TAS_RLS_MFENCE() _mm_sfence();
#else
#  define TAS_RLS_MFENCE()
#endif

#define LOCK_FREE   0
#define LOCK_UPDATE 1
#define LOCK_RESIZE 2

#if CLHT_USE_RTM == 1		/* USE RTM */
#  define LOCK_ACQ(lock, ht)			\
  lock_acq_rtm_chk_resize(lock, ht)
#  define LOCK_RLS(lock)			\
  if (likely(*(lock) == LOCK_FREE))		\
    {						\
      _xend();					\
      DPP(put_num_failed_on_new);		\
    }						\
  else						\
    {						\
      TAS_RLS_MFENCE();				\
     *lock = LOCK_FREE;				\
      DPP(put_num_failed_expand);		\
    }
#else  /* NO RTM */
#  define LOCK_ACQ(lock, ht)			\
  lock_acq_chk_resize(lock, ht)

#  define LOCK_RLS(lock)			\
  TAS_RLS_MFENCE();				\
 *lock = 0;	  

#endif	/* RTM */

#define LOCK_ACQ_RES(lock)			\
  lock_acq_resize(lock)

#define TRYLOCK_ACQ(lock)			\
  TAS_U8(lock)

#define TRYLOCK_RLS(lock)			\
  lock = LOCK_FREE

void ht_resize_help(clht_hashtable_t* h);

#if defined(DEBUG)
extern __thread uint32_t put_num_restarts;
#endif

static inline int
lock_acq_chk_resize(clht_lock_t* lock, clht_hashtable_t* h)
{
  char once = 1;
  clht_lock_t l;
  while ((l = CAS_U8(lock, LOCK_FREE, LOCK_UPDATE)) == LOCK_UPDATE)
    {
      if (once)
      	{
      	  DPP(put_num_restarts);
      	  once = 0;
      	}
      _mm_pause();
    }
  return 1;
}

static inline int
lock_acq_resize(clht_lock_t* lock)
{
  clht_lock_t l;
  while ((l = CAS_U8(lock, LOCK_FREE, LOCK_RESIZE)) == LOCK_UPDATE)
    {
      _mm_pause();
    }

  if (l == LOCK_RESIZE)
    {
      return 0;
    }

  return 1;
}


/* ******************************************************************************** */
#if CLHT_USE_RTM == 1  /* use RTM */
/* ******************************************************************************** */

#include <immintrin.h>		/*  */

static inline int
lock_acq_rtm_chk_resize(clht_lock_t* lock, clht_hashtable_t* h)
{

  int rtm_retries = 1;
  do 
    {
      /* while (unlikely(*lock == LOCK_UPDATE)) */
      /* 	{ */
      /* 	  _mm_pause(); */
      /* 	} */

      if (likely(_xbegin() == _XBEGIN_STARTED))
	{
	  clht_lock_t lv = *lock;
	  if (likely(lv == LOCK_FREE))
	    {
	      return 1;
	    }
	  else if (lv == LOCK_RESIZE)
	    {
	      _xend();
#  if CLHT_HELP_RESIZE == 1
	      ht_resize_help(h);
#  endif

	      while (h->table_new == NULL)
		{
		  _mm_mfence();
		}

	      return 0;
	    }

	  DPP(put_num_restarts);
	  _xabort(0xff);
	}
    } while (rtm_retries-- > 0);

  return lock_acq_chk_resize(lock, h);
}
#endif	/* RTM */


/* ******************************************************************************** */
/* intefance */
/* ******************************************************************************** */

/* Create a new hashtable. */
clht_hashtable_t* clht_hashtable_create(uint64_t num_buckets);
clht_t* clht_create(uint64_t num_buckets);

/* Insert a key-value pair into a hashtable. */
int clht_put(clht_t* hashtable, clht_addr_t key, clht_val_t val);

/* Retrieve a key-value pair from a hashtable. */
clht_val_t clht_get(clht_hashtable_t* hashtable, clht_addr_t key);

/* Remove a key-value pair from a hashtable. */
clht_val_t clht_remove(clht_t* hashtable, clht_addr_t key);

size_t clht_size(clht_hashtable_t* hashtable);
size_t clht_size_mem(clht_hashtable_t* hashtable);
size_t clht_size_mem_garbage(clht_hashtable_t* hashtable);

void clht_gc_thread_init(clht_t* hashtable, int id);
inline void clht_gc_thread_version(clht_hashtable_t* h);
inline int clht_gc_get_id();
int clht_gc_collect(clht_t* h);
int clht_gc_collect_all(clht_t* h);
int clht_gc_free(clht_hashtable_t* hashtable);
void clht_gc_destroy(clht_t* hashtable);

void clht_print(clht_hashtable_t* hashtable);
size_t ht_status(clht_t* hashtable, int resize_increase, int just_print);

bucket_t* clht_bucket_create();
int ht_resize_pes(clht_t* hashtable, int is_increase, int by);

const char* clht_type_desc();

#endif /* _CLHT_LOCK_INS_H_ */

