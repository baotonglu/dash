/*   
 *   File: sspfd.h
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: sspfd interface, structures, and helper functions
 *   sspfd.h is part of sspfd
 *
 * The MIT License (MIT)
 *
 * Copyright (C) 2013  Vasileios Trigonakis
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

#ifndef _SSPFD_H_
#define _SSPFD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <inttypes.h>
#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#define SSPFD_DO_TIMINGS 1

#if SSPFD_DO_TIMINGS != 1	/* empty macros when not benchmarkings */
/* 
 * initialize `num_stores` stores, with `num_entries` number of measurement entries each. 
 * Use `id` as the thread id for printing purposes. Each store should be used to measure a single event. 
 */
#  define SSPFDINIT(num_stores, num_entries, id) 

/* 
 * terminate (free) the initialized stores.
 */
#  define SSPFDTERM()

/* 
 * start measurement (i.e., take start timestamp at this point) for store `store`.
 */
#  define SSPFDI(store) 

/* 
 * only if sspfd_get_id() == id, start measurement (i.e., take start timestamp at this point) 
 * for store `store`,
 */
#  define SSPFDI_ID(store, id) 

/* 
 * start measurement (i.e., take start timestamp at this point) for any store. Need to do the calc
 * with SSPFDO_G() then.
 */
#  define SSPFDI_G() 

/* 
 * only if sspfd_get_id() == id, start measurement (i.e., take start timestamp at this point) 
 * for any store. Need to do the calc with SSPFDO_ID_G() then.
 */
#  define SSPFDI_ID_G(id) 

/* 
 * stop measuring (i.e., take stop timestamp at this point) for store `store` and store the duration
 * since `SSPFDI(store)` in entry `entry`.
 */
#  define SSPFDO(store, entry) 

/* 
 * if sspfd_get_id() == id, stop measuring (i.e., take stop timestamp at this point) for store `store` 
 * and store the duration since `SSPFDI(store)` in entry `entry`.
 */
#  define SSPFDO_ID(store, entry) 

/* 
 * stop measuring (i.e., take stop timestamp at this point) for SSPFDI_G and store the duration
 * since `SSPFDI(store)` in entry `entry`.
 */
#  define SSPFDO_G(store, entry) 

/* 
 * if sspfd_get_id() == id, stop measuring (i.e., take stop timestamp at this point) for SSPFDI_ID_G
 * and store the duration since `SSPFDI(store)` in entry `entry`.
 */
#  define SSPFDO_ID_G(store, entry, id) 

/* 
 * generate statistics and print them for the first `num_vals` values of store `store`.
 */
#  define SSPFDP(store, num_vals) 

/* 
 * generate statistics and print them for the first `num_vals` values of store `store`. 
 * Additionally, print the first `num_print` measurements of this store.
 */
#  define SSPFDPN(store, num_vals, num_print)

/* 
 * generate statistics and print them for the first `num_vals` values of store `store`. 
 * Additionally, print the first `num_print` measurements of this store, in a comma-
 * seperated format.
 */
#  define SSPFDPN_COMMA(store, num_vals, num_print)

/* 
 * prefetch entry `entry` for store `store`, so that the overheads are minimized (only 
 * necessary if the application has accessed a huge amount of data before using `sspfd`. 
 */
#  define SSPFDPREFTCH(store, entry) 

/* 
 * generate statistics for the measurements in store `store` for the first `num_ops` values. 
 * Store the results in `statsp` pointer to a `sspfd_stats_t` structure. 
 */
#  define SSPFDSTATS(store, num_ops, statsp)

/* 
 * print the statistics in `statsp` pointer to a `sspfd_stats_t` structure.  
 */
#  define SSPFDPRINT(statsp)

/* 
 * print the first `num_print` measurements from store `store`.
 */
#  define SSPFDPRINTV(num_store, num_print)
/* 
 * print the first `num_print` measurements from store `store`, comma separated.
 */
#  define SSPFDPRINTV_COMMA(num_store, num_vals, num_print)

#endif 

#define SSPFD_PRINT(args...) printf("[%02lu] ", sspfd_get_id()); printf(args); printf("\n"); fflush(stdout)

typedef uint64_t ticks;

#if !defined(_GETTICKS_H_) && !defined(_H_GETTICKS_)
#  if defined(__i386__)
static inline ticks 
getticks(void) 
{
  ticks ret;

  __asm__ __volatile__("rdtsc" : "=A" (ret));
  return ret;
}
#  elif defined(__x86_64__)
static inline ticks
getticks(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#  elif defined(__sparc__)
static inline ticks
getticks()
{
  ticks ret;
  __asm__ __volatile__ ("rd %%tick, %0" : "=r" (ret) : "0" (ret)); 
  return ret;
}
#  elif defined(__tile__)
#    include <arch/cycle.h>
static inline ticks getticks()
{
  return get_cycle_count();
}
#  endif
#endif	/* _H_GETTICKS_ */

#if !defined(PREFETCHW)
#  if defined(__x86_64__) | defined(__i386__)
#    define PREFETCHW(x) asm volatile("prefetchw %0" :: "m" (*(unsigned long *)x)) /* write */
#  elif defined(__sparc__)
#    define PREFETCHW(x) __builtin_prefetch((const void*) x, 1, 3)
#  elif defined(__tile__)
#    define PREFETCHW(x) tmc_mem_prefetch (x, 64)
#  else
#    warning "You need to define PREFETCHW(x) for your architecture"
#  endif
#endif

typedef struct sspfd_stats
{
  uint64_t num_vals;
  double avg;
  double avg_10p;
  double avg_25p;
  double avg_50p;
  double avg_75p;
  double avg_rst;
  double abs_dev_10p;
  double abs_dev_25p;
  double abs_dev_50p;
  double abs_dev_75p;
  double abs_dev_rst;
  double abs_dev;
  double abs_dev_perc;
  double std_dev_10p;
  double std_dev_25p;
  double std_dev_50p;
  double std_dev_75p;
  double std_dev_rst;
  double std_dev;
  double std_dev_perc;
  double min_val;
  uint64_t min_val_idx;
  double max_val;
  uint64_t max_val_idx;
  uint32_t num_dev_10p;
  uint32_t num_dev_25p;
  uint32_t num_dev_50p;
  uint32_t num_dev_75p;
  uint32_t num_dev_rst;
} sspfd_stats_t;

#define SSPFD_PRINT_MAX 200

extern __thread volatile size_t sspfd_num_stores;
extern __thread volatile ticks** sspfd_store;
extern __thread volatile ticks* _sspfd_s;
extern __thread ticks _sspfd_s_global;
extern __thread volatile ticks sspfd_correction;

#if SSPFD_DO_TIMINGS == 1

#  define SSPFDINIT(num_stores, num_entries, id) sspfd_store_init(num_stores, num_entries, id)

#  define SSPFDTERM() sspfd_store_term()

#  define SSPFDI(store)				\
  {						\
  asm volatile ("");				\
  _sspfd_s[store] = getticks();

#  define SSPFDI_ID(store, id)			\
  {						\
  asm volatile ("");				\
  if (sspfd_get_id() == id)			\
    {						\
      _sspfd_s[store] = getticks();		\
    }

#  define SSPFDI_G()				\
  asm volatile("" ::: "memory");		\
  _sspfd_s_global = getticks();			\
  asm volatile("" ::: "memory");		

#  define SSPFDI_ID_G(id)			\
  asm volatile ("");				\
  if (sspfd_get_id() == id)			\
    {						\
      asm volatile("" ::: "memory");		\
      _sspfd_s_global = getticks();		\
      asm volatile("" ::: "memory");		\
    }

#  define SSPFDO(store, entry)						\
  asm volatile ("");							\
  sspfd_store[store][entry] =  getticks() - _sspfd_s[store] - sspfd_correction; \
  }

#  define SSPFDO_ID(store, entry, id)					\
  asm volatile("" ::: "memory");					\
  if (sspfd_get_id() == id)						\
    {									\
      asm volatile("" ::: "memory");					\
      sspfd_store[store][entry] =  getticks() - _sspfd_s[store] - sspfd_correction; \
      asm volatile("" ::: "memory");					\
    }									\
  }

#  define SSPFDO_G(store, entry)					\
  asm volatile("" ::: "memory");					\
  sspfd_store[store][entry] =  getticks() - _sspfd_s_global - sspfd_correction;	\
  asm volatile("" ::: "memory");					

#  define SSPFDO_ID_G(store, entry, id)					\
  asm volatile("" ::: "memory");					\
  if (sspfd_get_id() == id)						\
    {									\
      asm volatile("" ::: "memory");					\
      sspfd_store[store][entry] =  getticks() - _sspfd_s_global - sspfd_correction; \
      asm volatile("" ::: "memory");					\
    }								       


#  define SSPFDP(store, num_vals)		\
  {						\
    sspfd_stats_t ad;				\
    sspfd_get_stats(store, num_vals, &ad);	\
    sspfd_print_stats(&ad);			\
  }

#  define SSPFDPN(store, num_vals, num_print)				\
  {									\
    uint32_t _i;							\
    uint32_t p = num_print;						\
    if (p > num_vals) { p = num_vals; }					\
    for (_i = 0; _i < p; _i++)						\
      {									\
	printf("[%3d: %4ld] ", _i, (long int) sspfd_store[store][_i]);	\
      }									\
    sspfd_stats_t ad;							\
    sspfd_get_stats(store, num_vals, &ad);				\
    sspfd_print_stats(&ad);						\
  }

#  define SSPFDPN_COMMA(store, num_vals, num_print)		\
  {								\
    uint32_t _i;						\
    uint32_t p = num_print;					\
    if (p > num_vals) { p = num_vals; }				\
    for (_i = 0; _i < p; _i++)					\
      {								\
	    printf("%ld,", (long int) sspfd_store[store][_i]);	\
      }								\
    printf("\n");						\
    sspfd_stats_t ad;						\
    sspfd_get_stats(store, num_vals, &ad);			\
    sspfd_print_stats(&ad);					\
  }

#  define SSPFDPREFTCH(store, entry)		\
  SSPFDI(store);				\
  SSPFDO(store, entry);

#  define SSPFDSTATS(num_stores, num_ops, statsp) sspfd_get_stats(num_stores, num_ops, &stats)
#  define SSPFDPRINT(statsp) sspfd_print_stats(statsp)

#  define SSPFDPRINTV(num_store, num_print)				\
  {									\
    uint32_t _i;							\
    uint32_t p = num_print;						\
    if (p > num_vals) { p = num_vals; }					\
    for (_i = 0; _i < p; _i++)						\
      {									\
	printf("[%3d: %4ld] ", _i, (long int) sspfd_store[store][_i]);	\
      }									\
  }

#  define SSPFDPRINTV_COMMA(store, num_vals, num_print)		\
  {								\
    uint32_t _i;						\
    uint32_t p = num_print;					\
    if (p > num_vals) { p = num_vals; }				\
    for (_i = 0; _i < p; _i++)					\
      {								\
	printf("%ld,", (long int) sspfd_store[store][_i]);	\
      }								\
    printf("\n");						\
  }

#endif /* !SSPFD_DO_TIMINGS */


void sspfd_set_id(size_t id);
size_t sspfd_get_id();

void sspfd_store_init(size_t num_stores, size_t num_entries, size_t id);
void sspfd_store_term();
void sspfd_get_stats(const size_t store, const size_t num_vals, sspfd_stats_t* sspfd_stats);
void sspfd_print_stats(const sspfd_stats_t* sspfd_stats);

#ifdef __cplusplus
}
#endif

#endif	/* _SSPFD_H_ */
