/*   
 *   File: latency.h
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: 
 *   latency.h is part of ASCYLIB
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

#ifndef _LATENCY_H_
#define _LATENCY_H_

#define EXEC_IN_DEC_ID_ORDER(id, nthr)		\
  { int __i;					\
  for (__i = nthr - 1; __i >= 0; __i--)		\
    {						\
  if (id == __i)				\
    {

#define EXEC_IN_DEC_ID_ORDER_END(barrier)	\
  }						\
    barrier_cross(barrier);			\
    }}

#if RETRY_STATS == 1
#  define RETRY_STATS_VARS						\
  __thread size_t __parse_try, __update_try, __cleanup_try, __lock_try, __lock_queue, __lock_try_once
#  define RETRY_STATS_VARS_GLOBAL					\
  size_t __parse_try_global, __update_try_global, __cleanup_try_global, __lock_try_global, __lock_queue_global

extern RETRY_STATS_VARS;
extern RETRY_STATS_VARS_GLOBAL;

#  define RETRY_STATS_ZERO()			\
  __parse_try = 0;				\
  __update_try = 0;				\
  __cleanup_try = 0;				\
  __lock_try = 0;				\
  __lock_queue = 0;				\
  __lock_try_once = 1;

#  define PARSE_TRY()        __parse_try++
#  define UPDATE_TRY()       __update_try++
#  define CLEANUP_TRY()      __cleanup_try++
#  define LOCK_TRY()         __lock_try++
#  define LOCK_TRY_ONCE()			\
  if (__lock_try_once)				\
    {						\
      LOCK_TRY();				\
    }
#  define LOCK_QUEUE(q)      __lock_queue += (q)
#  define LOCK_QUEUE_ONCE(q)			\
  if (__lock_try_once)				\
    {						\
      __lock_try_once = 0;			\
      __lock_queue += (q);			\
    }
#  define LOCK_TRY_ONCE_CLEAR()    __lock_try_once = 1
#  define RETRY_STATS_PRINT(thr, put, rem, upd_suc)   retry_stats_print(thr, put, rem, (upd_suc))
#  define RETRY_STATS_SHARE()			\
  __parse_try_global +=  __parse_try;		\
  __update_try_global += __update_try;		\
  __cleanup_try_global += __cleanup_try;	\
  __lock_try_global += __lock_try;		\
  __lock_queue_global += __lock_queue

static inline void 
retry_stats_print(size_t thr, size_t put, size_t rem, size_t upd_suc)
{
  double ratio_all = 100.0 * (((double) __parse_try_global / thr) - 1);
  printf("#parse_all:    %-10zu %f %%\n", __parse_try_global, ratio_all);
  size_t updates = put + rem;
  ratio_all = 100.0 * (((double) __update_try_global / updates) - 1);
  printf("#update_all:   %-10zu %f %%\n", __update_try_global, ratio_all);
  ratio_all = 100.0 * (__cleanup_try_global / (thr - (double) __cleanup_try_global));
  double ratio_upd = 100.0 * (__cleanup_try_global / (updates - (double) __cleanup_try_global));
  printf("#cleanup_all:  %-10zu %f %% %f %%\n", __cleanup_try_global, ratio_upd, ratio_all);
  ratio_all = (double) __lock_queue_global / __lock_try_global;
  if (__lock_try_global == 0)
    {
      ratio_all = 0;
    }
  double ratio_to_succ_upd = (double) __lock_try_global / upd_suc;
  printf("#lock_all:     %-10zu %f   %f\n", __lock_try_global, ratio_all, ratio_to_succ_upd);
}

#else  /* RETRY_STATS == 0 */
#  define RETRY_STATS_VARS
#  define RETRY_STATS_VARS_GLOBAL
#  define RETRY_STATS_ZERO()

#  define PARSE_TRY()     
#  define UPDATE_TRY()    
#  define CLEANUP_TRY()   
#  define RETRY_STATS_PRINT(thr, put, rem, upd_suc)
#  define RETRY_STATS_SHARE()
#  define LOCK_TRY()
#  define LOCK_QUEUE(q)
#  define LOCK_TRY_ONCE_CLEAR()
#endif	/* RETRY_STATS */

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* latency */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

#if defined(USE_SSPFD)
#  define PFD_TYPE 1
#  define SSPFD_DO_TIMINGS 1
#else
#  define PFD_TYPE 0
#  define SSPFD_NUM_ENTRIES 0
#  undef SSPFD_DO_TIMINGS
#  define SSPFD_DO_TIMINGS 0
#endif
#include "sspfd.h"

#ifdef __tile__
#  include <arch/atomic.h>
#  define LFENCE arch_atomic_read_barrier()
#elif defined(__sparc__)
#  define LFENCE  asm volatile("membar #LoadLoad"); 
#else 
#  define LFENCE asm volatile ("lfence")
#endif

#if !defined(COMPUTE_LATENCY)
#  define START_TS(s)
#  define END_TS(s, i)
#  define END_TS_ELSE(s, i, inc)
#  define ADD_DUR(tar)
#  define ADD_DUR_FAIL(tar)
#  define PF_INIT(s, e, id)
#  define PARSE_START_TS(s)
#  define PARSE_END_TS(s, i)
#  define PARSE_END_INC(i)
#elif PFD_TYPE == 0
#  define PARSE_START_TS(s)
#  define PARSE_END_TS(s, i)
#  define PARSE_END_INC(i)
#  define START_TS(s)				\
    asm volatile ("");				\
    start_acq = getticks();			\
    LFENCE;
#  define END_TS(s, i)				\
    asm volatile ("");				\
    end_acq = getticks();			\
    asm volatile ("");
#  define END_TS_ELSE(s, i, inc)		\
  else						\
    {						\
      END_TS(s, i);				\
      ADD_DUR(inc);				\
    }
#  define ADD_DUR(tar) tar += (end_acq - start_acq - correction)
#  define ADD_DUR_FAIL(tar)					\
  else								\
    {								\
      ADD_DUR(tar);						\
    }
#  define PF_INIT(s, e, id)
#else  /* PFD_TYPE == 1 */

#  define PF_NUM_STORES 6
#  define SSPFD_NUM_ENTRIES  (pf_vals_num + 1)
#  define PF_INIT(s, e, id) SSPFDINIT(PF_NUM_STORES, e, id)

#  if LATENCY_PARSING == 0
#    define PARSE_START_TS(s)
#    define PARSE_END_TS(s, i)
#    define PARSE_END_INC(i)
#    if LATENCY_ALL_CORES == 0
#      define START_TS(s)      SSPFDI_ID_G(0); LFENCE;
#      define END_TS(s, i)     SSPFDO_ID_G(s, (i) & pf_vals_num, 0)
#      define END_TS_ELSE(s, i, inc)     else { SSPFDO_ID_G(s, (i) & pf_vals_num, 0); }

#      define ADD_DUR(tar) 
#      define ADD_DUR_FAIL(tar)
#    else	 /* LATENCY_ALL_CORES == 1 */
#      define START_TS(s)      SSPFDI_G(); LFENCE;
#      define END_TS(s, i)     SSPFDO_G(s, (i) & pf_vals_num)
#      define END_TS_ELSE(s, i, inc)     else { SSPFDO_G(s, (i) & pf_vals_num); }

#      define ADD_DUR(tar) 
#      define ADD_DUR_FAIL(tar)
#    endif /* LATENCY_ALL_CORES  */
#  else /* LATENCY_PARSING == 1*/
extern size_t pf_vals_num;
#    define START_TS(s)
#    define END_TS(s, i)
#    define END_TS_ELSE(s, i, inc)
#    define ADD_DUR(tar)
#    define ADD_DUR_FAIL(tar)

#    define PARSE_START_TS(s)     SSPFDI_G(); LFENCE;
#    define PARSE_END_TS(s, i)    SSPFDO_G(s, (i) & pf_vals_num)
#    define PARSE_END_INC(i)      i++

#  endif	 /* LATENCY_PARSING */
#endif

static inline void
print_latency_stats(int ID, size_t num_entries, size_t num_entries_print)
{
#if (PFD_TYPE == 1) && defined(COMPUTE_LATENCY)
  if (ID == 0)
    {
#  if LATENCY_PARSING == 1
      printf("get ------------------------------------------------------------------------\n");
      printf("#latency_get_parse: ");
      SSPFDPN_COMMA(0, num_entries, num_entries_print);
      printf("put ------------------------------------------------------------------------\n");
      printf("#latency_put_parse: ");
      SSPFDPN_COMMA(1, num_entries, num_entries_print);
      printf("rem ------------------------------------------------------------------------\n");
      printf("#latency_rem_parse: ");
      SSPFDPN_COMMA(2, num_entries, num_entries_print);
#  else  /* LATENCY_PARSING == 0*/
      printf("get ------------------------------------------------------------------------\n");
      printf("#latency_get_suc: ");
      SSPFDPN_COMMA(0, num_entries, num_entries_print);
      printf("#latency_get_fal: ");
      SSPFDPN_COMMA(3, num_entries, num_entries_print);
      printf("put ------------------------------------------------------------------------\n");
      printf("#latency_put_suc: ");
      SSPFDPN_COMMA(1, num_entries, num_entries_print);
      printf("#latency_put_fal: ");
      SSPFDPN_COMMA(4, num_entries, num_entries_print);
      printf("rem ------------------------------------------------------------------------\n");
      printf("#latency_rem_suc: ");
      SSPFDPN_COMMA(2, num_entries, num_entries_print);
      printf("#latency_rem_fal: ");
      SSPFDPN_COMMA(5, num_entries, num_entries_print);
#  endif	/* LATENCY_PARSING */
    }
#  if LATENCY_ALL_CORES == 1
  else
    {
#    if LATENCY_PARSING == 1
      printf("get ------------------------------------------------------------------------\n");
      printf("#latency_get_parse: ");
      SSPFDPN_COMMA(0, num_entries, num_entries_print);
      printf("put ------------------------------------------------------------------------\n");
      printf("#latency_put_parse: ");
      SSPFDPN_COMMA(1, num_entries, num_entries_print);
      printf("rem ------------------------------------------------------------------------\n");
      printf("#latency_rem_parse: ");
      SSPFDPN_COMMA(2, num_entries, num_entries_print);
#    else  /* LATENCY_PARSING == 0*/
      printf("get ------------------------------------------------------------------------\n");
      printf("#latency_get_suc: ");
      SSPFDPN_COMMA(0, num_entries, num_entries_print);
      printf("#latency_get_fal: ");
      SSPFDPN_COMMA(3, num_entries, num_entries_print);
      printf("put ------------------------------------------------------------------------\n");
      printf("#latency_put_suc: ");
      SSPFDPN_COMMA(1, num_entries, num_entries_print);
      printf("#latency_put_fal: ");
      SSPFDPN_COMMA(4, num_entries, num_entries_print);
      printf("rem ------------------------------------------------------------------------\n");
      printf("#latency_rem_suc: ");
      SSPFDPN_COMMA(2, num_entries, num_entries_print);
      printf("#latency_rem_fal: ");
      SSPFDPN_COMMA(5, num_entries, num_entries_print);
#    endif	/* LATENCY_PARSING */
}
#  endif
#endif
}

#endif /* _LATENCY_H_ */
