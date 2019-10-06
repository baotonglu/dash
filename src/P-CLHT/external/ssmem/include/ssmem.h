/*   
 *   File: ssmem.h
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: ssmem interface and structures
 *   ssmem.h is part of ASCYLIB
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

/* December 10, 2013 */
#ifndef _SSMEM_H_
#define _SSMEM_H_

#include <stdio.h>
#include <stdint.h>

/* **************************************************************************************** */
/* parameters */
/* **************************************************************************************** */

#define SSMEM_TRANSPARENT_HUGE_PAGES 0 /* Use or not Linux transparent huge pages */
#define SSMEM_ZERO_MEMORY            0 /* Initialize allocated memory to 0 or not */
#define SSMEM_GC_FREE_SET_SIZE 507 /* mem objects to free before doing a GC pass */
#define SSMEM_GC_RLSE_SET_SIZE 3   /* num of released object before doing a GC pass */
#define SSMEM_DEFAULT_MEM_SIZE (32 * 1024 * 1024L) /* memory-chunk size that each threads
						    gives to the allocators */
#define SSMEM_MEM_SIZE_DOUBLE  1 /* if the allocator is out of memory, should it allocate
				  a 2x larger chunk than before? (in order to stop asking
				 for memory again and again */
#define SSMEM_MEM_SIZE_MAX     (4 * 1024 * 1024 * 1024LL) /* absolute max chunk size 
							   (e.g., if doubling is 1) */

/* increase the thread-local timestamp of activity on each ssmem_alloc() and/or ssmem_free() 
   call. If enabled (>0), after some memory is alloced and/or freed, the thread should not 
   access ANY ssmem-protected memory that was read (the reference were taken) before the
   current alloc or free invocation. If disabled (0), the program should employ manual 
   SSMEM_SAFE_TO_RECLAIM() calls to indicate when the thread does not hold any ssmem-allocated
   memory references. */

#define SSMEM_TS_INCR_ON_NONE   0
#define SSMEM_TS_INCR_ON_BOTH   1
#define SSMEM_TS_INCR_ON_ALLOC  2
#define SSMEM_TS_INCR_ON_FREE   3

#define SSMEM_TS_INCR_ON        SSMEM_TS_INCR_ON_FREE
/* **************************************************************************************** */
/* help definitions */
/* **************************************************************************************** */
#define ALIGNED(N) __attribute__ ((aligned (N)))
#define CACHE_LINE_SIZE 64

/* **************************************************************************************** */
/* data structures used by ssmem */
/* **************************************************************************************** */

/* an ssmem allocator */
typedef struct ALIGNED(CACHE_LINE_SIZE) ssmem_allocator
{
  union
  {
    struct
    {
      void* mem;		/* the actual memory the allocator uses */
      size_t mem_curr;		/* pointer to the next addrr to be allocated */
      size_t mem_size;		/* size of mem chunk */
      size_t tot_size;		/* total memory that the allocator uses */
      size_t fs_size;		/* size (in objects) of free_sets */
      struct ssmem_list* mem_chunks; /* list of mem chunks (used to free the mem) */

      struct ssmem_ts* ts;	/* timestamp object associated with the allocator */

      struct ssmem_free_set* free_set_list; /* list of free_set. A free set holds freed mem 
					     that has not yet been reclaimed */
      size_t free_set_num;	/* number of sets in the free_set_list */
      struct ssmem_free_set* collected_set_list; /* list of collected_set. A collected set
						  contains mem that has been reclaimed */
      size_t collected_set_num;	/* number of sets in the collected_set_list */
      struct ssmem_free_set* available_set_list; /* list of set structs that are not used
						  and can be used as free sets */
      size_t released_num;	/* number of released memory objects */
      struct ssmem_released* released_mem_list; /* list of release memory objects */
    };
    uint8_t padding[2 * CACHE_LINE_SIZE];
  };
} ssmem_allocator_t;

/* a timestamp used by a thread */
typedef struct ALIGNED(CACHE_LINE_SIZE) ssmem_ts
{
  union
  {
    struct
    {
      size_t version;
      size_t id;
      struct ssmem_ts* next;
    };
  };
  uint8_t padding[CACHE_LINE_SIZE];
} ssmem_ts_t;

/* 
 * a timestamped free_set. It holds:  
 *  1. the collection of timestamps at the point when the free_set gets full
 *  2. the array of freed pointers to be used by ssmem_free()
 *  3. a set_next pointer in order to be able to create linked lists of
 *   free_sets
 */
typedef struct ALIGNED(CACHE_LINE_SIZE) ssmem_free_set
{
  size_t* ts_set;		/* set of timestamps for GC */
  size_t size;
  long int curr;		
  struct ssmem_free_set* set_next;
  uintptr_t* set;
} ssmem_free_set_t;


/* 
 * a timestamped node of released memory. The memory will be returned to the OS
 * (free(node->mem)) when the current timestamp is greater than the one of the node
 */
typedef struct ssmem_released
{
  size_t* ts_set;
  void* mem;
  struct ssmem_released* next;
} ssmem_released_t;

/*
 * a generic list that keeps track of actual memory that has been allocated
 * (using malloc / memalign) and the different allocators that the list is using
 */
typedef struct ssmem_list
{
  void* obj;
  struct ssmem_list* next;
} ssmem_list_t;

/* **************************************************************************************** */
/* ssmem interface */
/* **************************************************************************************** */

/* initialize an allocator with the default number of objects */
void ssmem_alloc_init(ssmem_allocator_t* a, size_t size, int id);
/* initialize an allocator and give the number of objects in free_sets */
void ssmem_alloc_init_fs_size(ssmem_allocator_t* a, size_t size, size_t free_set_size, int id);
/* explicitely subscribe to the list of threads in order to used timestamps for GC */
void ssmem_gc_thread_init(ssmem_allocator_t* a, int id);
/* terminate the system (all allocators) and free all memory */
void ssmem_term();
/* terminate the allocator a and free all its memory
 * This function should NOT be used if the memory allocated by this allocator
 * might have been freed (and is still in use) by other allocators */
void ssmem_alloc_term(ssmem_allocator_t* a);

/* allocate some memory using allocator a */
void* ssmem_alloc(ssmem_allocator_t* a, size_t size);
/* free some memory using allocator a */
void ssmem_free(ssmem_allocator_t* a, void* obj);

/* release some memory to the OS using allocator a */
void ssmem_release(ssmem_allocator_t* a, void* obj);

/* increment the thread-local activity counter. Invoking this function suggests that
 no memory references to ssmem-allocated memory are held by the current thread beyond
this point. */
void ssmem_ts_next();
#define SSMEM_SAFE_TO_RECLAIM() ssmem_ts_next()


/* debug/help functions */
void ssmem_ts_list_print();
size_t* ssmem_ts_set_collect();
void ssmem_ts_set_print(size_t* set);

void ssmem_free_list_print(ssmem_allocator_t* a);
void ssmem_collected_list_print(ssmem_allocator_t* a);
void ssmem_available_list_print(ssmem_allocator_t* a);
void ssmem_all_list_print(ssmem_allocator_t* a, int id);


/* **************************************************************************************** */
/* platform-specific definitions */
/* **************************************************************************************** */

#if defined(__x86_64__)
#  define CAS_U64(a,b,c) __sync_val_compare_and_swap(a,b,c)
#  define FAI_U32(a) __sync_fetch_and_add(a,1)
#endif

#if defined(__sparc__)
#  include <atomic.h>
#  define CAS_U64(a,b,c) atomic_cas_64(a,b,c)
#  define FAI_U32(a) (atomic_inc_32_nv(a)-1)
#endif

#if defined(__tile__)
#  include <arch/atomic.h>
#  include <arch/cycle.h>
#  define CAS_U64(a,b,c) arch_atomic_val_compare_and_exchange(a,b,c)
#  define FAI_U32(a) arch_atomic_increment(a)
#endif


#endif /* _SSMEM_H_ */

