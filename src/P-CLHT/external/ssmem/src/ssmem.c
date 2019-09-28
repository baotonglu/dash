/*   
 *   File: ssmem.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: a simple object-based memory allocator with epoch-based 
 *                garbage collection
 *   ssmem.c is part of ASCYLIB
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

#include "ssmem.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>

ssmem_ts_t* ssmem_ts_list = NULL;
volatile uint32_t ssmem_ts_list_len = 0;
__thread volatile ssmem_ts_t* ssmem_ts_local = NULL;
__thread size_t ssmem_num_allocators = 0;
__thread ssmem_list_t* ssmem_allocator_list = NULL;

inline int 
ssmem_get_id()
{
  if (ssmem_ts_local != NULL)
    {
      return ssmem_ts_local->id;
    }
  return -1;
}

static ssmem_list_t* ssmem_list_node_new(void* mem, ssmem_list_t* next);

/* 
 * explicitely subscribe to the list of threads in order to used timestamps for GC
 */
void
ssmem_gc_thread_init(ssmem_allocator_t* a, int id)
{  
  a->ts = (ssmem_ts_t*) ssmem_ts_local;
  if (a->ts == NULL)
    {
      a->ts = (ssmem_ts_t*) memalign(CACHE_LINE_SIZE, sizeof(ssmem_ts_t));
      assert (a->ts != NULL);
      ssmem_ts_local = a->ts;

      a->ts->id = id;
      a->ts->version = 0;

      do
	{
	  a->ts->next = ssmem_ts_list;
	}
      while (CAS_U64((volatile uint64_t*) &ssmem_ts_list, 
		     (uint64_t) a->ts->next, (uint64_t) a->ts) 
	     != (uint64_t) a->ts->next);
  
       __attribute__ ((unused)) uint32_t null = FAI_U32(&ssmem_ts_list_len);
    }
}

ssmem_free_set_t* ssmem_free_set_new(size_t size, ssmem_free_set_t* next);


/* 
 * initialize allocator a with a custom free_set_size
 * If the thread is not subscribed to the list of timestamps (used for GC),
 * additionally subscribe the thread to the list
 */
void
ssmem_alloc_init_fs_size(ssmem_allocator_t* a, size_t size, size_t free_set_size, int id)
{
  ssmem_num_allocators++;
  ssmem_allocator_list = ssmem_list_node_new((void*) a, ssmem_allocator_list);

  if (id == 0)
    {
//      printf("[ALLOC] initializing allocator with fs size: %zu objects\n", free_set_size);
    }

#if SSMEM_TRANSPARENT_HUGE_PAGES
  int ret = posix_memalign(&a->mem, CACHE_LINE_SIZE, size);
  assert(ret == 0);
#else
  a->mem = (void*) memalign(CACHE_LINE_SIZE, size);
#endif
  assert(a->mem != NULL);
#if SSMEM_ZERO_MEMORY == 1
  memset(a->mem, 0, size);
#endif

  a->mem_curr = 0;
  a->mem_size = size;
  a->tot_size = size;
  a->fs_size = free_set_size;

  a->mem_chunks = ssmem_list_node_new(a->mem, NULL);

  ssmem_gc_thread_init(a, id);

  a->free_set_list = ssmem_free_set_new(a->fs_size, NULL);
  a->free_set_num = 1;

  a->collected_set_list = NULL;
  a->collected_set_num = 0;

  a->available_set_list = NULL;

  a->released_mem_list = NULL;
  a->released_num = 0;
}

/* 
 * initialize allocator a with the default SSMEM_GC_FREE_SET_SIZE
 * If the thread is not subscribed to the list of timestamps (used for GC),
 * additionally subscribe the thread to the list
 */
void
ssmem_alloc_init(ssmem_allocator_t* a, size_t size, int id)
{
  return ssmem_alloc_init_fs_size(a, size, SSMEM_GC_FREE_SET_SIZE, id);
}



/* 
 * 
 */
static ssmem_list_t*
ssmem_list_node_new(void* mem, ssmem_list_t* next)
{
  ssmem_list_t* mc;
  mc = (ssmem_list_t*) malloc(sizeof(ssmem_list_t));
  assert(mc != NULL);
  mc->obj = mem;
  mc->next = next;

  return mc;
}

/* 
 *
 */
inline ssmem_released_t*
ssmem_released_node_new(void* mem, ssmem_released_t* next)
{
  ssmem_released_t* rel;
  rel = (ssmem_released_t*) malloc(sizeof(ssmem_released_t) + (ssmem_ts_list_len * sizeof(size_t)));
  assert(rel != NULL);
  rel->mem = mem;
  rel->next = next;
  rel->ts_set = (size_t*) (rel + 1);

  return rel;
}

/* 
 * 
 */
ssmem_free_set_t*
ssmem_free_set_new(size_t size, ssmem_free_set_t* next)
{
  /* allocate both the ssmem_free_set_t and the free_set with one call */
  ssmem_free_set_t* fs = (ssmem_free_set_t*) memalign(CACHE_LINE_SIZE, sizeof(ssmem_free_set_t) + (size * sizeof(uintptr_t)));
  assert(fs != NULL);

  fs->size = size;
  fs->curr = 0;
  
  fs->set = (uintptr_t*) (((uintptr_t) fs) + sizeof(ssmem_free_set_t));
  fs->ts_set = NULL;	      /* will get a ts when it becomes full */
  fs->set_next = next;

  return fs;
}


/* 
 * 
 */
ssmem_free_set_t*
ssmem_free_set_get_avail(ssmem_allocator_t* a, size_t size, ssmem_free_set_t* next)
{
  ssmem_free_set_t* fs;
  if (a->available_set_list != NULL)
    {
      fs = a->available_set_list;
      a->available_set_list = fs->set_next;

      fs->curr = 0;
      fs->set_next = next;

      /* printf("[ALLOC] got free_set from available_set : %p\n", fs); */
    }
  else
    {
      fs = ssmem_free_set_new(size, next);
    }

  return fs;
}


/* 
 * 
 */
static void
ssmem_free_set_free(ssmem_free_set_t* set)
{
  free(set->ts_set);
  free(set);
}

/* 
 * 
 */
static inline void
ssmem_free_set_make_avail(ssmem_allocator_t* a, ssmem_free_set_t* set)
{
  /* printf("[ALLOC] added to avail_set : %p\n", set); */
  set->curr = 0;
  set->set_next = a->available_set_list;
  a->available_set_list = set;
}


/* 
 * terminated allocator a and free its memory
 */
void
ssmem_alloc_term(ssmem_allocator_t* a)
{
  /* printf("[ALLOC] term() : ~ total mem used: %zu bytes = %zu KB = %zu MB\n", */
  /* 	 a->tot_size, a->tot_size / 1024, a->tot_size / (1024 * 1024)); */
  ssmem_list_t* mcur = a->mem_chunks;
  do
    {
      ssmem_list_t* mnxt = mcur->next;
      free(mcur->obj);
      free(mcur);
      mcur = mnxt;
    }
  while (mcur != NULL);

  ssmem_list_t* prv = ssmem_allocator_list;
  ssmem_list_t* cur = ssmem_allocator_list;
  while (cur != NULL && (uintptr_t) cur->obj != (uintptr_t) a)
    {
      prv = cur;
      cur = cur->next;
    }

  if (cur == NULL)
    {
      printf("[ALLOC] ssmem_alloc_term: could not find %p in the ssmem_allocator_list\n", a);
    }
  else if (cur == prv)
    {
      ssmem_allocator_list = cur->next;
    }
  else
    {
      prv->next = cur->next;
    }

  if (--ssmem_num_allocators == 0)
    {
      free(a->ts);
    }


  /* printf("[ALLOC] free(free_set)\n"); fflush(stdout); */
  /* freeing free sets */
  ssmem_free_set_t* fs = a->free_set_list;
  while (fs != NULL)
    {
      ssmem_free_set_t* nxt = fs->set_next;
      ssmem_free_set_free(fs);
      fs = nxt;
    }

  /* printf("[ALLOC] free(collected_set)\n"); fflush(stdout); */
  /* freeing collected sets */
  fs = a->collected_set_list;
  while (fs != NULL)
    {
      ssmem_free_set_t* nxt = fs->set_next;
      ssmem_free_set_free(fs);
      fs = nxt;
    }

  /* printf("[ALLOC] free(available_set)\n"); fflush(stdout); */
  /* freeing available sets */
  fs = a->available_set_list;
  while (fs != NULL)
    {
      ssmem_free_set_t* nxt = fs->set_next;
      ssmem_free_set_free(fs);
      fs = nxt;
    }

  /* freeing the relased memory */
  ssmem_released_t* rel = a->released_mem_list;
  while (rel != NULL)
    {
      ssmem_released_t* next = rel->next;
      free(rel->mem);
      free(rel);
      rel = next;
    }

 }

/* 
 * terminate all allocators
 */
void
ssmem_term()
{
  while (ssmem_allocator_list != NULL)
    {
      ssmem_alloc_term((ssmem_allocator_t*) ssmem_allocator_list->obj);
    }
}

/* 
 * 
 */
inline void 
ssmem_ts_next()
{
  ssmem_ts_local->version++;
}

/* 
 * 
 */
size_t*
ssmem_ts_set_collect(size_t* ts_set)
{
  if (ts_set == NULL)
    {
      ts_set = (size_t*) malloc(ssmem_ts_list_len * sizeof(size_t));
      assert(ts_set != NULL);
    }

  ssmem_ts_t* cur = ssmem_ts_list;
  while (cur != NULL)
    {
      ts_set[cur->id] = cur->version;
      cur = cur->next;
    }

  return ts_set;
}

/* 
 * 
 */
void 
ssmem_ts_set_print(size_t* set)
{
  printf("[ALLOC] set: [");
  int i;
  for (i = 0; i < ssmem_ts_list_len; i++)
    {
      printf("%zu | ", set[i]);
    }
  printf("]\n");
}

#if !defined(PREFETCHW)
#  if defined(__x86_64__) | defined(__i386__)
#    define PREFETCHW(x) asm volatile("prefetchw %0" :: "m" (*(unsigned long *)(x))) /* write */
#  elif defined(__sparc__)
#    define PREFETCHW(x) __builtin_prefetch((const void*) (x), 1, 3)
#  elif defined(__tile__)
#    include <tmc/alloc.h>
#    include <tmc/udn.h>
#    include <tmc/sync.h>
#    define PREFETCHW(x) tmc_mem_prefetch ((x), 64)
#  else
#    warning "You need to define PREFETCHW(x) for your architecture"
#  endif
#endif

/* 
 * 
 */
void* 
ssmem_alloc(ssmem_allocator_t* a, size_t size)
{
  void* m = NULL;

  /* 1st try to use from the collected memory */
  ssmem_free_set_t* cs = a->collected_set_list;
  if (cs != NULL)
    {
      m = (void*) cs->set[--cs->curr];
      PREFETCHW(m);

      if (cs->curr <= 0)
	{
	  a->collected_set_list = cs->set_next;
	  a->collected_set_num--;

	  ssmem_free_set_make_avail(a, cs);
	}
    }
  else
    {
      if ((a->mem_curr + size) >= a->mem_size)
	{
#if SSMEM_MEM_SIZE_DOUBLE == 1
	  a->mem_size <<= 1;
	  if (a->mem_size > SSMEM_MEM_SIZE_MAX)
	    {
	      a->mem_size = SSMEM_MEM_SIZE_MAX;
	    }
#endif
	  /* printf("[ALLOC] out of mem, need to allocate (chunk = %llu MB)\n", */
	  /* 	 a->mem_size / (1LL<<20)); */
	  if (size > a->mem_size)
	    {
	      /* printf("[ALLOC] asking for large mem. chunk\n"); */
	      while (a->mem_size < size)
		{
		  if (a->mem_size > SSMEM_MEM_SIZE_MAX)
		    {
		      fprintf(stderr, "[ALLOC] asking for memory chunk larger than max (%llu MB) \n",
			      SSMEM_MEM_SIZE_MAX / (1024 * 1024LL));
		      assert(a->mem_size <= SSMEM_MEM_SIZE_MAX);
		    }
		  a->mem_size <<= 1;
		}
	      /* printf("[ALLOC] new mem size chunk is %llu MB\n", a->mem_size / (1024 * 1024LL)); */
	    }
#if SSMEM_TRANSPARENT_HUGE_PAGES
	  int ret = posix_memalign(&a->mem, CACHE_LINE_SIZE, a->mem_size);
	  assert(ret == 0);
#else
	  a->mem = (void*) memalign(CACHE_LINE_SIZE, a->mem_size);
#endif
	  assert(a->mem != NULL);
#if SSMEM_ZERO_MEMORY == 1
	  memset(a->mem, 0, a->mem_size);
#endif

	  a->mem_curr = 0;
      
	  a->tot_size += a->mem_size;

	  a->mem_chunks = ssmem_list_node_new(a->mem, a->mem_chunks);
	}

      m = a->mem + a->mem_curr;
      a->mem_curr += size;
    }

#if SSMEM_TS_INCR_ON == SSMEM_TS_INCR_ON_ALLOC || SSMEM_TS_INCR_ON == SSMEM_TS_INCR_ON_BOTH
  ssmem_ts_next();
#endif
  return m;
}


/* return > 0 iff snew is > sold for each entry */
static int			
ssmem_ts_compare(size_t* s_new, size_t* s_old)
{
  int is_newer = 1;
  int i;
  for (i = 0; i < ssmem_ts_list_len; i++)
    {
      if (s_new[i] <= s_old[i])
	{
	  is_newer = 0;
	  break;
	}
    }
  return is_newer;
}

/* return > 0 iff s_1 is > s_2 > s_3 for each entry */
static int __attribute__((unused))
ssmem_ts_compare_3(size_t* s_1, size_t* s_2, size_t* s_3)
{
  int is_newer = 1;
  int i;
  for (i = 0; i < ssmem_ts_list_len; i++)
    {
      if (s_1[i] <= s_2[i] || s_2[i] <= s_3[i])
	{
	  is_newer = 0;
	  break;
	}
    }
  return is_newer;
}

static void ssmem_ts_set_print_no_newline(size_t* set);

/* 
 *
 */
int
ssmem_mem_reclaim(ssmem_allocator_t* a)
{
  if (__builtin_expect(a->released_num > 0, 0))
    {
      ssmem_released_t* rel_cur = a->released_mem_list;
      ssmem_released_t* rel_nxt = rel_cur->next;

      if (rel_nxt != NULL && ssmem_ts_compare(rel_cur->ts_set, rel_nxt->ts_set))
	{
	  rel_cur->next = NULL;
	  a->released_num = 1;
	  /* find and collect the memory */
	  do
	    {
	      rel_cur = rel_nxt;
	      free(rel_cur->mem);
	      free(rel_cur);
	      rel_nxt = rel_nxt->next;
	    }
	  while (rel_nxt != NULL);
	}
    }

  ssmem_free_set_t* fs_cur = a->free_set_list;
  if (fs_cur->ts_set == NULL)
    {
      return 0;
    }
  ssmem_free_set_t* fs_nxt = fs_cur->set_next;
  int gced_num = 0;

  if (fs_nxt == NULL || fs_nxt->ts_set == NULL)		/* need at least 2 sets to compare */
    {
      return 0;
    }

  if (ssmem_ts_compare(fs_cur->ts_set, fs_nxt->ts_set))
    {
      gced_num = a->free_set_num - 1;
      /* take the the suffix of the list (all collected free_sets) away from the
	 free_set list of a and set the correct num of free_sets*/
      fs_cur->set_next = NULL;
      a->free_set_num = 1;

      /* find the tail for the collected_set list in order to append the new 
	 free_sets that were just collected */
      ssmem_free_set_t* collected_set_cur = a->collected_set_list; 
      if (collected_set_cur != NULL)
	{
	  while (collected_set_cur->set_next != NULL)
	    {
	      collected_set_cur = collected_set_cur->set_next;
	    }

	  collected_set_cur->set_next = fs_nxt;
	}
      else
	{
	  a->collected_set_list = fs_nxt;
	}
      a->collected_set_num += gced_num;
    }

  /* if (gced_num) */
  /*   { */
  /*     printf("//collected %d sets\n", gced_num); */
  /*   } */
  return gced_num;
}

/* 
 *
 */
inline void 
ssmem_free(ssmem_allocator_t* a, void* obj)
{
  ssmem_free_set_t* fs = a->free_set_list;
  if (fs->curr == fs->size)
    {
      fs->ts_set = ssmem_ts_set_collect(fs->ts_set);
      ssmem_mem_reclaim(a);

      /* printf("[ALLOC] free_set is full, doing GC / size of garbage pointers: %10zu = %zu KB\n", garbagep, garbagep / 1024); */
      ssmem_free_set_t* fs_new = ssmem_free_set_get_avail(a, a->fs_size, a->free_set_list);
      a->free_set_list = fs_new;
      a->free_set_num++;
      fs = fs_new;

    }
  
  fs->set[fs->curr++] = (uintptr_t) obj;
#if SSMEM_TS_INCR_ON == SSMEM_TS_INCR_ON_FREE || SSMEM_TS_INCR_ON == SSMEM_TS_INCR_ON_BOTH
  ssmem_ts_next();
#endif
}

/* 
 *
 */
inline void 
ssmem_release(ssmem_allocator_t* a, void* obj)
{
  ssmem_released_t* rel_list = a->released_mem_list;
  ssmem_released_t* rel = ssmem_released_node_new(obj, rel_list);
  rel->ts_set = ssmem_ts_set_collect(rel->ts_set);
  int rn = ++a->released_num;
  a->released_mem_list = rel;
  if (rn >= SSMEM_GC_RLSE_SET_SIZE)
    {
      ssmem_mem_reclaim(a);
    }
}


/* 
 *
 */
static void 
ssmem_ts_set_print_no_newline(size_t* set)
{
  printf("[");
  if (set != NULL)
    {
      int i;
      for (i = 0; i < ssmem_ts_list_len; i++)
	{
	  printf("%zu|", set[i]);
	}
    }
  else
    {
      printf(" no timestamp yet ");
    }
  printf("]");

}

/* 
 *
 */
void
ssmem_free_list_print(ssmem_allocator_t* a)
{
  printf("[ALLOC] free_set list (%zu sets): \n", a->free_set_num);

  int n = 0;
  ssmem_free_set_t* cur = a->free_set_list;
  while (cur != NULL)
    {
      printf("(%-3d | %p::", n++, cur);
      ssmem_ts_set_print_no_newline(cur->ts_set);
      printf(") -> \n");
      cur = cur->set_next;
    }
  printf("NULL\n");
}

/* 
 *
 */
void
ssmem_collected_list_print(ssmem_allocator_t* a)
{
  printf("[ALLOC] collected_set list (%zu sets): \n", a->collected_set_num);

  int n = 0;
  ssmem_free_set_t* cur = a->collected_set_list;
  while (cur != NULL)
    {
      printf("(%-3d | %p::", n++, cur);
      ssmem_ts_set_print_no_newline(cur->ts_set);
      printf(") -> \n");
      cur = cur->set_next;
    }
  printf("NULL\n");
}

/* 
 *
 */
void
ssmem_available_list_print(ssmem_allocator_t* a)
{
  printf("[ALLOC] avail_set list: \n");

  int n = 0;
  ssmem_free_set_t* cur = a->available_set_list;
  while (cur != NULL)
    {
      printf("(%-3d | %p::", n++, cur);
      ssmem_ts_set_print_no_newline(cur->ts_set);
      printf(") -> \n");
      cur = cur->set_next;
    }
  printf("NULL\n");
}

/* 
 *
 */
void
ssmem_all_list_print(ssmem_allocator_t* a, int id)
{
  printf("[ALLOC] [%-2d] free_set list: %-4zu / collected_set list: %-4zu\n",
	 id, a->free_set_num, a->collected_set_num);
}

/* 
 *
 */
void
ssmem_ts_list_print()
{
  printf("[ALLOC] ts list (%u elems): ", ssmem_ts_list_len);
  ssmem_ts_t* cur = ssmem_ts_list;
  while (cur != NULL)
    {
      printf("(id: %-2zu / version: %zu) -> ", cur->id, cur->version);
      cur = cur->next;
    }

  printf("NULL\n"); 
}
