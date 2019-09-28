/*   
 *   File: ssmem_test.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: 
 *   ssmem_test.c is part of ASCYLIB
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

#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sched.h>
#include <inttypes.h>
#include <sys/time.h>
#include <unistd.h>
#include <malloc.h>

#ifdef __sparc__
#  include <sys/types.h>
#  include <sys/processor.h>
#  include <sys/procset.h>
#endif

#include "ssmem.h"
#include "utils.h"

#define USE_MALLOC 0 		/* use malloc instead of ssmem */

/* atomic swap u64 */
#ifdef __sparc__
#  include <atomic.h>
#  define SWAP_U64(a, b) atomic_swap_64(a,b)
#elif defined(__x86_64__)
#  include <xmmintrin.h>
static inline uint64_t 
swap_uint64(volatile uint64_t* target,  uint64_t x) 
{
  __asm__ __volatile__("xchgq %0,%1"
		       :"=r" ((uint64_t) x)
		       :"m" (*(volatile uint64_t *)target), "0" ((uint64_t) x)
		       :"memory");

  return x;
}
#  define SWAP_U64(a, b) swap_uint64(a, b)
#elif defined(__tile__)
#  include <arch/atomic.h>
#  include <arch/cycle.h>
#  define SWAP_U64(a,b) arch_atomic_exchange(a,b)
#endif

/* ################################################################### *
 * GLOBALS
 * ################################################################### */


int num_allocs = 16;
int num_threads = 1;
int num_elements = 2048;
int duration = 1000;
int do_nothing = 0;
int do_releases = 0;
int seed = 0;
__thread unsigned long * seeds;
uint32_t rand_max;
#define rand_min 1

static volatile int stop;

/* ################################################################### *
 * BARRIER
 * ################################################################### */

typedef struct barrier 
{
  pthread_cond_t complete;
  pthread_mutex_t mutex;
  int count;
  int crossing;
} barrier_t;

void barrier_init(barrier_t *b, int n) 
{
  pthread_cond_init(&b->complete, NULL);
  pthread_mutex_init(&b->mutex, NULL);
  b->count = n;
  b->crossing = 0;
}

void barrier_cross(barrier_t *b) 
{
  pthread_mutex_lock(&b->mutex);
  /* One more thread through */
  b->crossing++;
  /* If not all here, wait */
  if (b->crossing < b->count) {
    pthread_cond_wait(&b->complete, &b->mutex);
  } else {
    pthread_cond_broadcast(&b->complete);
    /* Reset for next time */
    b->crossing = 0;
  }
  pthread_mutex_unlock(&b->mutex);
}
barrier_t barrier, barrier_global;

#define PFD_TYPE 0

#if defined(COMPUTE_THROUGHPUT)
#  define START_TS(s)
#  define END_TS(s, i)
#  define ADD_DUR(tar)
#  define ADD_DUR_FAIL(tar)
#  define PF_INIT(s, e, id)
#elif PFD_TYPE == 0
#  define START_TS(s)				\
  {						\
    asm volatile ("");				\
    start_acq = getticks();			\
    asm volatile ("");
#  define END_TS(s, i)				\
    asm volatile ("");				\
    end_acq = getticks();			\
    asm volatile ("");				\
    }

#  define ADD_DUR(tar) tar += (end_acq - start_acq - correction)
#  define ADD_DUR_FAIL(tar)					\
  else								\
    {								\
      ADD_DUR(tar);						\
    }
#  define PF_INIT(s, e, id)
#else
#  define SSPFD_NUM_ENTRIES  pf_vals_num
#  define START_TS(s)      SSPFDI(s)
#  define END_TS(s, i)     SSPFDO(s, i & SSPFD_NUM_ENTRIES)

#  define ADD_DUR(tar) 
#  define ADD_DUR_FAIL(tar)
#  define PF_INIT(s, e, id) SSPFDINIT(s, e, id)
#endif

typedef struct thread_data
{
  uint8_t id;
} thread_data_t;



volatile uint64_t total_ops = 0;
uintptr_t* array;
uintptr_t* array_obj;

void*
test(void* thread) 
{
  thread_data_t* td = (thread_data_t*) thread;
  uint8_t ID = td->id;

  seeds = seed_rand();

  /* printf("[%2d] starting:: %d allocs\n", ID, num_allocs); */

  ssmem_allocator_t* alloc = (ssmem_allocator_t*) memalign(CACHE_LINE_SIZE, 
							   num_allocs * sizeof(ssmem_allocator_t));
  assert(alloc != NULL);

  int i;
  for (i = 0; i < num_allocs; i++)
    {
      ssmem_alloc_init(alloc + i, SSMEM_DEFAULT_MEM_SIZE, ID);
    }

  size_t ops = 0;
  barrier_cross(&barrier_global);

  if (do_nothing)
    {
      while (stop == 0) 
	{
	  ops++;
	  int a    = (my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])) % num_allocs);
	  int b    = (my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])) % num_allocs);
	  
	  if (do_releases && (ops & 31) == 0)
	    {
	      size_t* obj_rel = (size_t*) malloc(sizeof(uintptr_t));
	      ssmem_release(alloc + b, (void*) obj_rel);
	    }

	  size_t* obj = (size_t*) ssmem_alloc(alloc + a, sizeof(uintptr_t));
	  ssmem_free(alloc + a, (void*) obj);
#if SSMEM_TS_INCR_ON == SSMEM_TS_INCR_ON_NONE
	  SSMEM_SAFE_TO_RECLAIM();
#endif
	}
    }
  else
    {
      while (stop == 0) 
	{
	  ops++;
	  int spot = (my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])) % rand_max);
	  int a    = (my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])) % num_allocs);

	  if (do_releases && (ops & 127) == 0)
	    {
	      size_t* obj_rel = (size_t*) malloc(sizeof(uintptr_t));
	      *obj_rel = *(size_t*) array_obj[spot];
	      size_t* old_rel = (size_t*) SWAP_U64((uint64_t*) (array_obj + spot), (uint64_t) obj_rel);
	      ssmem_release(alloc + a, (void*) old_rel);
	    }

	  /* Mem. reclaimation correctness test ------------------------------------------ */
	  // Idea: Array `array` contains memory objects that contain the value of the index
	  // of each array spot (i.e., *array[i] == i).
	  // On each iteration, threads read an array spot in `ref`.
	  // Their goal is to dereference (get the contents of `ref` and write
	  // the value in a newly allocated object `obj`. They then atomically swap
	  // `obj` with the current memory object in the array spot they are accessing.
	  // Naturally, the atomic swap return the previous memory object in the array
	  // spot. This `old` object is of course freed to avoid memory leaks. Due to
	  // this freeing, we might have a data race on dereferencing `ref` if the
	  // allocator does not handle memory reclaimation.

	  // Possible problem: Thread-1 first fetches `ref` from `spot` and later deferences 
	  // `ref`. Between these two memory accesses, another thread might use the same
	  // `spot`, thus freeing the memory that `ref` points at. This memory can now
	  // (without memory reclaiamation) be reused by an allocation, that will write a
	  // new, possibly different value in `*ref`. When Thread-1 dereferences `ref`, it
	  // it will get an incorrect value for its newly allocated `obj`.

	  /* ----------------------------------------------------------------------------- */
#if USE_MALLOC == 1
	  size_t* obj = (size_t*) malloc(sizeof(uintptr_t));
#else
	  size_t* obj = (size_t*) ssmem_alloc(alloc + a, sizeof(uintptr_t));
#endif
	  size_t* ref = (size_t*) array[spot];
	  *obj = *ref;

	  size_t* old = (size_t*) SWAP_U64((uint64_t*) (array + spot), (uint64_t) obj);

#if USE_MALLOC == 1
	  free((void*) old);
#else
	  ssmem_free(alloc + a, (void*) old);
#endif

#if SSMEM_TS_INCR_ON == SSMEM_TS_INCR_ON_NONE
	  SSMEM_SAFE_TO_RECLAIM();
#endif
	}
    }

  __sync_fetch_and_add(&total_ops, ops);

  /* printf("[%2d] stoping...\n", ID); */
  barrier_cross(&barrier);
  if (ID == 0)
    {
      int i;
      for (i = 0; i < rand_max; i++)
	{
	  if (array[i] != 0)
	    {
	      size_t v = *(size_t*) array[i];
	      if (v != i)
		{
		  printf("!! %d = %zu\n", i, v);
		}
	      v  = *(size_t*) array_obj[i];
	      if (v != i)
		{
		  printf("!! %d = %zu (release mem)\n", i, v);
		}
	    }
	}
    }
  barrier_cross(&barrier_global);
  ssmem_term();
  pthread_exit(NULL);
}

int
main(int argc, char **argv) 
{
  struct option long_options[] = {
    // These options don't set a flag
    {"help",                      no_argument,       NULL, 'h'},
    {"nothing",                   no_argument,       NULL, 't'},
    {"release",                   no_argument,       NULL, 'e'},
    {"duration",                  required_argument, NULL, 'd'},
    {"initial-size",              required_argument, NULL, 'i'},
    {"num-threads",               required_argument, NULL, 'n'},
    {"range",                     required_argument, NULL, 'r'},
    {NULL, 0, NULL, 0}
  };

  size_t initial = 8, range = 2048;

  int i, c;
  while(1) 
    {
      i = 0;
      c = getopt_long(argc, argv, "hed:i:n:r:t", long_options, &i);
		
      if(c == -1)
	break;
		
      if(c == 0 && long_options[i].flag == 0)
	c = long_options[i].val;
		
      switch(c) 
	{
	case 0:
	  /* Flag is automatically set */
	  break;
	case 'h':
	  printf("ssmem_test -- ssmem allocator stress test \n"
		 "Usage:\n"
		 "  ./ssmem_test [options...]\n"
		 "\n"
		 "Options:\n"
		 "  -h, --help\n"
		 "        Print this message\n"
		 "  -d, --duration <int>\n"
		 "        Test duration in milliseconds\n"
		 "  -i, --initial-size <int>\n"
		 "        Number of allocators\n"
		 "  -n, --num-threads <int>\n"
		 "        Number of threads\n"
		 "  -r, --range <int>\n"
		 "        Range of integer values inserted in set for testing\n"
		 "  -t, --nothing\n"
		 "        Do nothing but alloc/free\n"
		 "  -e, --release\n"
		 "        Additionally do mallocs and ssmem_releases in some rounds\n"
		 );
	  exit(0);
	case 'd':
	  duration = atoi(optarg);
	  break;
	case 'i':
	  initial = atoi(optarg);
	  break;
	case 'n':
	  num_threads = atoi(optarg);
	  break;
	case 't':
	  do_nothing = 1;
	  break;
	case 'e':
	  do_releases = 1;
	  break;
	case 'r':
	  range = atol(optarg);
	  break;
	case '?':
	default:
	  printf("Use -h or --help for help\n");
	  exit(1);
	}
    }

  num_allocs = initial;

  printf("# allocators: %d / threads: %d / range: %zu / do releases: %d / do nothing: %d\n", 
	 num_allocs, num_threads, range, do_releases, do_nothing);

  struct timeval start, end;
  struct timespec timeout;
  timeout.tv_sec = duration / 1000;
  timeout.tv_nsec = (duration % 1000) * 1000000;
    
  stop = 0;
    

  array = (uintptr_t*) calloc(range, sizeof(uintptr_t));
  assert(array != NULL);
  array_obj = (uintptr_t*) calloc(range, sizeof(uintptr_t));
  assert(array_obj != NULL);

  int j;
  for (j = 0; j < range; j++)
    {
      uintptr_t* obj = (uintptr_t*) malloc(sizeof(uintptr_t));
      assert(obj != NULL);
      *obj = j;
      array[j] = (uintptr_t) obj;
      size_t* s = (size_t*) malloc(sizeof(size_t));
      assert(s != NULL);
      *s = j;
      array_obj[j] = (uintptr_t) s;
    }

  rand_max = range;


  /* Initialize the hashtable */

  pthread_t threads[num_threads];
  pthread_attr_t attr;
  int rc;
  void *status;
    
  barrier_init(&barrier_global, num_threads + 1);
  barrier_init(&barrier, num_threads);
    
  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
  thread_data_t* tds = (thread_data_t*) malloc(num_threads * sizeof(thread_data_t));

  long t;
  for(t = 0; t < num_threads; t++)
    {
      tds[t].id = t;
      rc = pthread_create(&threads[t], &attr, test, tds + t);
      if (rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
        
    }
    
  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);
    
  barrier_cross(&barrier_global);
  gettimeofday(&start, NULL);
  nanosleep(&timeout, NULL);

  stop = 1;
  barrier_cross(&barrier_global);
  gettimeofday(&end, NULL);
  duration = (end.tv_sec * 1000 + end.tv_usec / 1000) - (start.tv_sec * 1000 + start.tv_usec / 1000);
    
  for(t = 0; t < num_threads; t++) 
    {
      rc = pthread_join(threads[t], &status);
      if (rc) 
	{
	  printf("ERROR; return code from pthread_join() is %d\n", rc);
	  exit(-1);
	}
    }

  free(tds);

  free(array);
  for (j = 0; j < range; j++)
    {
      free((void*) array_obj[j]);
    }
  free(array_obj);
  double throughput = total_ops * 1000.0 / duration;
  printf("#throughput with %-4d threads: %10.0f = %.3f M\n", num_threads, throughput, throughput / 1.e6);

  pthread_exit(NULL);
  return 0;
}
