/*   
 *   File: snap_stress.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: 
 *   snap_stress.c is part of ASCYLIB
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
#include "utils.h"
#include "atomic_ops.h"
#ifdef __sparc__
#  include <sys/types.h>
#  include <sys/processor.h>
#  include <sys/procset.h>
#endif

#include "sspfd.h"
#include "clht_lf_res.h"
#include "ssmem.h"

/* #define DETAILED_THROUGHPUT */

/* ################################################################### *
 * GLOBALS
 * ################################################################### */


clht_hashtable_t** hashtable;
int num_buckets = 256;
int num_threads = 1;
int num_elements = KEY_BUCKT;
int duration = 1000;
int run_correctness = 0;
int print_vals_num = 0;
size_t pf_vals_num = 8191;
size_t obj_size = 4;

int seed = 0;
__thread unsigned long * seeds;
uint32_t rand_max;
#define rand_min 2

static volatile int stop;
__thread uint32_t phys_id;

volatile ticks *putting_succ;
volatile ticks *putting_fail;
volatile ticks *getting_succ;
volatile ticks *getting_fail;
volatile ticks *removing_succ;
volatile ticks *removing_fail;
volatile ticks *putting_count;
volatile ticks *putting_count_succ;
volatile ticks *getting_count;
volatile ticks *getting_count_succ;
volatile ticks *removing_count;
volatile ticks *removing_count_succ;
volatile ticks *total;


/* ################################################################### *
 * LOCALS
 * ################################################################### */

#ifdef DEBUG
extern __thread uint32_t put_num_restarts;
extern __thread uint32_t put_num_failed_expand;
extern __thread uint32_t put_num_failed_on_new;
#endif

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

#if !defined(COMPUTE_LATENCY)
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
  clht_t* ht;
} thread_data_t;


size_t 
math_pow(size_t n, size_t pow)
{
  size_t base = n;
  int r;
  for (r = 0; r < pow; r++)
    {
      base *= n;
    }
  return base;
}

size_t 
hash_rep(size_t key, size_t times)
{
  int r;
  for (r = 0; r < times; r++)
    {
      key = __ac_Jenkins_hash_64(key);
    }
  return key;
}


volatile clht_snapshot_t* snap;
volatile size_t key[KEY_BUCKT] = {0};
volatile uintptr_t val[KEY_BUCKT] = {0};
#define GET_VAL(v) (*(size_t*) (v))


void*
test(void* thread) 
{
  size_t num_retry_cas1 = 0, num_retry_cas2 = 0, num_retry_cas3 = 0 , num_retry_cas4 = 0, num_retry_cas5 = 0;
  thread_data_t* td = (thread_data_t*) thread;
  uint8_t ID = td->id;
  phys_id = the_cores[ID % (NUMBER_OF_SOCKETS * CORES_PER_SOCKET)];
  set_cpu(phys_id);

  ssmem_allocator_t* alloc = (ssmem_allocator_t*) memalign(CACHE_LINE_SIZE, sizeof(ssmem_allocator_t));
  assert(alloc != NULL);
  ssmem_alloc_init(alloc, SSMEM_DEFAULT_MEM_SIZE, ID);

  ssmem_gc_thread_init(alloc, ID);

  PF_INIT(3, SSPFD_NUM_ENTRIES, ID);

#if defined(COMPUTE_LATENCY)
  volatile ticks my_putting_succ = 0;
  volatile ticks my_putting_fail = 0;
  volatile ticks my_getting_succ = 0;
  volatile ticks my_getting_fail = 0;
  volatile ticks my_removing_succ = 0;
  volatile ticks my_removing_fail = 0;
#endif
  uint64_t my_putting_count = 0;
  uint64_t my_getting_count = 0;
  uint64_t my_removing_count = 0;

  uint64_t my_putting_count_succ = 0;
  uint64_t my_getting_count_succ = 0;
  uint64_t my_removing_count_succ = 0;
    
#if defined(COMPUTE_LATENCY) && PFD_TYPE == 0
  volatile ticks start_acq, end_acq;
  volatile ticks correction = getticks_correction_calc();
#endif
    
  seeds = seed_rand();
    
  MEM_BARRIER;

  barrier_cross(&barrier);

  barrier_cross(&barrier_global);

  size_t obj_size_bytes = obj_size * sizeof(size_t);
  volatile size_t* dat = (size_t*) malloc(obj_size_bytes);
  assert(dat != NULL);

  size_t* obj = NULL;

  while (stop == 0) 
    {
      size_t rand = (my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])));
      size_t k = (rand & 1) + 2;
      rand &= 1023;

      /* search baby! */

      int i;
      for (i = 0; i < KEY_BUCKT; i++)
	{
	  volatile uintptr_t v = val[i];
	  if (snap->map[i] == MAP_VALID && key[i] == k)
	    {
	      if (val[i] == v)
		{
		  if (GET_VAL(v) != k)
		    {
		      printf("[%02d] :get: key != val for %zu\n", ID, k);
		    }
		  break;
		}
	    }
	}

      if (rand > 513)
	{
	  my_putting_count++;

	  if (obj != NULL)
	    {
	      ssmem_free(alloc, (void*) obj);
	    }
	  obj = ssmem_alloc(alloc, 8);
	  *obj = k;


	  int empty_index = -2;
	  clht_snapshot_t s;

	retry:
	  s.snapshot = snap->snapshot;

	  int i;
	  for (i = 0; i < KEY_BUCKT; i++)
	    {
	      volatile uintptr_t v = val[i];
	      if (snap->map[i] == MAP_VALID && key[i] == k)
		{
		  if (val[i] == v)
		    {
		      if (empty_index > 0)
			{
			  snap->map[empty_index] = MAP_INVLD;
			}
		      goto end;
		    }
		}
	    }

	  clht_snapshot_all_t s1;
	  if (empty_index < 0)
	    {
	      empty_index = snap_get_empty_index(s.snapshot);
	      if (empty_index < 0)
		{
		  num_retry_cas1++;
		  goto end;
		}

	      s1 = snap_set_map(s.snapshot, empty_index, MAP_INSRT);
	      if (CAS_U64(&snap->snapshot, s.snapshot, s1) != s.snapshot)
		{
		  empty_index = -2;
		  num_retry_cas2++;
		  goto retry;
		}

	      val[empty_index] = (uintptr_t) obj;
	      key[empty_index] = k;
	    }
	  else
	    {
	      s1 = snap_set_map(s.snapshot, empty_index, MAP_INSRT);
	    }

	  clht_snapshot_all_t s2 = snap_set_map_and_inc_version(s1, empty_index, MAP_VALID);
	  if (CAS_U64(&snap->snapshot, s1, s2) != s1)
	    {
	      num_retry_cas3++;
	      /* key[empty_index] = 0; */
	      /* val[empty_index] = 0; */
	      goto retry;
	    }

	  obj = NULL;
	  my_putting_count_succ++;
	end:
	  ;
	}
      else
	{
	  my_removing_count++;
	  clht_snapshot_t s;

	retry_rem:
	  s.snapshot = snap->snapshot;

	  volatile uintptr_t v; 
	  int i, removed = 0;
	  for (i = 0; i < KEY_BUCKT && !removed; i++)
	    {
	      if (key[i] == k && s.map[i] == MAP_VALID)
		{
		  v = val[i];
		  clht_snapshot_all_t s1 = snap_set_map(s.snapshot, i, MAP_INVLD);
		  if (CAS_U64(&snap->snapshot, s.snapshot, s1) == s.snapshot)
		    {
		      /* snap->map[i] = MAP_INVLD; */
		      removed = 1;
		    }
		  else
		    {
		      num_retry_cas4++;
		      goto retry_rem;
		    }
		}
	    }
	  if (removed)
	    {
	      ssmem_free(alloc, (void*) v);
	      my_removing_count_succ++;
	    }
	}
    }

  free((void*) dat);
   
#if defined(DEBUG)
  if (put_num_restarts | put_num_failed_expand | put_num_failed_on_new)
    {
      /* printf("put_num_restarts = %3u / put_num_failed_expand = %3u / put_num_failed_on_new = %3u \n", */
      /* 	     put_num_restarts, put_num_failed_expand, put_num_failed_on_new); */
    }
#endif

  if (ID < 2)
    {
      printf("#retry-stats-thread-%d: #cas1: %-8zu / #cas2: %-8zu /"
	     "#cas3: %-8zu / #cas4: %-8zu / #cas5: %-8zu\n", 
	     ID, num_retry_cas1, num_retry_cas2, num_retry_cas3, num_retry_cas4, num_retry_cas5);
    }

  /* printf("gets: %-10llu / succ: %llu\n", num_get, num_get_succ); */
  /* printf("rems: %-10llu / succ: %llu\n", num_rem, num_rem_succ); */
  barrier_cross(&barrier);

#if defined(COMPUTE_LATENCY)
  putting_succ[ID] += my_putting_succ;
  putting_fail[ID] += my_putting_fail;
  getting_succ[ID] += my_getting_succ;
  getting_fail[ID] += my_getting_fail;
  removing_succ[ID] += my_removing_succ;
  removing_fail[ID] += my_removing_fail;
#endif
  putting_count[ID] += my_putting_count;
  getting_count[ID] += my_getting_count;
  removing_count[ID]+= my_removing_count;

  putting_count_succ[ID] += my_putting_count_succ;
  getting_count_succ[ID] += my_getting_count_succ;
  removing_count_succ[ID]+= my_removing_count_succ;

#if (PFD_TYPE == 1) && defined(COMPUTE_LATENCY)
  if (ID == 0)
    {
      printf("get ----------------------------------------------------\n");
      SSPFDPN(0, SSPFD_NUM_ENTRIES, print_vals_num);
      printf("put ----------------------------------------------------\n");
      SSPFDPN(1, SSPFD_NUM_ENTRIES, print_vals_num);
      printf("rem ----------------------------------------------------\n");
      SSPFDPN(2, SSPFD_NUM_ENTRIES, print_vals_num);

    }
#endif

  /* SSPFDTERM(); */

  pthread_exit(NULL);
}

int
main(int argc, char **argv) 
{
  set_cpu(the_cores[0]);
    
  assert(sizeof(clht_hashtable_t) == 2*CACHE_LINE_SIZE);

  struct option long_options[] = {
    // These options don't set a flag
    {"help",                      no_argument,       NULL, 'h'},
    {"duration",                  required_argument, NULL, 'd'},
    {"num-threads",               required_argument, NULL, 'n'},
    {"range",                     required_argument, NULL, 'r'},
    {"correctness",               no_argument, NULL, 'c'},
    {"num-buckets",               required_argument, NULL, 'b'},
    {"print-vals",                required_argument, NULL, 'v'},
    {"vals-pf",                   required_argument, NULL, 'f'},
    {"obj-size",                  required_argument, NULL, 's'},
    {NULL, 0, NULL, 0}
  };

  int i, c;
  while(1) 
    {
      i = 0;
      c = getopt_long(argc, argv, "hAf:d:i:n:r:s:u:m:a:l:p:b:v:f:c", long_options, &i);
		
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
	  printf("intset -- STM stress test "
		 "(linked list)\n"
		 "\n"
		 "Usage:\n"
		 "  intset [options...]\n"
		 "\n"
		 "Options:\n"
		 "  -h, --help\n"
		 "        Print this message\n"
		 "  -d, --duration <int>\n"
		 "        Test duration in milliseconds\n"
		 "  -n, --num-threads <int>\n"
		 "        Number of threads\n"
		 "  -r, --range <int>\n"
		 "        Range of integer values inserted in set\n"
		 "  -s, --obj-size <int>\n"
		 "        Size of the objects stored in the hash table\n"
		 "  -b, --num-buckets <int>\n"
		 "        Number of initial buckets (stronger than -l)\n"
		 "  -v, --print-vals <int>\n"
		 "        When using detailed profiling, how many values to print.\n"
		 "  -f, --val-pf <int>\n"
		 "        When using detailed profiling, how many values to keep track of.\n"
		 );
	  exit(0);
	case 'd':
	  duration = atoi(optarg);
	  break;
	case 'n':
	  num_threads = atoi(optarg);
	  break;
	case 's':
	  obj_size = atol(optarg);
	  break;
	case 'v':
	  print_vals_num = atoi(optarg);
	  break;
	case 'f':
	  pf_vals_num = pow2roundup(atoi(optarg)) - 1;
	  break;
	case '?':
	default:
	  printf("Use -h or --help for help\n");
	  exit(1);
	}
    }


  rand_max = num_elements;
    
  struct timeval start, end;
  struct timespec timeout;
  timeout.tv_sec = duration / 1000;
  timeout.tv_nsec = (duration % 1000) * 1000000;
    
  printf("//duration: sec: %lu, ns: %lu\n", timeout.tv_sec, timeout.tv_nsec);
  stop = 0;
    
  /* Initialize the hashtable */

  snap = (clht_snapshot_t*) memalign(CACHE_LINE_SIZE, CACHE_LINE_SIZE);
  assert(snap != NULL);

  /* Initializes the local data */
  putting_succ = (ticks *) calloc(num_threads , sizeof(ticks));
  putting_fail = (ticks *) calloc(num_threads , sizeof(ticks));
  getting_succ = (ticks *) calloc(num_threads , sizeof(ticks));
  getting_fail = (ticks *) calloc(num_threads , sizeof(ticks));
  removing_succ = (ticks *) calloc(num_threads , sizeof(ticks));
  removing_fail = (ticks *) calloc(num_threads , sizeof(ticks));
  putting_count = (ticks *) calloc(num_threads , sizeof(ticks));
  putting_count_succ = (ticks *) calloc(num_threads , sizeof(ticks));
  getting_count = (ticks *) calloc(num_threads , sizeof(ticks));
  getting_count_succ = (ticks *) calloc(num_threads , sizeof(ticks));
  removing_count = (ticks *) calloc(num_threads , sizeof(ticks));
  removing_count_succ = (ticks *) calloc(num_threads , sizeof(ticks));
    
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
    
  volatile ticks putting_suc_total = 1;
  volatile ticks putting_fal_total = 1;
  volatile ticks getting_suc_total = 1;
  volatile ticks getting_fal_total = 1;
  volatile ticks removing_suc_total = 1;
  volatile ticks removing_fal_total = 1;
  volatile uint64_t putting_count_total = 1;
  volatile uint64_t putting_count_total_succ = 1;
  volatile uint64_t getting_count_total = 1;
  volatile uint64_t getting_count_total_succ = 1;
  volatile uint64_t removing_count_total = 1;
  volatile uint64_t removing_count_total_succ = 1;
    
  for(t=0; t < num_threads; t++) 
    {
      putting_suc_total += putting_succ[t];
      putting_fal_total += putting_fail[t];
      getting_suc_total += getting_succ[t];
      getting_fal_total += getting_fail[t];
      removing_suc_total += removing_succ[t];
      removing_fal_total += removing_fail[t];
      putting_count_total += putting_count[t];
      putting_count_total_succ += putting_count_succ[t];
      getting_count_total += getting_count[t];
      getting_count_total_succ += getting_count_succ[t];
      removing_count_total += removing_count[t];
      removing_count_total_succ += removing_count_succ[t];
    }

  if(putting_count_total == 0) 
    {
      putting_suc_total = 0;
      putting_fal_total = 0;
      putting_count_total = 1;
      putting_count_total_succ = 2;
    }
    
  if(getting_count_total == 0) 
    {
      getting_suc_total = 0;
      getting_fal_total = 0;
      getting_count_total = 1;
      getting_count_total_succ = 2;
    }
    
  if(removing_count_total == 0) 
    {
      removing_suc_total = 0;
      removing_fal_total = 0;
      removing_count_total = 1;
      removing_count_total_succ = 2;
    }
    
#if defined(COMPUTE_LATENCY)
#  if defined(DEBUG)
  printf("#thread get_suc get_fal put_suc put_fal rem_suc rem_fal\n"); fflush(stdout);
#  endif
  printf("%d\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
	 num_threads,
	 getting_suc_total / getting_count_total_succ,
	 getting_fal_total / (getting_count_total - getting_count_total_succ),
	 putting_suc_total / putting_count_total_succ,
	 putting_fal_total / (putting_count_total - putting_count_total_succ),
	 removing_suc_total / removing_count_total_succ,
	 removing_fal_total / (removing_count_total - removing_count_total_succ)
	 );
#endif

    
#define LLU long long unsigned int

  int pr = (int) (putting_count_total_succ - removing_count_total_succ);
  int size_after = 0;

  for (i = 0; i < rand_max; i++)
    {
      size_after += (snap->map[i] == MAP_VALID);
    }
  printf("\n");
  printf("#Maps     | ");
  for (i = 0; i < rand_max; i++)
    {
      printf("#%d = %-5d | ", i, snap->map[i]);
    }
  printf("\n");
  printf("#Keys     | ");
  for (i = 0; i < rand_max; i++)
    {
      printf("#%d = %-5jd | ", i, key[i]);
    }
  printf("\n");
  printf("#Vals     | ");
  for (i = 0; i < rand_max; i++)
    {
      printf("#%d = %-5jd | ", i, (val[i]) ? *(size_t*) val[i] : -1);
    }
  printf("\n");

#if defined(DEBUG)
  printf("puts - rems  : %d\n", pr);
#endif
  /* assert(size_after == (pr)); */
  if (size_after != pr)
    {
      printf("######                                                                                      SIZE missmatch\n");
    }


  printf("    : %-10s | %-10s | %-11s | %s\n", "total", "success", "succ %", "total %");
  uint64_t total = putting_count_total + getting_count_total + removing_count_total;
  double putting_perc = 100.0 * (1 - ((double)(total - putting_count_total) / total));
  double removing_perc = 100.0 * (1 - ((double)(total - removing_count_total) / total));
  printf("puts: %-10llu | %-10llu | %10.1f%% | %.1f%%\n", (LLU) putting_count_total, 
	 (LLU) putting_count_total_succ,
	 (1 - (double) (putting_count_total - putting_count_total_succ) / putting_count_total) * 100,
	 putting_perc);
  printf("rems: %-10llu | %-10llu | %10.1f%% | %.1f%%\n", (LLU) removing_count_total, 
	 (LLU) removing_count_total_succ,
	 (1 - (double) (removing_count_total - removing_count_total_succ) / removing_count_total) * 100,
	 removing_perc);


  size_t all_total = putting_count_total + getting_count_total + removing_count_total;
  double throughput = (all_total) * 1000.0 / duration;
  printf("#txs tot (%zu\n", all_total);
  printf("#txs %-4d(%-10.0f = %.3f M\n", num_threads, throughput, throughput / 1e6);
    
    
  /* Last thing that main() should do */
  //printf("Main: program completed. Exiting.\n");
  pthread_exit(NULL);
    
  return 0;
}
