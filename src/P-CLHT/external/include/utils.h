#ifndef _UTILS_H_INCLUDED_
#define _UTILS_H_INCLUDED_
//some utility functions
//#define _GNU_SOURCE
//#define USE_MUTEX_LOCKS
//#define ADD_PADDING
//#define OPTERON
//#define OPTERON_OPTIMIZE

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
#elif defined(__tile__)
#  include <arch/atomic.h>
#  include <arch/cycle.h>
#  include <tmc/cpus.h>
#  include <tmc/task.h>
#  include <tmc/spin.h>
#  include <sched.h>
#else
#  include <emmintrin.h>
#  include <xmmintrin.h>
#  ifdef NUMA
#    include <numa.h>
#  endif
#endif
#include <pthread.h>

#include "sspfd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ALIGNED(N) __attribute__ ((aligned (N)))

#ifdef __sparc__
#  define PAUSE    asm volatile("rd    %%ccr, %%g0\n\t" \
				::: "memory")

#elif defined(__tile__)
#  define PAUSE cycle_relax()
#else
#  define PAUSE _mm_pause()
#endif
  static inline void
  pause_rep(uint32_t num_reps)
  {
    uint32_t i;
    for (i = 0; i < num_reps; i++)
      {
	PAUSE;
	/* PAUSE; */
	/* asm volatile ("NOP"); */
      }
  }

  static inline void
  nop_rep(uint32_t num_reps)
  {
    uint32_t i;
    for (i = 0; i < num_reps; i++)
      {
	asm volatile ("NOP");
      }
  }


  //machine dependent parameters
#  define NUMBER_OF_SOCKETS 1
#  define CORES_PER_SOCKET 8
#  define CACHE_LINE_SIZE 64
#  define NOP_DURATION 2
#  define USE_HYPERTRHEADS 1
  static uint8_t __attribute__ ((unused)) the_cores[] = {
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 9, 10, 11, 12, 13, 14, 15,
  };

  /* PLATFORM specific -------------------------------------------------------------------- */
#if defined(__x86_64__)
#  define PREFETCHW(x)		     asm volatile("prefetchw %0" :: "m" (*(unsigned long *)x))
#  define PREFETCH(x)		     asm volatile("prefetch %0" :: "m" (*(unsigned long *)x))
#elif defined(__sparc__)
#  define PREFETCHW(x)		
#elif defined(XEON)
#  define PREFETCHW(x)		
#else
#  define PREFETCHW(x)		
#endif

  //debugging functions
#ifdef DEBUG
#  define DPRINT(args...) fprintf(stderr,args);
#  define DDPRINT(fmt, args...) printf("%s:%s:%d: "fmt, __FILE__, __FUNCTION__, __LINE__, args)
#else
#  define DPRINT(...)
#  define DDPRINT(fmt, ...)
#endif



  static inline int get_cluster(int thread_id) {
#ifdef __solaris__
    if (thread_id>64){
      perror("Thread id too high");
      return 0;
    }
    return thread_id/CORES_PER_SOCKET;
    //    return the_sockets[thread_id];
#elif XEON
    if (thread_id>=80){
      perror("Thread id too high");
      return 0;
    }

    return the_sockets[thread_id];    
#elif defined(__tile__)
    return 0;
#else
    return thread_id/CORES_PER_SOCKET;
#endif
  }

  typedef uint64_t ticks;

  static inline double wtime(void)
  {
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)t.tv_sec + ((double)t.tv_usec)/1000000.0;
  }

  static inline void set_cpu(int cpu) {
#if !defined(NO_SET_CPU)
#  ifdef __sparc__
    processor_bind(P_LWPID,P_MYID, cpu, NULL);
#  elif defined(__tile__)
    if (cpu>=tmc_cpus_grid_total()) {
      perror("Thread id too high");
    }
    // cput_set_t cpus;
    if (tmc_cpus_set_my_cpu(cpu)<0) {
      tmc_task_die("tmc_cpus_set_my_cpu() failed."); 
    }    
#  else
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu, &mask);
    /* numa_set_preferred(get_cluster(cpu)); */
    pthread_t thread = pthread_self();
    if (pthread_setaffinity_np(thread, sizeof(cpu_set_t), &mask) != 0) {
      fprintf(stderr, "Error setting thread affinity\n");
    }
#  endif
#endif
  }

  /* #if defined(__i386__) */
  /*   static inline ticks getticks(void) { */
  /*     ticks ret; */

  /*     __asm__ __volatile__("rdtsc" : "=A" (ret)); */
  /*     return ret; */
  /*   } */
  /* #elif defined(__x86_64__) */
  /*   static inline ticks getticks(void) */
  /*   { */
  /*     unsigned hi, lo; */
  /*     __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi)); */
  /*     return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 ); */
  /*   } */
  /* #elif defined(__sparc__) */
  /*   static inline ticks getticks(){ */
  /*     ticks ret; */
  /*     __asm__ __volatile__ ("rd %%tick, %0" : "=r" (ret) : "0" (ret));  */
  /*     return ret; */
  /*   } */
  /* #elif defined(__tile__) */
  /*   static inline ticks getticks(){ */
  /*     return get_cycle_count(); */
  /*   } */
  /* #endif */

  static inline void cdelay(ticks cycles){
    ticks __ts_end = getticks() + (ticks) cycles;
    while (getticks() < __ts_end);
  }

  static inline void cpause(ticks cycles){
#if defined(XEON)
    cycles >>= 3;
    ticks i;
    for (i=0;i<cycles;i++) {
      _mm_pause();
    }
#else
    ticks i;
    for (i=0;i<cycles;i++) {
      __asm__ __volatile__("nop");
    }
#endif
  }

  static inline void udelay(unsigned int micros)
  {
    double __ts_end = wtime() + ((double) micros / 1000000);
    while (wtime() < __ts_end);
  }

  //getticks needs to have a correction because the call itself takes a
  //significant number of cycles and skewes the measurement
  static inline ticks getticks_correction_calc() {
#define GETTICKS_CALC_REPS 5000000
    ticks t_dur = 0;
    uint32_t i;
    for (i = 0; i < GETTICKS_CALC_REPS; i++) {
      ticks t_start = getticks();
      ticks t_end = getticks();
      t_dur += t_end - t_start;
    }
    //    printf("corr in float %f\n", (t_dur / (double) GETTICKS_CALC_REPS));
    ticks getticks_correction = (ticks)(t_dur / (double) GETTICKS_CALC_REPS);
    return getticks_correction;
  }

  static inline ticks get_noop_duration() {
#define NOOP_CALC_REPS 1000000
    ticks noop_dur = 0;
    uint32_t i;
    ticks corr = getticks_correction_calc();
    ticks start;
    ticks end;
    start = getticks();
    for (i=0;i<NOOP_CALC_REPS;i++) {
      __asm__ __volatile__("nop");
    }
    end = getticks();
    noop_dur = (ticks)((end-start-corr)/(double)NOOP_CALC_REPS);
    return noop_dur;
  }

  /// Round up to next higher power of 2 (return x if it's already a power
  /// of 2) for 32-bit numbers
  static inline uint64_t pow2roundup (uint64_t x){
    if (x==0) return 1;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return x+1;
  }
#define my_random xorshf96

  /* 
   * Returns a pseudo-random value in [1;range).
   * Depending on the symbolic constant RAND_MAX>=32767 defined in stdlib.h,
   * the granularity of rand() could be lower-bounded by the 32767^th which might 
   * be too high for given values of range and initial.
   */
  static inline long rand_range(long r) {
    int m = RAND_MAX;
    long d, v = 0;

    do {
      d = (m > r ? r : m);
      v += 1 + (long) (d * ((double) rand() / ((double) (m) + 1.0)));
      r -= m;
    } while (r > 0);
    return v;
  }
    
  //fast but weak random number generator for the sparc machine
  static inline uint32_t fast_rand() {
    return ((getticks()&4294967295)>>4);
  }


  static inline unsigned long* seed_rand() {
    unsigned long* seeds;
    seeds = (unsigned long*) memalign(64, 64);
    seeds[0] = getticks() % 0x992123E456789LL;
    seeds[1] = getticks() % 0x22136D2436069LL;
    seeds[2] = getticks() % 0x2119F521288629LL;
    return seeds;
  }

  //Marsaglia's xorshf generator
  static inline unsigned long xorshf96(unsigned long* x, unsigned long* y, unsigned long* z) {          //period 2^96-1
    unsigned long t;
    (*x) ^= (*x) << 16;
    (*x) ^= (*x) >> 5;
    (*x) ^= (*x) << 1;

    t = *x;
    (*x) = *y;
    (*y) = *z;
    (*z) = t ^ (*x) ^ (*y);

    return *z;
  }

#ifdef __cplusplus
}

#endif


#endif
