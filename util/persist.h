#ifndef UTIL_PERSIST_H_
#define UTIL_PERSIST_H_

#include <cstdlib>
#include <stdint.h>

#define CPU_FREQ_MHZ (1994)  // cat /proc/cpuinfo
#define CAS(_p, _u, _v)  (__atomic_compare_exchange_n (_p, _u, _v, false, __ATOMIC_ACQUIRE, __ATOMIC_ACQUIRE))
#define CACHE_LINE (64)

extern uint64_t kWriteLatencyInNS;
extern uint64_t clflushCount;

static inline void CPUPause(void) {
  __asm__ volatile("pause":::"memory");
}

static inline unsigned long ReadTSC(void) {
  unsigned long var;
  unsigned int hi, lo;
  asm volatile("rdtsc":"=a"(lo),"=d"(hi));
  var = ((unsigned long long int) hi << 32) | lo;
  return var;
}

inline void mfence(void) {
  asm volatile("mfence":::"memory");
}

inline void flush(char* ptr){
  asm volatile("clflush %0" : "+m" (*(volatile char*)ptr));
}
/*
inline void clflush(char* data, size_t len) {
  volatile char *ptr = (char*)((unsigned long)data & (~(CACHE_LINE-1)));
  mfence();
  for (; ptr < data+len; ptr+=CACHE_LINE) {
    unsigned long etcs = ReadTSC() + (unsigned long) (kWriteLatencyInNS*CPU_FREQ_MHZ/1000);
    asm volatile("clflush %0" : "+m" (*(volatile char*)ptr));
    while (ReadTSC() < etcs) CPUPause();
    clflushCount++;
  }
  mfence();
}
*/

inline void clflush(char* data, size_t len) {
  volatile char *ptr = (char*)((unsigned long)data & (~(CACHE_LINE-1)));
  mfence();
  for (; ptr < data+len; ptr+=CACHE_LINE) {
    asm volatile("clflush %0" : "+m" (*(volatile char*)ptr));
    //clflushCount++;
  }
  mfence();
}

#endif  // UTIL_PERSIST_H_
