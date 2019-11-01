#pragma once
#include <immintrin.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <cstdint>
#include <iostream>
#include <time.h>

#ifdef PMEM
#include "libpmem.h"
#include "libpmemobj.h"
#endif

static constexpr const uint32_t kCacheLineSize = 64;

static bool FileExists(const char *pool_path) {
  struct stat buffer;
  return (stat(pool_path, &buffer) == 0);
}

#ifdef PMEM
#define CREATE_MODE_RW (S_IWUSR | S_IRUSR)
#ifndef POBJ_LAYOUT_TOID(allocator, char)
POBJ_LAYOUT_BEGIN(allocator);
POBJ_LAYOUT_TOID(allocator, char)
POBJ_LAYOUT_END(allocator)
#endif
#endif

#define LOG_FATAL(msg)      \
  std::cout << msg << "\n"; \
  exit(-1)

#define LOG(msg) std::cout << msg << "\n"

// ADD and SUB return the value after add or sub
#define ADD(_p, _v) (__atomic_add_fetch(_p, _v, __ATOMIC_SEQ_CST))
#define SUB(_p, _v) (__atomic_sub_fetch(_p, _v, __ATOMIC_SEQ_CST))
#define LOAD(_p) (__atomic_load_n(_p, __ATOMIC_SEQ_CST))
#define STORE(_p, _v) (__atomic_store_n(_p, _v, __ATOMIC_SEQ_CST))

int msleep(uint64_t msec){
  struct timespec ts;
  int res;

  if (msec < 0)
  {
      errno = EINVAL;
      return -1;
  }

  ts.tv_sec = msec / 1000;
  ts.tv_nsec = (msec % 1000) * 1000000;

  do {
      res = nanosleep(&ts, &ts);
  } while (res && errno == EINTR);

  return res;
}
