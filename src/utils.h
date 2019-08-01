#pragma once
#include <stdlib.h>
#include <sys/stat.h>
#include <cstdint>
#include <iostream>

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
POBJ_LAYOUT_BEGIN(allocator);
POBJ_LAYOUT_TOID(allocator, char)
POBJ_LAYOUT_END(allocator)
#endif

#define LOG_FATAL(msg)      \
  std::cout << msg << "\n"; \
  exit(-1)

#define LOG(msg) std::cout << msg << "\n"
