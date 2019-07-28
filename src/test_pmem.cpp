#include "libpmemobj.h"
#include "pm_allocator.h"
#include "utils.h"

static const char *pool_name = "pmem_hash.data";
static const uint32_t pool_size = 1024 * 1024 * 1024;

int main() {
  Allocator::Initialize(pool_name, pool_size);
  return 0;
}
