
#pragma once
#include "utils.h"

static const char* layout_name = "hashtable";

struct Allocator {
 public:
#ifdef PMEM
  static void Initialize(const char* pool_name, uint32_t pool_size) {
    instance_ = new Allocator(pool_name, pool_size);
  }

  Allocator(const char* pool_name, uint32_t pool_size) {
    if (!FileExists(pool_name)) {
      LOG("creating a new pool");
      pm_pool_ =
          pmemobj_create(pool_name, layout_name, pool_size, CREATE_MODE_RW);
      if (pm_pool_ == nullptr) {
        LOG_FATAL("failed to create a pool;");
      }
    } else {
      LOG("opening an existing pool");
      pm_pool_ = pmemobj_open(pool_name, layout_name);
      if (pm_pool_ == nullptr) {
        LOG_FATAL("failed to open the pool;");
      }
    }
  }

  PMEMobjpool* pm_pool_{nullptr};
  static Allocator* instance_;
  static Allocator* Get() { return instance_; }
#endif

  static void ZAllocate(void** ptr, uint32_t alignment, size_t size) {
#ifdef PMEM
    PMEMoid pm_ptr;
    auto ret =
        pmemobj_zalloc(instance_->pm_pool_, &pm_ptr, size, TOID_TYPE_NUM(char));

    if (ret) {
      LOG_FATAL("allocation error");
    }
    /* this should happen in a transaction
     */
    *ptr = pmemobj_direct(pm_ptr);
#else
    posix_memalign(ptr, alignment, size);
    memset(*ptr, 0, size);
#endif
  }

  static void Free(void* ptr) {
#ifdef PMEM
    auto oid_ptr = pmemobj_oid(ptr);
    TOID(char) ptr_cpy;
    TOID_ASSIGN(ptr_cpy, oid_ptr);
    POBJ_FREE(&ptr_cpy);
#else
    free(ptr);
#endif
  }
};

#ifdef PMEM
Allocator* Allocator::instance_ = nullptr;
#endif
