
#pragma once
#include "utils.h"

static const char* layout_name = "hashtable";

struct Allocator {
 public:
#ifdef PMEM
  static void Initialize(const char* pool_name, size_t pool_size) {
    instance_ = new Allocator(pool_name, pool_size);
  }

  Allocator(const char* pool_name, size_t pool_size) {
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

  static void Allocate(void** ptr, uint32_t alignment, size_t size,
                       int (*alloc_constr)(PMEMobjpool* pool, void* ptr,
                                           void* arg),
                       void* arg) {
    PMEMoid pm_ptr;
    auto ret = pmemobj_alloc(instance_->pm_pool_, &pm_ptr, size,
                             TOID_TYPE_NUM(char), alloc_constr, arg);
    if (ret) {
      LOG_FATAL("allocation error");
    }
    *ptr = pmemobj_direct(pm_ptr);
  }

  static void* GetRoot(size_t size) {
    return pmemobj_direct(pmemobj_root(instance_->pm_pool_, size));
  }

  static void Persist(void* ptr, size_t size) {
    pmemobj_persist(instance_->pm_pool_, ptr, size);
  }

#endif

  static void Allocate(void** ptr, uint32_t alignment, size_t size) {
    posix_memalign(ptr, alignment, size);
  }

  static void ZAllocate(void** ptr, uint32_t alignment, size_t size) {
#ifdef PMEM
    PMEMoid pm_ptr;
    auto ret =
        pmemobj_zalloc(instance_->pm_pool_, &pm_ptr, size, TOID_TYPE_NUM(char));

    if (ret) {
      LOG_FATAL("allocation error");
    }
    /* FIXME: this should happen in a transaction
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
