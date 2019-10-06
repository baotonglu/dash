
#pragma once
#include <garbage_list.h>
#include <sys/mman.h>
#include "utils.h"
#include "x86intrin.h"

static const char* layout_name = "hashtable";
static const constexpr uint64_t pool_addr = 0x7ff700000000;

struct Allocator {
 public:
#ifdef PMEM
  static void Initialize(const char* pool_name, size_t pool_size) {
    instance_ = new Allocator(pool_name, pool_size);
    instance_->epoch_manager_.Initialize();
    instance_->garbage_list_.Initialize(&instance_->epoch_manager_, 1024);
    std::cout << "pool opened at: " << std::hex << instance_->pm_pool_
              << std::dec << std::endl;
  }

  static void ReInitialize_test_only(const char* pool_name, size_t pool_size) {
    pmemobj_close(instance_->pm_pool_);
    delete instance_;
    Allocator::Initialize(pool_name, pool_size);
  }

  Allocator(const char* pool_name, size_t pool_size) {
    if (!FileExists(pool_name)) {
      LOG("creating a new pool");
      pm_pool_ = pmemobj_create_addr(pool_name, layout_name, pool_size,
                                     CREATE_MODE_RW, (void*)pool_addr);
      if (pm_pool_ == nullptr) {
        LOG_FATAL("failed to create a pool;");
      }
      return;
    }
    LOG("opening an existing pool, and trying to map to same address");
    /* Need to open an existing persistent pool */
    pm_pool_ = pmemobj_open_addr(pool_name, layout_name, (void*)pool_addr);
    if (pm_pool_ == nullptr) {
      LOG_FATAL("failed to open the pool");
    }
  }

  PMEMobjpool* pm_pool_{nullptr};
  EpochManager epoch_manager_{};
  GarbageList garbage_list_{};

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

  static void NTWrite64(uint64_t* ptr, uint64_t val) {
    _mm_stream_si64((long long*)ptr, val);
  }

  static void NTWrite32(uint32_t* ptr, uint32_t val) {
    _mm_stream_si32((int*)ptr, val);
  }

  static PMEMobjpool* GetPool() { return instance_->pm_pool_; }

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
    // #ifdef PMEM
    //     auto oid_ptr = pmemobj_oid(ptr);
    //     TOID(char) ptr_cpy;
    //     TOID_ASSIGN(ptr_cpy, oid_ptr);
    //     POBJ_FREE(&ptr_cpy);
    // #else
    //     free(ptr);
    // #endif
  }
};

#ifdef PMEM
Allocator* Allocator::instance_ = nullptr;
#endif