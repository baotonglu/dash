
#pragma once
#include <sys/mman.h>
#include "utils.h"

static const char* layout_name = "hashtable";

struct allocator_root {
  uint64_t obj{0};
  uint64_t pool_addr{0};
};

struct Allocator {
 public:
#ifdef PMEM
  static void Initialize(const char* pool_name, size_t pool_size) {
    instance_ = new Allocator(pool_name, pool_size);
    auto root_oid = pmemobj_root(instance_->pm_pool_, sizeof(allocator_root));
    std::cout << instance_->pm_pool_ << std::endl;
    struct _pobj_pcache* cache = &_pobj_cached_pool;
    cache->pop = instance_->pm_pool_;
    cache->uuid_lo = root_oid.pool_uuid_lo;
    // cache->invalidate = 1;
    // std::cout << instance_->pm_pool_->root_offset << std::endl;
  }

  Allocator(const char* pool_name, size_t pool_size) {
    if (!FileExists(pool_name)) {
      LOG("creating a new pool");
      pm_pool_ =
          pmemobj_create(pool_name, layout_name, pool_size, CREATE_MODE_RW);
      if (pm_pool_ == nullptr) {
        LOG_FATAL("failed to create a pool;");
      }
      auto root_obj = reinterpret_cast<allocator_root*>(
          pmemobj_direct(pmemobj_root(pm_pool_, sizeof(allocator_root))));
      root_obj->pool_addr = reinterpret_cast<uint64_t>(pm_pool_);

      std::cout << std::hex << (uint64_t)pm_pool_ << std::endl;
      return;
    }

    LOG("opening an existing pool, and trying to map to same address");
    /* Need to open an existing persistent pool */
    // pm_pool_ = pmemobj_open(pool_name, layout_name);
    // if (pm_pool_ == nullptr) {
    // LOG_FATAL("failed to open the pool;");
    // }
    // auto root_obj_ptr = reinterpret_cast<allocator_root*>(
    // pmemobj_direct(pmemobj_root(pm_pool_, sizeof(allocator_root))));
    // auto root_obj = *root_obj_ptr;

    // if ((uint64_t)pm_pool_ == root_obj.pool_addr) {
    /* new pool address is the same as the old one, we are good */
    // LOG("successfully mapping to the same address");
    // return;
    // }
    // std::cout << std::hex << (uint64_t)pm_pool_ << std::endl;
    // std::cout << std::hex
    // << reinterpret_cast<allocator_root*>(pmemobj_direct(
    //  pmemobj_root(pm_pool_, sizeof(allocator_root))))
    // << std::endl;

    LOG("failed to map to the previous address, retrying...");
    // pmemobj_close(pm_pool_);

    int fd{0};
    if ((fd = open(pool_name, O_RDWR, 0666)) < 0) {
      perror("open");
      exit(1);
    }

    if ((errno = posix_fallocate(fd, 0, pool_size)) != 0) {
      perror("allocate");
      exit(1);
    }

    auto pool_addr = mmap(reinterpret_cast<void*>(0x7f4d80000000ull), pool_size,
                          PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
    if (reinterpret_cast<uint64_t>(pool_addr) != 0x7f4d80000000ull) {
      perror("mmap");
      exit(1);
    }
    auto root_addr = reinterpret_cast<allocator_root*>(pool_addr + 0x3c0550);
    pm_pool_ = reinterpret_cast<PMEMobjpool*>(pool_addr);
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

  static allocator_root* GetRoot() {
    auto root_obj = pmemobj_root(instance_->pm_pool_, sizeof(allocator_root));
    return reinterpret_cast<allocator_root*>(pmemobj_direct(root_obj));
  }

  static void Persist(void* ptr, size_t size) {
    pmemobj_persist(instance_->pm_pool_, ptr, size);
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
