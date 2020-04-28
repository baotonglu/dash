
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
#pragma once
#include <garbage_list.h>
#include <sys/mman.h>

#include "../util/utils.h"
#include "x86intrin.h"

static const char* layout_name = "hashtable";
static const constexpr uint64_t pool_addr = 0x5f0000000000;

typedef void (*DestroyCallback)(void* callback_context, void* object);

struct Allocator {
 public:
#ifdef PMEM
  static void Initialize(const char* pool_name, size_t pool_size) {
    instance_ = new Allocator(pool_name, pool_size);
    instance_->epoch_manager_.Initialize();
    instance_->garbage_list_.Initialize(&instance_->epoch_manager_,
                                        instance_->pm_pool_, 1024 * 8);
    std::cout << "pool opened at: " << std::hex << instance_->pm_pool_
              << std::dec << std::endl;
  }

  static void Close_pool() {
    pmemobj_close(instance_->pm_pool_);
    delete instance_;
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

  /* Must ensure that this pointer is in persistent memory*/
  static void Allocate(void** ptr, uint32_t alignment, size_t size,
                       int (*alloc_constr)(PMEMobjpool* pool, void* ptr,
                                           void* arg),
                       void* arg) {
    TX_BEGIN(instance_->pm_pool_) {
      pmemobj_tx_add_range_direct(ptr, sizeof(*ptr));
      *ptr = pmemobj_direct(pmemobj_tx_alloc(size, TOID_TYPE_NUM(char)));
      alloc_constr(instance_->pm_pool_, *ptr, arg);
    }
    TX_ONABORT { LOG_FATAL("Allocate: TXN Allocation Error"); }
    TX_END
  }

  static void Allocate(PMEMoid* pm_ptr, uint32_t alignment, size_t size,
                       int (*alloc_constr)(PMEMobjpool* pool, void* ptr,
                                           void* arg),
                       void* arg) {
    auto ret = pmemobj_alloc(instance_->pm_pool_, pm_ptr, size,
                             TOID_TYPE_NUM(char), alloc_constr, arg);
    if (ret) {
      LOG_FATAL("Allocate: Allocation Error in PMEMoid");
    }
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

  /*Must ensure that this pointer is in persistent memory*/
  static void ZAllocate(void** ptr, uint32_t alignment, size_t size) {
#ifdef PMEM
    TX_BEGIN(instance_->pm_pool_) {
      pmemobj_tx_add_range_direct(ptr, sizeof(*ptr));
      *ptr = pmemobj_direct(pmemobj_tx_zalloc(size, TOID_TYPE_NUM(char)));
    }
    TX_ONABORT { LOG_FATAL("ZAllocate: TXN Allocation Error"); }
    TX_END
#else
    posix_memalign(ptr, alignment, size);
    memset(*ptr, 0, size);
#endif
  }

  static void ZAllocate(PMEMoid* pm_ptr, uint32_t alignment, size_t size) {
    auto ret =
        pmemobj_zalloc(instance_->pm_pool_, pm_ptr, size, TOID_TYPE_NUM(char));

    if (ret) {
      std::cout << "Allocation size = " << size << std::endl;
      LOG_FATAL("allocation error");
    }
  }

  static void DefaultCallback(void* callback_context, void* ptr) {
#ifdef PMEM
    auto oid_ptr = pmemobj_oid(ptr);
    TOID(char) ptr_cpy;
    TOID_ASSIGN(ptr_cpy, oid_ptr);
    POBJ_FREE(&ptr_cpy);
#else
    free(ptr);
#endif
  }

  static void Free(void* ptr, DestroyCallback callback = DefaultCallback,
                   void* context = nullptr) {
    instance_->garbage_list_.Push(ptr, callback, context);
  }

  static void Free(GarbageList::Item* item, void* ptr,
                   DestroyCallback callback = DefaultCallback,
                   void* context = nullptr) {
    item->SetValue(ptr, instance_->epoch_manager_.GetCurrentEpoch(), callback,
                   context);
  }

  static EpochGuard AquireEpochGuard() {
    return EpochGuard{&instance_->epoch_manager_};
  }

  static void Protect() { instance_->epoch_manager_.Protect(); }

  static void Unprotect() { instance_->epoch_manager_.Unprotect(); }

  static GarbageList::Item* ReserveItem() {
    return instance_->garbage_list_.ReserveItem();
  }

  static void ResetItem(GarbageList::Item* mem) {
    instance_->garbage_list_.ResetItem(mem);
  }

  static void EpochRecovery() {
    instance_->garbage_list_.Recovery(&instance_->epoch_manager_,
                                      instance_->pm_pool_);
  }
};

#ifdef PMEM
Allocator* Allocator::instance_ = nullptr;
#endif
