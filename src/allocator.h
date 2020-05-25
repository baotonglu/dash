
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
#pragma once
#include <garbage_list.h>
#include <sys/mman.h>
#include <cstring>

#include "../util/utils.h"
#include "x86intrin.h"

typedef void (*DestroyCallback)(void* callback_context, void* object);

struct Allocator {
 public:
  static void Initialize() {
    instance_ = new Allocator();
    instance_->epoch_manager_.Initialize();
    instance_->garbage_list_.Initialize(&instance_->epoch_manager_, 1024 * 8);
  }

  static void Close_pool() {
    delete instance_;
  }

  Allocator() {
  }

  EpochManager epoch_manager_{};
  GarbageList garbage_list_{};

  static Allocator* instance_;
  static Allocator* Get() { return instance_; }

  static void Allocate(void** ptr, uint32_t alignment, size_t size) {
    posix_memalign(ptr, alignment, size);
  }

  /*Must ensure that this pointer is in persistent memory*/
  static void ZAllocate(void** ptr, uint32_t alignment, size_t size) {
    posix_memalign(ptr, alignment, size);
    memset(*ptr, 0, size);
  }

  static void DefaultCallback(void* callback_context, void* ptr) {
    free(ptr);
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
};

Allocator* Allocator::instance_ = nullptr;