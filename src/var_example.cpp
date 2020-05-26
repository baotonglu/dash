// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong.
// All rights reserved. Licensed under the MIT license.

// This is an example program about how to use Dash hash table in your
// application. We take Dash-EH as an exmaple. This example emphasizes the
// following aspects: (1) Which header files to include (2) How/where to set up
// pools (3) Which shared libraries to be linked: pthread; (4) How to
// use the basic operations supported by Dash

// Please check the CMakeLists.txt to know the details how to link the
// libraries.

#include "Hash.h"
#include "allocator.h"
#include "ex_finger.h"

/*variable-length key experiment*/
int main() {
  // Step 1: initialize the allocator
  Allocator::Initialize();

  // Step 2: Initialize the hash table
  // During initialization phase, allocate 64 segments for Dash-EH
  size_t segment_number = 64;
  Hash<string_key *> *hash_table = new extendible::Finger_EH<string_key *>(segment_number);
  //Hash<uint64_t> *hash_table = new new linear::Linear<T>();

  // Step 3: Generate string_key workload
  uint64_t variable_length = 16;
  uint64_t generate_num = 1024*1024;
  char* workload;
  workload = reinterpret_cast<char*>(malloc(generate_num * (variable_length + sizeof(string_key))));

  int word_num = (variable_length / 8) + (((variable_length % 8) != 0) ? 1 : 0);
  uint64_t *_key = reinterpret_cast<uint64_t *>(malloc(word_num * sizeof(uint64_t)));

  for (uint64_t i = 0; i < generate_num; ++i) {
    string_key* var_key = reinterpret_cast<string_key *>(workload +
                                             i * (variable_length + sizeof(string_key)));
    var_key->length = variable_length;
    for (int j = 0; j < word_num; ++j) {
      _key[j] = i;
    }
    memcpy(var_key->key, _key, variable_length);
  }

  // Step 4: Operate on the hash table
  // If using multi-threads, we need to use epoch for correct memory
  // reclamation; To make it simple, this example program only use one thread
  // but we still show how to use epoch mechanism; We enter into the epoch for
  // every 1024 operations to reduce the overhead; The following example inserts
  // 1024 * 1024 key-value pairs to the table, and then do the search and
  // delete operations
 uint64_t string_key_size = sizeof(string_key) + variable_length;

  // Insert
  int alread_exist = 0;
  for (uint64_t i = 0; i < 1024; ++i) {
    // Enroll into the epoch, if using one thread, epoch mechanism is actually
    // not needed
    auto epoch_guard = Allocator::AquireEpochGuard();
    for (uint64_t j = 0; j < 1024; ++j) {
      string_key* var_key = reinterpret_cast<string_key*>(workload + string_key_size * (i * 1024 + j));
      auto ret = hash_table->Insert(var_key, DEFAULT, true);
      if (ret == -1){
        alread_exist++;
      }
    }
  }

  std::cout << "already exist for first insert = " << alread_exist << std::endl;
  
  // testing for duplicate insert
  alread_exist = 0;
  for (uint64_t i = 0; i < 1024; ++i) {
    // Enroll into the epoch, if using one thread, epoch mechanism is actually
    // not needed
    auto epoch_guard = Allocator::AquireEpochGuard();
    for (uint64_t j = 0; j < 1024; ++j) {
      string_key* var_key = reinterpret_cast<string_key*>(workload + string_key_size * (i * 1024 + j));
      auto ret = hash_table->Insert(var_key, DEFAULT, true);
      if (ret == -1){
        alread_exist++;
      }
    }
  }

  std::cout << "already exist for second insert = " << alread_exist << std::endl;

  // Search
  Value_t value;
  uint64_t not_found = 0;
  uint64_t not_match_value = 0;
  for (uint64_t i = 0; i < 1024; ++i) {
    auto epoch_guard = Allocator::AquireEpochGuard();
    for (uint64_t j = 0; j < 1024; ++j) {
    string_key* var_key = reinterpret_cast<string_key*>(workload + string_key_size * (i * 1024 + j));
      if (hash_table->Get(var_key, &value, true) == NONE) {
        not_found++;
      }

      if (value != DEFAULT){
        not_match_value++;
      }
    }
  }
  std::cout << "The number of keys not found: " << not_found << std::endl;
  std::cout << "The number of value not match: " << not_match_value << std::endl;

  // Delete
  for (uint64_t i = 0; i < 1024; ++i) {
    auto epoch_guard = Allocator::AquireEpochGuard();
    for (uint64_t j = 0; j < 1024; ++j) {
    string_key* var_key = reinterpret_cast<string_key*>(workload + string_key_size * (i * 1024 + j));
      hash_table->Delete(var_key, true);
    }
  }

  return 0;
}
