#include <sys/time.h>
#include <cstring>
#include <thread>
#include "../util/System.hpp"
#include "../util/random.h"
#include "allocator.h"
#include "ex_finger.h"
#include "libpmemobj.h"
#include "utils.h"

/*
* This program is to test the correctness of the variable length key of extendible hashing
* To test whether the test program has the problem or the implementation of extendible hashing has problems
*/

//static const char *pool_name = "/mnt/pmem0/pmem_hash.data";
static const char *pool_name = "pmem_hash.data";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;
/*
struct string_arr{
  uint64_t key_[2];
  uint64_t length;
};
*/
Finger_EH<string_key *> *eh;
uint8_t *workload;
struct timeval tv1, tv2;

int main(int argc, char const *argv[]) {
  assert(argc >= 4);
  int initCap = atoi(argv[1]);
  uint64_t insert_num = atoi(argv[2]);
  int thread_num = atoi(argv[3]);

  std::cout << "InitCap = " << initCap << "; insert_num = " << insert_num << "; thread_num = " << thread_num << std::endl;

  Allocator::Initialize(pool_name, pool_size);
  assert(argc >= 4);

  eh = reinterpret_cast<Finger_EH<string_key *> *>(
      Allocator::GetRoot(sizeof(Finger_EH<string_key *>)));
  new (eh) Finger_EH<string_key *>(initCap);

  eh->pool_addr = Allocator::Get()->pm_pool_;

  /* Generate Workload for variable_length key*/
  //Allocator::ZAllocate((void **)&workload, kCacheLineSize,
  //                     20 * insert_num * 2);

  std::cout << "Allocate " << (double)(20ul * insert_num * 2ul) / 1024ul /1024ul / 1024ul <<  " GB Memory!"<< std::endl;
  workload = (uint8_t *)malloc(20ul * insert_num * 2ul);
  if(workload == NULL){
    std::cout << "Allocation Error!" << std::endl;
    return 0;
  }else{
    std::cout << "The initial addrs :" << std::hex << workload << std::endl; 
  }

  int i;
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL};
  int length = 4;
  init_by_array64(init, length);

  /* Generate Workload for variable length key*/
  uint64_t key_[2];
  for (uint64_t i = 0; i < insert_num * 2; ++i) {
    string_key * var_key = reinterpret_cast<string_key *>(workload + i * 20);
    var_key->length = 16;
    key_[0] = genrand64_int64();
    key_[1] = genrand64_int64();
    memcpy(var_key->key, key_, 16);
    if(i % 1000000 == 0){
      //std::cout <<" Current memory address is " << var_key << std::endl;
      std::cout << "Finish the initialization process at " << std::dec << i << std::endl;
      //std::cout << (double)((uint64_t)(&var_array[i]) - (uint64_t)(&var_array[0]))/1024/1024/1024 << "GB" << std::endl;
    }
  }

  std::cout << "Insertion begin" << std::endl;
  clock_t begin = clock();
  for (uint64_t i = 0; i < insert_num; ++i) {
    string_key * var_key = reinterpret_cast<string_key *>(workload + i * 20);
    eh->Insert(var_key, DEFAULT);
  }
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Insertion throughput: " << insert_num / time_spent << " ops/s"
            << std::endl;

  std::cout << "Positive search begin" << std::endl;
  uint64_t not_found = 0;
  begin = clock();
  for (uint64_t i = 0; i < insert_num; ++i) {
    string_key * var_key = reinterpret_cast<string_key *>(workload + i * 20);
    if (eh->Get(var_key) == NONE){
      not_found++;
    }
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Not found = " << not_found << std::endl;
  std::cout << "Postive search throughput: " << insert_num / time_spent
            << " ops/s" << std::endl;

  std::cout << "Negative search begin" << std::endl;
  not_found = 0;
  begin = clock();
  uint64_t total_num = insert_num * 2;
  for (uint64_t i = insert_num; i < total_num; ++i) {
    string_key * var_key = reinterpret_cast<string_key *>(workload + i * 20);
    if (eh->Get(var_key) == NONE) {
      not_found++;
    }
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Not found = " << not_found << std::endl;
  std::cout << "Negative search throughput: " << insert_num / time_spent
            << " ops/s" << std::endl;
  return 0;
}