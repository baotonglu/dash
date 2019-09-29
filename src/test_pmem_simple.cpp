#include <sys/time.h>
#include <cstring>
#include <thread>
#include "../util/System.hpp"
#include "../util/random.h"
#include "allocator.h"
#include "ex_finger.h"
#include "libpmemobj.h"
#include "utils.h"

static const char *pool_name = "/mnt/pmem0/pmem_hash.data";
// static const char *pool_name = "pmem_hash.data";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;

Finger_EH<Key_t> *eh;
uint64_t *workload;
struct timeval tv1, tv2;

int main(int argc, char const *argv[]) {
  assert(argc >= 4);
  int initCap = atoi(argv[1]);
  int insert_num = atoi(argv[2]);
  int thread_num = atoi(argv[3]);

  Allocator::Initialize(pool_name, pool_size);
  assert(argc >= 4);

  eh = reinterpret_cast<Finger_EH<Key_t> *>(
      Allocator::GetRoot(sizeof(Finger_EH<Key_t>)));
  new (eh) Finger_EH<Key_t>(initCap);

  eh->pool_addr = Allocator::Get()->pm_pool_;

  /* Generate Workload for fixed_length key or variable_length key*/
  Allocator::ZAllocate((void **)&workload, kCacheLineSize,
                       sizeof(uint64_t) * (insert_num + 100) * 2);
  int i;
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL};
  int length = 4;
  init_by_array64(init, length);

  /* Generate Workload*/
  for (int i = 0; i < insert_num * 2 + 2; ++i) {
    uint64_t _key = genrand64_int64();
    workload[i] = _key;
  }

  std::cout << "Insertion begin" << std::endl;
  Key_t _key;
  Value_t _value;
  clock_t begin = clock();
  for (int i = 0; i < insert_num; ++i) {
    _key = workload[i];
    _value = reinterpret_cast<Value_t>(workload[i]);
    eh->Insert(_key, _value);
  }
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Insertion throughput: " << insert_num / time_spent << " ops/s"
            << std::endl;

  std::cout << "Postive search begin" << std::endl;
  uint64_t not_found = 0;
  begin = clock();
  for (int i = 0; i < insert_num; ++i) {
    _key = workload[i];
    if (eh->Get(_key) == NONE) {
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
  for (int i = insert_num; i < insert_num * 2; ++i) {
    _key = workload[i];
    if (eh->Get(_key) == NONE) {
      not_found++;
    }
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Not found = " << not_found << std::endl;
  std::cout << "Negative search throughput: " << insert_num / time_spent
            << " ops/s" << std::endl;

  std::cout << "Delete begin" << std::endl;
  not_found = 0;
  begin = clock();
  for (int i = 0; i < insert_num; ++i) {
    _key = workload[i];
    if (eh->Delete(_key) == false) {
      not_found++;
    }
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Not found = " << not_found << std::endl;
  std::cout << "Delete throughput: " << insert_num / time_spent << " ops/s"
            << std::endl;

  return 0;
}