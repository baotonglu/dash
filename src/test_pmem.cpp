#include <sys/time.h>
#include <thread>
#include "../util/System.hpp"
#include "../util/random.h"
#include "allocator.h"
#include "ex_finger.h"
#include "libpmemobj.h"
#include "utils.h"

static const char *pool_name = "pmem_hash.data";
static const uint32_t pool_size = 1024 * 1024 * 1024;

Finger_EH *eh;
uint64_t **workload;
struct timeval tv1, tv2;

struct range {
  uint64_t index;
  uint64_t begin;
  uint64_t end;
};

void concurr_insert(struct range *_range) {
  size_t key;
  char arr[64];
  Value_t value = (Value_t)arr;
  uint64_t *array = workload[_range->index];

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
    key = i;
    eh->Insert(key, value);
  }
}

void concurr_get(struct range *_range) {
  size_t key;
  uint64_t *array = workload[_range->index];

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
    key = i;
    eh->Get(key);
    /*
    if (eh->Get(key) == NONE)
    {
            std::cout<<"Search the key "<< i << ": ERROR!"<<std::endl;
            eh->FindAnyway(key);
    }*/
  }
}

void concurr_delete(struct range *_range) {
  size_t key;
  for (uint64_t i = _range->begin; i < _range->end; ++i) {
    key = i;
    // eh->Delete(key);

    if (eh->Delete(key) == false) {
      std::cout << "Delete the key " << i << ": ERROR!" << std::endl;
    }
  }
}

int main(int argc, char const *argv[]) {
  assert(argc >= 4);
  int initCap = atoi(argv[1]);
  int insert_num = atoi(argv[2]);
  int thread_num = atoi(argv[3]);

  Allocator::Initialize(pool_name, pool_size);

  std::cout << "The initCap is " << initCap << std::endl;
  std::cout << "The inserted number is " << insert_num << std::endl;
  std::cout << "The thread number is " << thread_num << std::endl;

  double duration;

  eh = reinterpret_cast<Finger_EH *>(Allocator::GetRoot(sizeof(Finger_EH)));
  new (eh) Finger_EH(initCap);

  workload = new uint64_t *[thread_num];

  int chunk_size = insert_num / thread_num;
  struct range rarray[thread_num];
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].index = i;
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

  int i;
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL},
                     length = 4;
  init_by_array64(init, length);

  /* Genereate Workload*/
  for (int i = 0; i < thread_num; ++i) {
    auto num = rarray[i].end - rarray[i].begin;
    workload[i] = new uint64_t[num];
    for (int j = 0; j < num; ++j) {
      workload[i][j] = genrand64_int64();
    }
  }

  /*-----------------------------------------------Concurrent Insertion
   * Test-----------------------------------------------------------------------*/
  std::thread *thread_array[thread_num];
  /* The muli-thread execution begins*/
  eh->getNumber();

  LOG("Concurrent insertion "
      "begin-----------------------------------------------------------------");
  // System::profile("Insertion", [&](){
  gettimeofday(&tv1, NULL);
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(concurr_insert, &rarray[i]);
  }
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
  gettimeofday(&tv2, NULL);
  duration = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double)(tv2.tv_sec - tv1.tv_sec);
  printf(
      "For %d threads, Insert Total time = %f seconds, the throughput is %f "
      "options/s\n",
      thread_num,
      (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double)(tv2.tv_sec - tv1.tv_sec),
      insert_num / duration);
  //});

  // printf("the collison check is %d\n", eh->count);
  LOG("Concurrent insertion "
      "end------------------------------------------------------------------");
  eh->getNumber();
  // eh->CheckDepthCount();
  return 0;
}
