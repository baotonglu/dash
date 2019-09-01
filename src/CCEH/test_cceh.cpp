#include <sys/time.h>
#include <thread>
#include "../../util/System.hpp"
#include "../../util/random.h"
#include "../allocator.h"
#include "CCEH_baseline.h"
#include "libpmemobj.h"
#include "../utils.h"

#define FIXED 1

//static const char *pool_name = "/mnt/pmem0/pmem_cceh.data";
static const char *pool_name = "pmem_cceh.data";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 10ul;

#ifdef FIXED
CCEH<Key_t> *eh;
#else
CCEH<char *> *eh;
#endif

uint64_t *workload;
struct timeval tv1, tv2;

struct range {
  uint64_t index;
  uint64_t begin;
  uint64_t end;
};

/*fixed length 16-byte key*/
struct string_key{
  char key[16];
};

void concurr_insert(struct range *_range) {
#ifdef FIXED
  size_t key;
#else
  char* key;
  string_key *var_workload = reinterpret_cast<string_key *>(workload);
#endif
  char arr[64];
  Value_t value = (Value_t)arr;

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
#else
    key = reinterpret_cast<char *>(var_workload + i);
#endif
    eh->Insert(key, value);
  }
}

void concurr_get(struct range *_range) {
#ifdef FIXED
  size_t key;
#else
  char* key;
  string_key *var_workload = reinterpret_cast<string_key *>(workload);
#endif
  
  uint32_t not_found = 0;
  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
#else
    key = reinterpret_cast<char *>(var_workload + i);
#endif
    if (eh->Get(key) == NONE)
    {
      not_found++;
    }
  }
  std::cout <<"not_found = "<<not_found<<std::endl;
}

void concurr_delete(struct range *_range) {
#ifdef FIXED
  size_t key;
#else
  char* key;
  string_key *var_workload = reinterpret_cast<string_key *>(workload);
#endif

  uint32_t not_found = 0;
  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
#else
    key = reinterpret_cast<char *>(var_workload + i);
#endif
    if (eh->Delete(key) == false) {
	    not_found++;
    } 
  }
  std::cout<<"not found = "<<not_found<<std::endl;
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

#ifdef FIXED
  eh = reinterpret_cast<CCEH<Key_t> *>(Allocator::GetRoot(sizeof(CCEH<Key_t>)));
  new (eh) CCEH<Key_t>(initCap);
#else
  eh = reinterpret_cast<CCEH<char *> *>(Allocator::GetRoot(sizeof(CCEH<char *>)));
  new (eh) CCEH<char *>(initCap);
#endif
  eh->pool_addr = Allocator::Get()->pm_pool_;

  int chunk_size = insert_num / thread_num;
  struct range rarray[thread_num];
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].index = i;
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

 /* Generate Workload for fixed_length key or variable_length key*/
  Allocator::ZAllocate((void **)&workload, kCacheLineSize, sizeof(uint64_t) * (insert_num + 100) * 4);
  int i;
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL},
  length = 4;
  init_by_array64(init, length);

  /* Generate Workload*/
  char var_key[24];
  for(int i = 0; i < insert_num*2 + 2; ++i){
    uint64_t _key = genrand64_int64();
#ifdef FIXED
    workload[i] = _key;
#else
    snprintf(var_key, 24, "%lld", _key);
    var_key[15] = '\0';
    string_key *var_workload = reinterpret_cast<string_key *>(workload);
    strcpy(reinterpret_cast<char *>(var_workload + i), var_key);
#endif
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
  /*-----------------------------------------------Concurrent Get
   * Test-----------------------------------------------------------------------*/
  Allocator::ReInitialize_test_only(pool_name, pool_size);
  LOG("Concurrent positive get "
      "begin!------------------------------------------------------------");
  // System::profile("NP_search", [&](){
  gettimeofday(&tv1, NULL);
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(concurr_get, &rarray[i]);
  }

  for (int i = 0; i < thread_num; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
  gettimeofday(&tv2, NULL);
  duration = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double)(tv2.tv_sec - tv1.tv_sec);
  printf(
      "For %d threads, Get Total time = %f seconds, the throughput is %f "
      "options/s\n",
      thread_num,
      (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double)(tv2.tv_sec - tv1.tv_sec),
      insert_num / duration);
  //});
  LOG("Concurrent positive get "
      "end!---------------------------------------------------------------");

  LOG("Concurrent negative get "
      "begin!-------------------------------------------------------------");
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = insert_num + i * chunk_size + 1;
    rarray[i].end = insert_num + (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + insert_num + 1;

  gettimeofday(&tv1, NULL);
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(concurr_get, &rarray[i]);
  }

  for (int i = 0; i < thread_num; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
  gettimeofday(&tv2, NULL);
  duration = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double)(tv2.tv_sec - tv1.tv_sec);
  printf(
      "For %d threads, Get Total time = %f seconds, the throughput is %f "
      "options/s\n",
      thread_num,
      (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double)(tv2.tv_sec - tv1.tv_sec),
      insert_num / duration);
  LOG("Concurrent negative get "
      "end!---------------------------------------------------------------");

  LOG("Concurrent deletion begin-----------------------------------------------------------------");
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

  //System::profile("Delete", [&](){
  gettimeofday(&tv1, NULL);
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(concurr_delete, &rarray[i]);
  }

  for (int i = 0; i < thread_num; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
  gettimeofday(&tv2, NULL);
 //});
  duration = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double)(tv2.tv_sec - tv1.tv_sec);
  printf(
      "For %d threads, Delete Total time = %f seconds, the throughput is %f "
      "options/s\n",
      thread_num,
      (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double)(tv2.tv_sec - tv1.tv_sec),
      insert_num / duration);
  LOG("Concurrent deletion end-------------------------------------------------------------------");
  return 0;
}
