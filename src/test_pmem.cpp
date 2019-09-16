#include <sys/time.h>
#include <thread>
#include <cstring>
#include "../util/System.hpp"
#include "../util/random.h"
#include "allocator.h"
#include "libpmemobj.h"
#include "utils.h"

//#define LINEAR 1
#define FIXED 1
//#define TEST_BANDWIDTH 1

#ifndef LINEAR
#include "ex_finger.h"
#else
#include "lh_finger.h"
#endif

static const char *pool_name = "/mnt/pmem0/pmem_hash.data";
//static const char *pool_name = "pmem_hash.data";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;

#ifndef LINEAR

#ifdef FIXED
Finger_EH<Key_t> *eh;
#else
Finger_EH<char*> *eh;
#endif

#else

#ifdef FIXED
Linear<Key_t> *eh;
#else
Linear<char*> *eh;
#endif

#endif
uint64_t *workload;
uint64_t *persist_workload;
uint64_t *value_workload;
struct timeval tv1, tv2;

/*fixed length 16-byte key*/
struct string_key{
  char key[16];
};

struct range {
  uint64_t index;
  uint64_t begin;
  uint64_t end;
};

void concurr_insert(struct range *_range) {
#ifdef FIXED
  size_t key;
#else
  char* key;
  //string_key *var_workload = reinterpret_cast<string_key *>(workload);
  string_key *var_workload = reinterpret_cast<string_key *>(persist_workload);
#endif
  char arr[64];
  Value_t value;
  auto _value_workload = reinterpret_cast<Value_t *>(value_workload);

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
    //key = i;
#else
    key = reinterpret_cast<char *>(var_workload + i);
#endif
    value = _value_workload[i];
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
    //key = i;
#else
    key = reinterpret_cast<char *>(var_workload + i);
#endif
    if (eh->Get(key) == NONE)
    {
      not_found++;
    }
  }
  //std::cout <<"not_found = "<<not_found<<std::endl;
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
    //key = i;
#else
    key = reinterpret_cast<char *>(var_workload + i);
#endif
    if (eh->Delete(key) == false) {
	    not_found++;
    } 
  }
  //std::cout<<"not found = "<<not_found<<std::endl;
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
#ifndef LINEAR  

#ifdef FIXED
  eh = reinterpret_cast<Finger_EH<Key_t> *>(Allocator::GetRoot(sizeof(Finger_EH<Key_t>)));
  new (eh) Finger_EH<Key_t>(initCap);
#else
  eh = reinterpret_cast<Finger_EH<char *> *>(Allocator::GetRoot(sizeof(Finger_EH<char *>)));
  new (eh) Finger_EH<char *>(initCap);
#endif

#else

#ifdef FIXED
  eh = reinterpret_cast<Linear<Key_t> *>(Allocator::GetRoot(sizeof(Linear<Key_t>)));
  new (eh) Linear<Key_t>();
#else
  eh = reinterpret_cast<Linear<char *> *>(Allocator::GetRoot(sizeof(Linear<char *>)));
  new (eh) Linear<char *>();
#endif

#endif
  eh->pool_addr = Allocator::Get()->pm_pool_;

  /* Description of the workload*/
  int chunk_size = insert_num / thread_num;
  struct range rarray[thread_num];
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].index = i;
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

#ifdef TEST_BANDWIDTH
  struct range rarray_insert[24];
  int insert_chunk = insert_num / 24;
  for(int i = 0; i < 24; ++i){
    rarray_insert[i].index = i;
    rarray_insert[i].begin = i * insert_chunk + 1;
    rarray_insert[i].end = (i + 1) * insert_chunk + 1;
  }
  rarray_insert[23].end = insert_num + 1;
#endif
  /* Generate Workload for fixed_length key or variable_length key*/
  //Allocator::ZAllocate((void **)&workload, kCacheLineSize, sizeof(uint64_t) * (insert_num + 100) * 4);
  workload = (uint64_t*)malloc((insert_num + 100)*sizeof(uint64_t)*4);
  value_workload = (uint64_t*)malloc((insert_num + 100)*sizeof(uint64_t));
#ifndef FIXED 
  Allocator::ZAllocate((void **)&persist_workload, kCacheLineSize, sizeof(uint64_t) * (insert_num + 100) * 2);
#endif
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

  for(int i = 0; i < insert_num + 1; ++i){
    uint64_t _value = genrand64_int64();
    value_workload[i] = _value;
  }

#ifndef FIXED
  string_key *var_workload = reinterpret_cast<string_key *>(workload);
  string_key *p_var_workload = reinterpret_cast<string_key *>(persist_workload);
  for(int i = 0; i < insert_num + 1; ++i){
    strcpy(reinterpret_cast<char *>(p_var_workload + i), reinterpret_cast<char *>(var_workload + i));
  }
#endif
  /*-----------------------------------------------Concurrent Insertion
   * Test-----------------------------------------------------------------------*/
  std::thread *thread_array[1024];
  /* The muli-thread execution begins*/
  //eh->getNumber();

  LOG("Concurrent insertion "
      "begin-----------------------------------------------------------------");
std::string insertion = "Insertion_";
insertion = insertion + std::to_string(thread_num);
System::profile(insertion, [&](){
    gettimeofday(&tv1, NULL);
#ifdef TEST_BANDWIDTH
  for (int i = 0; i < 24; ++i) {
    thread_array[i] = new std::thread(concurr_insert, &rarray_insert[i]);
  }
  for (int i = 0; i < 24; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
#else
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(concurr_insert, &rarray[i]);
  }
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
#endif

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
  });

  // printf("the collison check is %d\n", eh->count);
  LOG("Concurrent insertion "
      "end------------------------------------------------------------------");
  //eh->getNumber();
  /*-----------------------------------------------Concurrent Get
   * Test-----------------------------------------------------------------------*/
 
  Allocator::ReInitialize_test_only(pool_name, pool_size);
  LOG("Concurrent positive get "
      "begin!------------------------------------------------------------");
  std::string pos = "Pos_search_";
 pos = pos + std::to_string(thread_num);
 System::profile(pos, [&](){
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
  });
  
  LOG("Concurrent negative get "
      "begin!-------------------------------------------------------------");
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = insert_num + i * chunk_size + 1;
    rarray[i].end = insert_num + (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + insert_num + 1;

  std::string neg= "NP_search_";
 neg = neg + std::to_string(thread_num);
  System::profile(neg, [&](){
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
  });
  LOG("Concurrent negative get "
      "end!---------------------------------------------------------------");
  
  /*-----------------------------------------------Concurrent Delete
   * Test--------------------------------------------------------------------*/
  
  LOG("Concurrent deletion begin-----------------------------------------------------------------");
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

  std::string del = "Delete_";
  del = del + std::to_string(thread_num);
  System::profile(del, [&](){
  gettimeofday(&tv1, NULL);
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(concurr_delete, &rarray[i]);
  }

  for (int i = 0; i < thread_num; ++i) {
    thread_array[i]->join();
    delete thread_array[i];
  }
  gettimeofday(&tv2, NULL);
 });
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
