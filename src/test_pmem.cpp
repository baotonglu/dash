#include <immintrin.h>
#include <sys/time.h>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <thread>
#include "../util/System.hpp"
#include "../util/random.h"
#include "../util/uniform.hpp"
#include "allocator.h"
#include "libpmemobj.h"
#include "utils.h"

//#define LINEAR 1
//#define FIXED 1
//#define MIXED_TEST 1
//#define TEST_BANDWIDTH 1

#ifndef LINEAR
#include "ex_finger.h"
#else
#include "lh_finger.h"
#endif

static const char *pool_name = "/mnt/pmem0/pmem_hash.data";
// static const char *pool_name = "pmem_hash.data";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;
#define BILLION 1E9

#ifndef LINEAR

#ifdef FIXED
Finger_EH<Key_t> *eh;
#else
Finger_EH<char *> *eh;
#endif

#else

#ifdef FIXED
Linear<Key_t> *eh;
#else
Linear<char *> *eh;
#endif

#endif
uint64_t *workload;
uint64_t *persist_workload;
uint64_t *value_workload;
struct timeval tv1, tv2;
clock_t start, finish;
struct timespec requestStart, requestEnd;
int insert_num;
int bar_a;
int bar_b;
int bar_c;  // use to notify the
std::mutex mtx;
std::condition_variable cv;
bool finished = false;
// ADD and SUB return the value after add or sub
#define ADD(_p, _v) (__atomic_add_fetch(_p, _v, __ATOMIC_SEQ_CST))
#define SUB(_p, _v) (__atomic_sub_fetch(_p, _v, __ATOMIC_SEQ_CST))
#define LOAD(_p) (__atomic_load_n(_p, __ATOMIC_SEQ_CST))
#define STORE(_p, _v) (__atomic_store_n(_p, _v, __ATOMIC_SEQ_CST))

/*fixed length 16-byte key*/
// struct string_key{
//  char key[16];
//  int length;
//};

struct range {
  uint32_t index;
  uint64_t begin;
  uint64_t end;
  uint64_t random_num;
};

void clear_cache(int insert_num) {
  uint32_t not_found = 0;
  auto _value_workload = reinterpret_cast<Value_t *>(value_workload);
  for (int i = 0; i < insert_num; ++i) {
    if (_value_workload[i] == NONE) {
      not_found++;
    }
  }
  printf("clear cache: %u\n", not_found);
}

void mixed(struct range *_range) {
  cpu_set_t my_set;
  CPU_ZERO(&my_set);
  CPU_SET(_range->index, &my_set);
  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);

#ifdef FIXED
  size_t key;
#else
  char *key;
  // string_key *var_workload = reinterpret_cast<string_key *>(workload);
  string_key *var_workload = reinterpret_cast<string_key *>(persist_workload);
#endif
  UniformRandom rng(_range->random_num);
  char arr[64];
  Value_t value;
  uint32_t random;
  uint32_t not_found = 0;
  auto _value_workload = reinterpret_cast<Value_t *>(value_workload);

  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
#else
    // key = reinterpret_cast<char *>(var_workload + i);
    key = var_workload[i].key;
#endif
    random = rng.next_uint32() % 10;
    if (random <= 1) {
      value = _value_workload[i];
      eh->Insert(key, value);
    } else {
      if (eh->Get(key) == NONE) {
        not_found++;
      }
    }
  }
  std::cout << "not_found = " << not_found << std::endl;
  /*the last thread notify the main thread to wake up*/
  if (SUB(&bar_c, 1) == 0) {
    std::unique_lock<std::mutex> lck(mtx);
    finished = true;
    cv.notify_one();
  }
}

void concurr_insert(struct range *_range) {
  cpu_set_t my_set;
  CPU_ZERO(&my_set);
  CPU_SET(_range->index, &my_set);
  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);

#ifdef FIXED
  size_t key;
#else
  char *key;
  string_key *var_workload = reinterpret_cast<string_key *>(persist_workload);
  // string_key *var_workload = reinterpret_cast<string_key *>(workload);
#endif
  char arr[64];
  Value_t value;
  auto _value_workload = reinterpret_cast<Value_t *>(value_workload);

  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
    // key = i;
#else
    // key = reinterpret_cast<char *>(var_workload + i);
    key = var_workload[i].key;
#endif
    value = _value_workload[i];
    eh->Insert(key, value);
  }

  /*the last thread notify the main thread to wake up*/
  if (SUB(&bar_c, 1) == 0) {
    std::unique_lock<std::mutex> lck(mtx);
    finished = true;
    cv.notify_one();
  }
}

void concurr_get(struct range *_range) {
  cpu_set_t my_set;
  CPU_ZERO(&my_set);
  CPU_SET(_range->index, &my_set);
  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);

#ifdef FIXED
  size_t key;
#else
  char *key;
  string_key *var_workload = reinterpret_cast<string_key *>(workload);
#endif
  uint32_t not_found = 0;
  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/

  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
#else
    // key = reinterpret_cast<char *>(var_workload + i);
    key = (char *)&(var_workload[i].key);
#endif
    if (eh->Get(key) == NONE) {
      not_found++;
    }
  }
  std::cout << "not_found = " << not_found << std::endl;
  /*the last thread notify the main thread to wake up*/
  if (SUB(&bar_c, 1) == 0) {
    std::unique_lock<std::mutex> lck(mtx);
    finished = true;
    cv.notify_one();
  }
}

void concurr_delete(struct range *_range) {
  cpu_set_t my_set;
  CPU_ZERO(&my_set);
  CPU_SET(_range->index, &my_set);
  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);

#ifdef FIXED
  size_t key;
#else
  char *key;
  string_key *var_workload = reinterpret_cast<string_key *>(workload);
#endif
  uint32_t not_found = 0;

  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/
  for (uint64_t i = _range->begin; i < _range->end; ++i) {
#ifdef FIXED
    key = workload[i];
#else
    // key = reinterpret_cast<char *>(var_workload + i);
    key = var_workload[i].key;
#endif
    if (eh->Delete(key) == false) {
      not_found++;
    }
  }
  std::cout << "not found = " << not_found << std::endl;
  /*the last thread notify the main thread to wake up*/
  if (SUB(&bar_c, 1) == 0) {
    std::unique_lock<std::mutex> lck(mtx);
    finished = true;
    cv.notify_one();
  }
}

/*Initialize the PM space for the tested data structure*/
void initialize_index(int initCap) {
#ifndef LINEAR

#ifdef FIXED
  eh = reinterpret_cast<Finger_EH<Key_t> *>(
      Allocator::GetRoot(sizeof(Finger_EH<Key_t>)));
  new (eh) Finger_EH<Key_t>(initCap);
#else
  eh = reinterpret_cast<Finger_EH<char *> *>(
      Allocator::GetRoot(sizeof(Finger_EH<char *>)));
  new (eh) Finger_EH<char *>(initCap);
#endif

#else

#ifdef FIXED
  eh = reinterpret_cast<Linear<Key_t> *>(
      Allocator::GetRoot(sizeof(Linear<Key_t>)));
  new (eh) Linear<Key_t>();
#else
  eh = reinterpret_cast<Linear<char *> *>(
      Allocator::GetRoot(sizeof(Linear<char *>)));
  new (eh) Linear<char *>();
#endif

#endif
  eh->pool_addr = Allocator::Get()->pm_pool_;
}

/*generate random 8-byte number and store it in the memory_region*/
void generate_8B(void *memory_region, int generate_num, bool persist) {
  uint64_t *array = reinterpret_cast<uint64_t *>(memory_region);

  for (int i = 0; i < generate_num; ++i) {
    array[i] = genrand64_int64();
  }

  if (persist) {
    Allocator::Persist(memory_region, generate_num * sizeof(uint64_t));
  }
}

/*generate random 16-byte string and store it in the memory_region*/
void generate_16B(void *memory_region, int generate_num, bool persist) {
  string_key *var_workload = reinterpret_cast<string_key *>(memory_region);
  char var_key[24];

  for (int i = 0; i < generate_num; ++i) {
    // uint64_t _key = genrand64_int64();
    // snprintf(var_key, 24, "%lld", _key);
    // memcpy(var_workload + i, var_key, 16);
    uint64_t *key = reinterpret_cast<uint64_t *>(var_workload[i].key);
    key[0] = genrand64_int64();
    key[1] = genrand64_int64();
    var_workload[i].length = 16;
  }

  if (persist) {
    Allocator::Persist(memory_region, generate_num * sizeof(string_key));
  }
}

void generalBench(range *rarray, int thread_num, std::string profile_name,
                  void (*test_func)(struct range *)) {
  std::thread *thread_array[1024];
  profile_name = profile_name + std::to_string(thread_num);
  double duration;
  finished = false;
  clear_cache(thread_num);
  bar_a = 1;
  bar_b = thread_num;
  bar_c = thread_num;
  System::profile(profile_name, [&]() {
    for (int i = 0; i < thread_num; ++i) {
      thread_array[i] = new std::thread(*test_func, &rarray[i]);
    }

    while (LOAD(&bar_b) != 0)
      ;  // Spin
    std::unique_lock<std::mutex> lck(
        mtx);  // get the lock of condition variable
    // printf("All threads are ready\n");

    gettimeofday(&tv1, NULL);
    STORE(&bar_a, 0);  // start test
    while (!finished)
      cv.wait(lck);  // go to sleep and wait for the wake-up from child threads
    gettimeofday(&tv2, NULL);  // test end

    for (int i = 0; i < thread_num; ++i) {
      thread_array[i]->join();
      delete thread_array[i];
    }
    duration = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
               (double)(tv2.tv_sec - tv1.tv_sec);
    printf(
        "For %d threads,Total time = %f seconds, the throughput is %f "
        "options/s\n",
        thread_num,
        (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
            (double)(tv2.tv_sec - tv1.tv_sec),
        insert_num / duration);
  });
}

int main(int argc, char const *argv[]) {
  assert(argc >= 4);
  int initCap = atoi(argv[1]);
  insert_num = atoi(argv[2]);
  int thread_num = atoi(argv[3]);
  int mixed_num = 200000000;
  int i;
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL},
                     length = 4;
  init_by_array64(init, length); /*initialize random number generation*/

  std::cout << "The string key length = " << sizeof(struct string_key)
            << std::endl;

  Allocator::Initialize(pool_name, pool_size);

  std::cout << "The initCap is " << initCap << std::endl;
  std::cout << "The inserted number is " << insert_num << std::endl;
  std::cout << "The thread number is " << thread_num << std::endl;

  srand((unsigned)time(NULL));
  initialize_index(initCap);
  /*************************************************Generate
   * Workload************************************************************/
  /* Description of the workload*/
  int chunk_size = insert_num / thread_num;
  struct range rarray[thread_num];
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].index = i;
    rarray[i].random_num = rand();
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

#ifdef TEST_BANDWIDTH
  struct range rarray_insert[24];
  int insert_chunk = insert_num / 24;
  for (int i = 0; i < 24; ++i) {
    rarray_insert[i].index = i;
    rarray_insert[i].random_num = rand();
    rarray_insert[i].begin = i * insert_chunk + 1;
    rarray_insert[i].end = (i + 1) * insert_chunk + 1;
  }
  rarray_insert[23].end = insert_num + 1;
#endif
  /* Generate Workload for fixed_length key or variable_length key*/
  int generate_num = insert_num;
#ifdef MIXED_TEST
  generate_num = mixed_num;
#endif

#ifdef FIXED
  workload = (uint64_t *)malloc((generate_num + 100) * sizeof(uint64_t) * 2);
#else
  workload = (uint64_t *)malloc((generate_num + 100) * sizeof(string_key) * 2);
#endif
  value_workload = (uint64_t *)malloc((generate_num + 100) * sizeof(uint64_t));

#ifndef FIXED
#ifdef MIXED_TEST
  Allocator::ZAllocate((void **)&persist_workload, kCacheLineSize,
                       sizeof(string_key) * (generate_num + 100) * 2);
#else
  Allocator::ZAllocate((void **)&persist_workload, kCacheLineSize,
                       sizeof(string_key) * (generate_num + 100));
#endif
#endif

#ifdef FIXED
  generate_8B(workload, generate_num * 2 + 2, false);
#else
  generate_16B(workload, generate_num * 2 + 2, false);
#endif

  generate_8B(value_workload, generate_num + 1, false);

#ifndef FIXED
  memcpy(persist_workload, workload, (generate_num + 1) * sizeof(string_key));
  Allocator::Persist(persist_workload, (generate_num + 1) * sizeof(string_key));
#ifdef MIXED_TEST
  string_key *_persist = reinterpret_cast<string_key *>(persist_workload);
  generate_16B(&_persist[generate_num + 1], generate_num + 1, true);
#endif
  // string_key *var_workload = reinterpret_cast<string_key *>(workload);
  // string_key *p_var_workload = reinterpret_cast<string_key
  // *>(persist_workload); for(int i = 0; i < generate_num + 1; ++i){
  //  strcpy(reinterpret_cast<char *>(p_var_workload + i), reinterpret_cast<char
  //  *>(var_workload + i));
  //}
#endif

  /**************************************************Benchmark***********************************************************/

  /******************Benchmark for insert***********************/
  printf("Insert workload begin\n");
#ifdef TEST_BANDWIDTH
  generalBench(rarray_insert, 24, "Insertion_", &concurr_insert);
#else
  /* normal benchmark for insert operation*/
  generalBench(rarray, thread_num, "Insertion_", &concurr_insert);
#endif

  /******************Benchmark for mixed workload************************/
#ifdef MIXED_TEST
  printf("Mixed workload begin\n");
  chunk_size = mixed_num / thread_num;
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = insert_num + i * chunk_size + 1;
    rarray[i].end = insert_num + (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + mixed_num + 1;
  generalBench(rarray, thread_num, "Mixed_", &mixed);
#endif

  /******************Benchmark for negative search***********************/

  printf("Neg search workload begin\n");

  chunk_size = insert_num / thread_num;
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = insert_num + i * chunk_size + 1;
    rarray[i].end = insert_num + (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + insert_num + 1;
  generalBench(rarray, thread_num, "Neg_get_", &concurr_get);

  /******************Benchmark for positive search***********************/
  printf("Pos search workload begin\n");
  chunk_size = insert_num / thread_num;
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;

  generalBench(rarray, thread_num, "Pos_get_", &concurr_get);

  /*********************Benchmark for delete ****************************/
  /*
  printf("Delete workload begin\n");

  for (int i = 0; i < thread_num; ++i) {
    rarray[i].begin = i * chunk_size + 1;
    rarray[i].end = (i + 1) * chunk_size + 1;
  }
  rarray[thread_num - 1].end = insert_num + 1;
  generalBench(rarray, thread_num, "Delete_", &concurr_delete);*/
  return 0;
}
