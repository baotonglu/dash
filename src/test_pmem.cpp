#include <gflags/gflags.h>
#include <immintrin.h>
#include <sys/time.h>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <thread>
#include "../util/System.hpp"
#include "../util/random.h"
#include "../util/uniform.hpp"
#include "allocator.h"
#include "ex_finger.h"
#include "libpmemobj.h"
#include "utils.h"

static const char *pool_name = "/mnt/pmem0/pmem_hash.data";
// static const char *pool_name = "pmem_hash.data";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;
DEFINE_string(k, "fixed", "the type of stored keys: fixed/variable");
DEFINE_int32(i, 64, "the initial number of segments in extendible hashing");
DEFINE_int32(t, 1, "the number of concurrent threads");
DEFINE_int32(n, 0, "the number of load");
DEFINE_int32(p, 20000000,
             "the number of operations(insert/search/deletion) to execute");
DEFINE_string(op, "full",
              "which type of operation to execute:insert/pos/neg/delete/mixed");
DEFINE_double(r, 1, "read ratio for mixed workload");
DEFINE_double(s, 0, "insert ratio for mixed workload");
DEFINE_double(d, 0, "delete ratio for mixed workload");

// ADD and SUB return the value after add or sub
#define ADD(_p, _v) (__atomic_add_fetch(_p, _v, __ATOMIC_SEQ_CST))
#define SUB(_p, _v) (__atomic_sub_fetch(_p, _v, __ATOMIC_SEQ_CST))
#define LOAD(_p) (__atomic_load_n(_p, __ATOMIC_SEQ_CST))
#define STORE(_p, _v) (__atomic_store_n(_p, _v, __ATOMIC_SEQ_CST))

int32_t initCap, thread_num, load_num, operation_num;
std::string operation;
std::string key_type;
int bar_a, bar_b, bar_c;
double read_ratio, insert_ratio, delete_ratio;
std::mutex mtx;
std::condition_variable cv;
bool finished = false;
struct timeval tv1, tv2;

struct range {
  int index;
  int begin;
  int end;
  int length; /*if this is the variable lenght key, use this parameter to
                 indicate the length of the key*/
  void *workload;
  uint64_t random_num;
};

void set_affinity(uint32_t idx) {
  cpu_set_t my_set;
  CPU_ZERO(&my_set);
  if (idx < 24) {
    CPU_SET(idx, &my_set);
  } else {
    CPU_SET(idx + 24, &my_set);
  }
  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);
}

template <class T>
Finger_EH<T> *InitializeIndex(int seg_num) {
  Finger_EH<T> *eh = reinterpret_cast<Finger_EH<T> *>(
      Allocator::GetRoot(sizeof(Finger_EH<T>)));
  new (eh) Finger_EH<T>(seg_num);
  return eh;
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
void generate_16B(void *memory_region, int generate_num, int length,
                  bool persist) {
  string_key *var_key;
  uint64_t *_key;
  char *workload = reinterpret_cast<char *>(memory_region);

  int word_num = (length / 8) + (((length % 8) != 0) ? 1 : 0);
  _key = reinterpret_cast<uint64_t *>(malloc(word_num * sizeof(uint64_t)));

  for (int i = 0; i < generate_num; ++i) {
    var_key = reinterpret_cast<string_key *>(workload +
                                             i * (length + sizeof(string_key)));
    var_key->length = length;
    for (int i = 0; i < word_num; ++i) {
      _key[i] = genrand64_int64();
    }
    memcpy(var_key->key, _key, length);
  }

  if (persist) {
    Allocator::Persist(memory_region,
                       generate_num * (sizeof(string_key) + length));
  }
}

template <class T>
void Load(int kv_num, Finger_EH<T> *index, int length) {
  if (kv_num == 0) return;
  std::string fixed("fixed");
  T key;
  if constexpr (!std::is_pointer_v<T>) {
    for (int i = 0; i < kv_num; ++i) {
      key = genrand64_int64();
      index->Insert(key, DEFAULT);
    }
  } else { /*genereate 16B key*/
    char *persist_workload;
    Allocator::ZAllocate((void **)&persist_workload, kCacheLineSize,
                         (sizeof(string_key) + length) * kv_num);
    generate_16B(persist_workload, kv_num, length, true);

    int string_key_size = sizeof(string_key) + length;
    for (int i = 0; i < kv_num; ++i) {
      key = reinterpret_cast<T>(persist_workload + i * string_key_size);
      index->Insert(key, DEFAULT);
    }
  }
  std::cout << "Finish loading " << kv_num << " keys" << std::endl;
}

template <class T>
void concurr_insert(struct range *_range, Finger_EH<T> *index) {
  set_affinity(_range->index);
  int begin = _range->begin;
  int end = _range->end;
  char *workload = reinterpret_cast<char *>(_range->workload);
  T key;

  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/

  // if(key_type == "fixed"){/*fixed length key benchmark*/
  if constexpr (!std::is_pointer_v<T>) {
    T *key_array = reinterpret_cast<T *>(workload);
    for (int i = begin; i < end; ++i) {
      index->Insert(key_array[i], DEFAULT);
    }
  } else {
    T var_key;
    int string_key_size = sizeof(string_key) + _range->length;
    for (int i = begin; i < end; ++i) {
      var_key = reinterpret_cast<T>(workload + string_key_size * i);
      index->Insert(var_key, DEFAULT);
    }
  }

  if (SUB(&bar_c, 1) == 0) {
    std::unique_lock<std::mutex> lck(mtx);
    finished = true;
    cv.notify_one();
  }
}

template <class T>
void mixed(struct range *_range, Finger_EH<T> *index) {
  set_affinity(_range->index);
  int begin = _range->begin;
  int end = _range->end;
  char *workload = reinterpret_cast<char *>(_range->workload);
  T *key_array = reinterpret_cast<T *>(_range->workload);
  T key;
  int string_key_size = sizeof(string_key) + _range->length;

  UniformRandom rng(_range->random_num);
  uint32_t random;
  uint32_t not_found = 0;

  uint32_t insert_sign = (uint32_t)(insert_ratio * 100);
  uint32_t read_sign = (uint32_t)(read_ratio * 100) + insert_sign;
  uint32_t delete_sign = (uint32_t)(delete_ratio * 100) + read_sign;
  std::cout << "insert sign = " << insert_sign << std::endl;
  std::cout << "read sign = " << read_sign << std::endl;

  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/

  for (uint64_t i = begin; i < end; ++i) {
    if constexpr (std::is_pointer_v<T>) { /* variable length*/
      key = reinterpret_cast<T>(workload + string_key_size * i);
    } else {
      key = key_array[i];
    }

    random = rng.next_uint32() % 100;
    if (random < insert_sign) { /*insert*/
      index->Insert(key, DEFAULT);
    } else if (random < read_sign) { /*get*/
      if (index->Get(key) == NONE) {
        not_found++;
      }
    } else { /*delete*/
      index->Delete(key);
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

template <class T>
void GeneralBench(range *rarray, Finger_EH<T> *index, int thread_num,
                  int operation_num, std::string profile_name,
                  void (*test_func)(struct range *, Finger_EH<T> *)) {
  std::thread *thread_array[1024];
  profile_name = profile_name + std::to_string(thread_num);
  double duration;
  finished = false;
  bar_a = 1;
  bar_b = thread_num;
  bar_c = thread_num;

  // System::profile(profile_name, [&]() {
  for (int i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(*test_func, &rarray[i], index);
  }

  while (LOAD(&bar_b) != 0)
    ;                                     // Spin
  std::unique_lock<std::mutex> lck(mtx);  // get the lock of condition variable

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
      "%d threads, Time = %f s, throughput = %f "
      "ops/s\n",
      thread_num,
      (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double)(tv2.tv_sec - tv1.tv_sec),
      operation_num / duration);
  //});
}

void *GenerateWorkload(int generate_num, int length) {
  /*Since there are both positive search and negative search, it should generate
   * 2 * generate_num workload*/
  void *workload;
  if (key_type == "fixed") {
    workload = malloc(generate_num * sizeof(uint64_t));
    generate_8B(workload, generate_num, false);
  } else { /*Generate the variable lengh workload*/
    workload = malloc(generate_num * (length + sizeof(string_key)));
    generate_16B(workload, generate_num, length, false);
  }
  return workload;
}

template <class T>
void Run() {
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL},
                     length = 4;
  init_by_array64(init, length); /*initialize random number generation*/
  Allocator::Initialize(pool_name, pool_size); /* allocator initialization*/

  /* Initialize Index*/
  Finger_EH<T> *index = InitializeIndex<T>(initCap);

  /* Load Phase, first generate workload, if the key is variable length key, we
   * use 16B key to represent*/
  Load<T>(load_num, index, 16);

  /* Generate the workload and corresponding range array*/
  void *workload = GenerateWorkload(operation_num * 2, 16);
  void *insert_workload;
  if (key_type == "variable") {
    Allocator::ZAllocate((void **)&insert_workload, kCacheLineSize,
                         (sizeof(string_key) + 16) * operation_num * 2);
    memcpy(insert_workload, workload,
           (sizeof(string_key) + 16) * operation_num * 2);
  } else {
    insert_workload = workload;
  }

  /* Description of the workload*/
  srand((unsigned)time(NULL));
  struct range *rarray;
  int chunk_size = operation_num / thread_num;
  rarray = reinterpret_cast<range *>(malloc(thread_num * (sizeof(range))));
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].index = i;
    rarray[i].random_num = rand();
    rarray[i].begin = i * chunk_size;
    rarray[i].end = (i + 1) * chunk_size;
    rarray[i].length = 16;
    rarray[i].workload = workload;
  }
  rarray[thread_num - 1].end = operation_num;

  /* Benchmark Phase */
  if (operation == "insert") {
    for (int i = 0; i < thread_num; ++i) {
      rarray[i].workload = insert_workload;
    }
    GeneralBench<T>(rarray, index, thread_num, operation_num, "Insert",
                    &concurr_insert);
  } else if (operation == "pos") {
  } else if (operation == "neg") {
  } else if (operation == "delete") {
  } else if (operation == "mixed") { /*mixed workload, the default setting is
                                        read 80%, insert 20%*/
    for (int i = 0; i < thread_num; ++i) {
      rarray[i].workload = insert_workload;
    }
    GeneralBench<T>(rarray, index, thread_num, operation_num, "Insert", &mixed);
  } else { /*do the benchmark for all single operations*/
    for (int i = 0; i < thread_num; ++i) {
      rarray[i].workload = insert_workload;
    }
    GeneralBench<T>(rarray, index, thread_num, operation_num, "Insert",
                    &concurr_insert);
    for (int i = 0; i < thread_num; ++i) {
      rarray[i].workload = workload;
    }
  }
}

bool check_ratio() {
  int read_portion = (int)(read_ratio * 100);
  int insert_portion = (int)(insert_ratio * 100);
  int delete_portion = (int)(delete_ratio * 100);
  if ((read_portion + insert_portion + delete_portion) != 100) return false;
  return true;
}

int main(int argc, char *argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  initCap = FLAGS_i;
  thread_num = FLAGS_t;
  load_num = FLAGS_n;
  operation_num = FLAGS_p;
  key_type = FLAGS_k;
  std::string fixed("fixed");
  operation = FLAGS_op;

  read_ratio = FLAGS_r;
  insert_ratio = FLAGS_s;
  delete_ratio = FLAGS_d;
  if (!check_ratio()) {
    std::cout << "The ratio is wrong!" << std::endl;
    return 0;
  }

  if (key_type.compare(fixed) == 0) {
    Run<uint64_t>();
  } else {
    Run<string_key *>();
  }
}