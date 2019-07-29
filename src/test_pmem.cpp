#include "ex_finger.h"
#include "libpmemobj.h"
#include "allocator.h"
#include "utils.h"

static const char *pool_name = "pmem_hash.data";
static const uint32_t pool_size = 1024 * 1024 * 1024;

int main(int argc, char const *argv[]) {
  assert(argc >= 4);
  int initCap = atoi(argv[1]);
  int insert_num = atoi(argv[2]);
  int thread_num = atoi(argv[3]);

  Allocator::Initialize(pool_name, pool_size);

  std::cout << "The initCap is " << initCap << std::endl;
  std::cout << "The inserted number is " << insert_num << std::endl;
  std::cout << "The thread number is " << thread_num << std::endl;
  auto root_finger_eh =
      reinterpret_cast<Finger_EH *>(Allocator::GetRoot(sizeof(Finger_EH)));
  new (root_finger_eh) Finger_EH(initCap);
  return 0;
}
