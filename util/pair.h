#ifndef UTIL_PAIR_H_
#define UTIL_PAIR_H_

#include <cstdlib>
#include <immintrin.h>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

typedef size_t Key_t;
typedef const char* Value_t;

const Key_t SENTINEL = -2;  // 11111...110
const Key_t INVALID = -1;   // 11111...111

const Value_t NONE = 0x0;
const Value_t DEFAULT = reinterpret_cast<Value_t>(1);

/*variable length key*/
struct string_key{
    int length;
    char key[0];
};

struct Pair {
  Key_t key;
  Value_t value;

  Pair(void) : key{INVALID} {}

  Pair(Key_t _key, Value_t _value) : key{_key}, value{_value} {}

  Pair& operator=(const Pair& other) {
    key = other.key;
    value = other.value;
    return *this;
  }

  void* operator new(size_t size) {
    void* ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }

  void* operator new[](size_t size) {
    void* ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }
};

#endif  // UTIL_PAIR_H_
