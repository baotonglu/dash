#ifndef UNIFORM_HPP
#define UNIFORM_HPP
class UniformRandom {
 public:
  UniformRandom() : seed_(0) {}
  explicit UniformRandom(uint64_t seed) : seed_(seed) {}

  uint64_t get_current_seed() const { return seed_; }
  void set_current_seed(uint64_t seed) { seed_ = seed; }

  /**
   * @brief Fill up the give memory with random data.
   * @details
   * Call this to pre-calculate many random numbers. When the function we test
   * is very fast, generating random numbers might become the bottleneck. This
   * method is to avoid it. Use this like following:
   * @code{.cpp}
   * UniformRandom random(context->get_thread_id());  // or any other random
   * seed AlignedMemory memory;
   * CHECK_ERROR(context->get_thread_memory()->get_node_memory()->allocate_numa_memory(
   *      (1 << 16) * sizeof(uint32_t), &memory));
   * random.fill_memory(&memory);
   * const uint32_t *randoms = reinterpret_cast<const
   * uint32_t*>(memory.get_block()); <start the experiment timer> for (int i =
   * 0; i < TRIALS; ++i) { uint32_t value = randoms[i & 0xFFFF];
   *   ...
   * @endcode
   */

  uint64_t next_uint64() {
    return (static_cast<uint64_t>(next_uint32()) << 32) | next_uint32();
  }
  uint32_t next_uint32() {
    seed_ = seed_ * 0xD04C3175 + 0x53DA9022;
    return (seed_ >> 32) ^ (seed_ & 0xFFFFFFFF);
  }

 private:
  uint64_t seed_;

  /**
   * C is a run-time constant randomly chosen within [0 .. A] that can be
   * varied without altering performance. The same C value, per field
   * (C_LAST, C_ID, and OL_I_ID), must be used by all emulated terminals.
   * constexpr, but let's not bother C++11.
   */
  uint32_t get_c(uint32_t A) const {
    // yes, I'm lazy. but this satisfies the spec.
    const uint64_t kCSeed = 0x734b00c6d7d3bbdaULL;
    return kCSeed % (A + 1);
  }
};
#endif