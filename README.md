# Dash: Scalable Hashing on Persistent Memory

Persistent memory friendly hashing index.

## What's included

- Dash EH - Proposed DASH extendable hashing
- Dash LH - Proposed DASH linear Hashing
- CCEH - PMDK patched CCEH variant used in our benchmark
- Level Hashing - PMDK patched level hashing variant used in our benchmark
- Mini benchmark framework
- All benchmark results

Fully open-sourced under MIT license.


## Build using CMake

We tested our build with Linux Kernel 5.2.11 and gcc 9.2.0.

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON .. 
make -j
```

## Running

```bash
./build/test_pmem --helpshort
Usage: 
    ./build/test_pmem [OPTION...]

-index      the index to evaluate:dash-ex/dash-lh/cceh/level (default: "dash-ex")
-op         the type of operation to execute:insert/pos/neg/delete/mixed (default: "full")
-n          the number of warm-up workload (default: 0)
-p          the number of operations(insert/search/delete) to execute (default: 20000000)
-t          the number of concurrent threads (default: 1)
-r          search ratio for mixed workload: 0.0~1.0 (default: 1.0)
-s          insert ratio for mixed workload: 0.0~1.0 (default: 0.0)
-d          delete ratio for mixed workload: 0.0~1.0 (default: 0.0)
-e          whether to register epoch in application level: 0/1 (default: 0)
-k          (the type of stored keys: fixed/variable (default: "fixed")
-vl         the length of the variable length key (default: 16)
```

## Miscellaneous

We noticed a possible `mmap` bug on our testing environment: `MAP_SHARED_VALIDATE` is incompatible with `MAP_FIXED_NOREPLACE`.
To ensure safe memory mapping, we modified the original PMDK to use `MAP_SHARED` rather than `MAP_SHARED_VALIDATE`, which has the same functionality as the former one except for extra flag validation.
For a more detailed explanation and minimal reproducible code, please check out our [blog post](https://blog.haoxp.xyz/posts/mmap-bug/) about this issue.
