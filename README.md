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
make test_pmem
```

### Run on a real device

```bash
./build/test_pmem 64 1000000 1
```

### Simulate on DRAM 

```bash
PMEM_IS_PMEM_FORCE=1 ./build/test_pmem 64 1000000 1
```


## Miscellaneous

We noticed a possible `mmap` bug on our testing environment: `MAP_SHARED_VALIDATE` is incompatible with `MAP_FIXED_NOREPLACE`.
To ensure safe memory mapping, we modified the original PMDK to use `MAP_SHARED` rather than `MAP_SHARED_VALIDATE`, which has the same functionality as the former one except for extra flag validation.
For a more detailed explanation and minimal reproducible code, please check out our [blog post](https://blog.haoxp.xyz/posts/mmap-bug/) about this issue.
