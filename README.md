# Dash: Dynamic and Scalable Hashing on Persistent Memory 

[![Build Status](https://dev.azure.com/haoxiangpeng/VeryPM/_apis/build/status/XiangpengHao.VeryPM?branchName=master)](https://dev.azure.com/haoxiangpeng/VeryPM/_build/latest?definitionId=2&branchName=master)

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
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_PMEM=ON .. 
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
