# PM_Hashing
persistent memory friendly hashing index


## Build using CMake

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_PMEM=ON .. 
make -j
```

### Run on a real device

```bash
./build/test_pmem 64 1000000 1
```

### Simulate on DRAM 

```bash
PMEM_IS_PMEM_FORCE=1 ./build/test_pmem 64 1000000 1
```
