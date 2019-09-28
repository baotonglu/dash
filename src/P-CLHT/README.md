## P-CLHT: Persistent Cache-Line Hash Table

`P-CLHT` is a crash consistent version of [Cache-Line Hash Table](https://dl.acm.org/citation.cfm?id=2694359) (CLHT).
CLHT is a cache-friendly hash table which restricts each bucket to be of the size of a cache line. 
CLHT is an unordered index only supporting point queries.

**Conversion**. `CLHT-LB` using lock-based writes for concurrency is converted into `P-CLHT` by adding cache 
line flushes and memory fences after each critical volatile store. 

**Performance**. Compared with [CCEH](https://www.usenix.org/conference/fast19/presentation/nam) that is a 
state-of-the-art unordered index, `P-CLHT` shows **2.38x**, **1.35x**, and **1.25x** better performance in 
YCSB workload A, B, C respectively using random integer keys while **0.37x** worse in Load workload.

**Support**. `P-CLHT` supports Insert, Delete, and Point Lookup operations. Each operation works for only integer keys.

**Use Case**. `P-CLHT` provides the superior performance of insertion and point lookup, even if not supporting
range scans. Therefore, it would be appropriate to be used for the applications only consisting of point queries.

## Build & Run

#### Build

```
$ bash compile.sh lb
$ mkdir build
$ cd build
$ cmake ..
$ make -j
```

#### Run

```
$ ./example 10000 4

usage: ./example [n] [nthreads]
n: number of keys (integer)
nthreads: number of threads (integer)
```
