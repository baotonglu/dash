#!/bin/bash
# For 100% read operation based 100 millions records inserted into the dataset

rm -rf _CCEH

thread_num=(0 1 2 4 10 20 30 40 80)
workload=(0 20000000)
base=(0 1000000)

#{1..6}
for i in 1
do 
	for j in 1
	do
		#OMP_PLACES=threads OMP_PROC_BIND=true OMP_NESTED=true PMEM_IS_PMEM_FORCE=1 LD_PRELOAD=./build/pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 ./build/test_pmem 64 ${workload[$i]} ${thread_num[$j]}
		OMP_PLACES=threads OMP_PROC_BIND=true OMP_NESTED=true PMEM_IS_PMEM_FORCE=1 LD_PRELOAD=./build/pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 ./src/CCEH/test_cceh 64 ${workload[$i]} ${thread_num[$j]}
		printf "Done for cceh dm uni: %d %d\n" ${workload[$i]} ${thread_num[$j]}
	done
done
