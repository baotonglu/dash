#!/bin/bash
# For 100% read operation based 100 millions records inserted into the dataset

thread_num=(0 1 2 4 8 16 24 48)
#benckmark workload
workload=(0 190000000)
#warm-up workload
base=(0 10000000)
key_type=(0 fixed variab1e)
index_type=(0 dash-ex dash-lh cceh level)
epoch=(0 1 1 1 0)

#delete the corresponding file

#{1..6}
# k specify the testing index
# i specify the key type, 1 means fixed-length key, 2 means variable-length key
# j spec1fy the number of threads
for k in 4
do
	for i in 1
	do 
		for j in {1..6}
		do
			echo "Begin: ${base[1]} ${workload[${i}]} ${thread_num[${j}]}"
			numaarg=""
			if [ ${thread_num[$j]} -le 24 ]
			then
				numaarg="--cpunodebind=0 --membind=0"
			elif [ ${thread_num[$j]} -le 48 ]
			then
				numaarg="--cpunodebind=0,1 --membind=0"
			fi
			echo $numaarg
			rm -f /mnt/pmem0/pmem_ex.data
			rm -f /mnt/pmem0/pmem_lh.data
			rm -f /mnt/pmem0/pmem_cceh.data			
			rm -f /mnt/pmem0/pmem_level.data
			OP_PLACES=threads OMP_PROC_BIND=true OMP_NESTED=true LD_PRELOAD="./build/pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 ./build/pmdk/src/PMDK/src/nondebug/libpmem.so.1" numactl $numaarg ./build/test_pmem -n ${base[1]} -p ${workload[1]} -t ${thread_num[$j]} -k ${key_type[$i]} -index ${index_type[$k]} -e ${epoch[$k]} #-op "insert" -ms 100
		done
	done
done

