#!/bin/bash
# For 100% read operation based 100 millions records inserted into the dataset

thread_num=(0 1 2 4 8 16 24 48)
workload=(0 190000000)
base=(0 10000000)
key_type=(0 fixed variable)
index_type=(0 dash-ex dash-lh cceh level)
epoch=(0 1 1 1 0)

#delete the corresponding file

#{1..6}
for k in 2
do
	for i in 1
	do 
		for j in 1
		do
			echo "Begin: ${base[1]} ${workload[${i}]} ${thread_num[${j}]}"
			numaarg=""
			if [ ${thread_num[$j]} -le 24 ]
			then
				#numaarg="--cpunodebind=0 --membind=0 --physcpubind=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
				numaarg="--cpunodebind=0 --membind=0"
			elif [ ${thread_num[$j]} -le 48 ]
			then
				#numaarg="--cpunodebind=0,1 --membind=0,1 --physcpubind=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,2728,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47"
				numaarg="--cpunodebind=0,1 --membind=0"
			fi
			echo $numaarg
			rm -f /mnt/pmem0/pmem_ex.data
			rm -f /mnt/pmem0/pmem_lh.data
			rm -f /mnt/pmem0/pmem_cceh.data			
			rm -f /mnt/pmem0/pmem_level.data
			OP_PLACES=threads OMP_PROC_BIND=true OMP_NESTED=true PMEM_IS_PMEM_FORCE=1 LD_PRELOAD="./build/pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 ./build/pmdk/src/PMDK/src/nondebug/libpmem.so.1" numactl $numaarg ./build/test_pmem -n ${base[1]} -p ${workload[1]} -t ${thread_num[$j]} -k ${key_type[$i]} -index ${index_type[$k]} -e ${epoch[$k]} #-op "recovery" -ms 100 #-distribution "zipfian"
		done
	done
done

#./extract_plot2.sh
# OMP_PROC_BIND=true

