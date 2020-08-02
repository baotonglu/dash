#!/bin/bash

# number of threads to run
thread_num=(0 1 2 4 8 16 24 48)
# benckmark workload, number of opeartions to run
workload=(0 190000000)
# warm-up workload, number of key-value to insert for warm-up
base=(0 10000000)
# type of keys
key_type=(0 fixed variab1e)
# which index to evaluate
index_type=(0 dash-ex dash-lh cceh level)
# whether to use to epcoh manager
epoch=(0 1 1 1 0)

# k specify the testing index, 1 = dash-ex, 2 = dash-lh, 3 = cceh, 4 = level
# i specify the key type, 1 means fixed-length key, 2 means variable-length key
# j spec1fy the number of threads, 1 means one threads, 2 means two threads, 3 means four threads...
for k in 1
do
	for i in 1
	do
		for j in 1
		do
			echo "Begin: ${base[1]} ${workload[1]} ${thread_num[${j}]}"
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
      LD_PRELOAD="./build/pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 \
      ./build/pmdk/src/PMDK/src/nondebug/libpmem.so.1" \
      numactl $numaarg ./build/test_pmem \
      -n ${base[1]} \
      -loadType 0 \
      -p ${workload[1]} \
      -t ${thread_num[$j]} \
      -k ${key_type[$i]} \
      -distribution "uniform" \
      -index ${index_type[$k]} \
      -e ${epoch[$k]} \
      -ed 1000 \
      -op "full" \
      -r 0.8 \
      -s 0.2 \
      -ms 100 \
      -ps 60
		done
	done
done
#done
