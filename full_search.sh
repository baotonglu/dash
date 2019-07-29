#!/bin/bash
# For 100% read operation based 100 millions records inserted into the dataset

thread_num=(0 1 2 5 10 20 30 40 80)
workload=(0 100000000)
base=(0 1000000)

#delete the corresponding file

#rm -f dm_cceh_search.txt
#rm -f dm_cceh_search_base.txt
#rm -f dm_level_search.txt
#rm -f pm_level_search.txt
#rm -f pm_cceh_search_base.txt
#rm pm_spinlock2.txt

#{1..6}
for i in 1
do 
	for j in 1
	do
		echo "Begin: ${base[1]} ${workload[${i}]} ${thread_num[${j}]}"
		numaarg=""
		if [ ${thread_num[$j]} -le 10 ]
		then
			numaarg="--cpunodebind=0 --membind=0 --physcpubind=0,4,8,12,16,20,24,28,32,36"
		elif [ ${thread_num[$j]} -le 20 ]
		then
			numaarg="--cpunodebind=0,1 --membind=0,1 --physcpubind=0,4,8,12,16,20,24,28,32,36,1,5,9,13,17,21,25,29,33,37"
		elif [ ${thread_num[$j]} -le 30 ]
                then
                        numaarg="--cpunodebind=0,1,2 --membind=0,1,2 --physcpubind=0,4,8,12,16,20,24,28,32,36,1,5,9,13,17,21,25,29,33,37,2,6,10,14,18,22,26,30,34,38"
		elif [ ${thread_num[$j]} -le 40 ]
		then
			numaarg="--cpunodebind=0,1,2,3 --membind=0,1,2,3 --physcpubind=0,4,8,12,16,20,24,28,32,36,1,5,9,13,17,21,25,29,33,37,2,6,10,14,18,22,26,30,34,38,3,7,11,15,19,23,27,31,35,39"
		elif [ ${thread_num[$j]} -le 80 ]
		then 
			numaarg="--cpunodebind=0,1,2,3 --membind=0,1,2,3"
		fi
		echo $numaarg
#		OMP_PLACES=threads OMP_PROC_BIND=true OMP_NESTED=true PMEM_IS_PMEM_FORCE=1 LD_PRELOAD=/usr/lib/libjemalloc.so numactl $numaarg ./new-src/cceh_bench 2 ${workload[$i]} ${thread_num[$j]} 
		OMP_PLACES=threads OMP_PROC_BIND=true OMP_NESTED=true PMEM_IS_PMEM_FORCE=1 LD_PRELOAD=/usr/lib/libjemalloc.so numactl $numaarg ./build/test_pmem 2 ${workload[$i]} ${thread_num[$j]} #>> cuckoo_finger.txt
		printf "Done for cceh dm uni: %d %d\n" ${workload[$i]} ${thread_num[$j]}
	done
done

#./extract_plot2.sh
# OMP_PROC_BIND=true

