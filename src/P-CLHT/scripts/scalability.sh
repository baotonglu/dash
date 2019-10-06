#!/bin/bash

cores=$1;
shift;

source scripts/lock_exec;
source scripts/config;

prog=$1;
shift;
params="$@";


echo "#cores  throughput  %linear scalability"

printf "%-8d" 1;
thr1=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-12d" $thr1;
printf "%-8.2f" 100.00;
printf "%-8d\n" 1;

for c in $cores
do
    if [ $c -eq 1 ]
    then
	continue;
    fi;

    printf "%-8d" $c;
    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-12d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-8.2f" $linear_p;
    printf "%-8.2f\n" $scl;

done;

source scripts/unlock_exec;
