#!/bin/bash

reps=$1;
shift;

for r in $(seq 1 1 $reps);
do
    printf "%-5d : " $r
    res=$(./hyht -b1 -i81920 -n80 -d1000 | grep "garbage:");
    echo $res;
done;


