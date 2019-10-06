#!/bin/bash

reps=$1;
shift;

if [ $(uname -n) = "smal1.sics.se" ];
then
    run=./run
elif [ $(uname -n) = "parsasrv1.epfl.ch" ];
then
    run=./run
fi


execute=$@;

for rep in $(seq 1 1 $reps)
do
     $run ./$execute
done;


