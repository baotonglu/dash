#!/bin/bash

pattern=$1;
shift;

num_bins=$(ls -1 bin/*$pattern* | wc -l);
echo "# of executables found: $num_bins";

if [ $num_bins -gt 10 ];
then
    echo "** Cannot execute more than 9 benchmarks at a time";
    exit -1;
fi;

if [ $num_bins -eq 0 ];
then
    echo "** No executables found";
    exit -1;
fi;

bins=$(ls bin/*$pattern*);

./scripts/scalability$num_bins.sh socket $bins $@;
