#!/bin/bash

initials="16 64 128 512 1024 2048 4096 8192";
updates="1 10 20 50 100";
load_factors="2 4 8 16";
reps=20;

for i in $initials
do
    r=$((2*$i));	
    for u in $updates
    do
	for l in $load_factors
	do
	    settings="-i$i -r$r -u$u -l$l";
	    echo "## $settings";
	    for r in $(seq 1 1 $reps)
	    do
		echo "# rep: $r";
		./scripts/scalability.sh socket "./throughput_ticket"  $settings;
	    done;
	done;
    done;
done;


