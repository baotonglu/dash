#!/bin/bash

out_file="data/compare_lf_hy_noise_"$(date | gawk '// {print $2"_"$3}')".dat";
echo "Output file: $out_file";
printf "" > $out_file;

initials="64 512 2048 8192";
updates="1 10 100";
load_factors="2 4"


for i in $initials
do
    r=$((2*$i));	
    for u in $updates
    do
	for l in $load_factors
	do
	    settings="-i$i -r$r -u$u -l$l";
	    echo "## $settings" | tee -a $out_file;
	    ./scripts/scalability2.sh socket "./lf-ht" "./throughput_ticket"  $settings | tee -a $out_file;
	done;
    done;
done;
