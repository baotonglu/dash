#!/bin/bash

out_file="data/compare_ht_"$(date | gawk '// {print $2"_"$3}')".dat";
echo "Output file: $out_file";
printf "" > $out_file;

initials="64 512 2048 8192";
updates="0 10 100";
load_factors="4 8"


for i in $initials
do
    r=$((2*$i));	
    for u in $updates
    do
	for l in $load_factors
	do
	    settings="-i$i -r$r -u$u -l$l";
	    echo "## $settings" | tee -a $out_file;
	    ./scripts/scalability5.sh socket "./bin/lf-ht" "./bin/lb-ht_ticket" "./bin/lb-ht_gl_ticket" "./bin/lb-ht_ticket -x2" "./bin/lb-ht_gl_ticket -x2"  $settings | tee -a $out_file;
	done;
    done;
done;
