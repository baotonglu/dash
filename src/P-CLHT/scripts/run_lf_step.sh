#!/bin/bash

out_file="data/lf_step_"$(date | awk '// {print $2"_"$3}')".dat";
echo "Output file: $out_file";
printf "" > $out_file;

initials="8 16 32 64 128 256 512 1024 2048 4096 8192";
updates="0 1 10 20 50 100";

for i in $initials
do
    for u in $updates
    do
	r=$((2*$i));	
	settings="-i$i -r$r -u$u";
	echo "## $settings" | tee -a $out_file;
	./scripts/scalability3.sh "9 10 11 12" ./bin/lf-ll ./bin/lf-sl ./bin/lf-ht $settings -l4 | tee -a $out_file;
    done;
done;
