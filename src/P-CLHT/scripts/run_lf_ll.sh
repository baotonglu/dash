#!/bin/bash


initials="8 16 32 64 128 256 512 1024 2048 4096 8192";
updates="1 10 20 50 100";

for i in $initials
do
    for u in $updates
    do
	r=$((2*$i));	
	settings="-i$i -r$r -u$u";
	echo "## $settings";
	./scripts/scalability2.sh socket ./bin/lf-ll ./bin/lf-ll_padded $settings
    done;
done;
