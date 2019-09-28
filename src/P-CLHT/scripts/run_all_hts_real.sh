#!/bin/bash

initials="16 128 1024 8192 65536";
updates="0 1 10 20 50 100";
duration=1000;


for i in $initials;
do
    r=$((2*$i));
    for u in $updates;
    do
	params="-i$i -u$u      -b$i -r$r -d$duration";
	echo "## PARAMS: $params";
	./scripts/scalability4.sh socket hyht hyht_lock_ins lfht_dup lfht  $params;
    done;
done;
