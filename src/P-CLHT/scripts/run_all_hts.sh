#!/bin/bash

initials="1024 8192 65536";
updates="0 1 10 20";
duration=1000;


for i in $initials;
do
    r=$((2*$i));
    for u in $updates;
    do
	params="-i$i -u$u      -b$i -r$r -d$duration";
	echo "## PARAMS: $params";
	./scripts/scalability8.sh socket $(find . -name "clht*" -type f -executable -print) $params;
    done;
done;
