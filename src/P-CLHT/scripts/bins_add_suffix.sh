#!/bin/bash

if [ $# -gt 1 ];
then
    suffix=$1;
    base_name=$2;
    echo  "Using suffix: $suffix on pattern: \""$base_name"-..$\"";

    for f in $(ls bin/* | grep $base_name"-..$");
    do
	echo "Renaming: $f to $f"_$suffix;
	mv $f $f"_"$suffix;
    done;

else
    echo "Please provide the: suffix base_pattern";
    echo "Example: $0 aligned lf";
fi;
