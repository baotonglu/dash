#!/bin/bash

## global  --------------------------------------------------------------------------------------------------------
out_folder="./data";
cores="all";
duration="1000";


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## executables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lls='./bin/lb-ll "./bin/lb-ll -x2" ./bin/lf-ll' # -x2 does not work
sls='./bin/lb-sl ./bin/lf-sl'
hts='./bin/lb-ht "./bin/lb-ht -x2" ./bin/lf-ht' # -x2 does not work
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## settings --------------------------------------------------------------------------------------------------------
initials="1024 16384";
updates="0 20 50";
load_factors="2 8"

for initial in $initials; 
do
    echo "* -i$initial";

    range=$((2*$initial));

    for update in $updates;
    do
	echo "** -u$update";

	out="$out_folder/ll.i$initial.u$update.dat";
	./scripts/scalability3.sh "$cores" ./bin/lb-ll "./bin/lb-ll -x2" ./bin/lf-ll -d$duration -i$initial -r$range -u$update | tee $out;
	# out="$out_folder/sl.i$initial.u$update.dat";
	# ./scripts/scalability2.sh "$cores" $sls -d$duration -i$initial -r$range -u$update | tee $out;

	for load_factor in $load_factors; 
	do
	    echo "*** -l$load_factor";
	    out="$out_folder/ht.i$initial.l$load_factor.u$update.dat";
	    ./scripts/scalability3.sh "$cores" ./bin/lb-ht "./bin/lb-ht -x2" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
	done
    done
done

