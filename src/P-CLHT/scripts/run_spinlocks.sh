#!/bin/bash

## global  --------------------------------------------------------------------------------------------------------
dat_folder="./data";
cores="all";
duration="1000";

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
update=10;
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## settings --------------------------------------------------------------------------------------------------------
app_all="lb-ll lb-sl lb-ht";
update_all="10 20";
initial_all="1000 10000 20000 50000";
load_factor_all="1 2 10";

for update in $update_all
do
    echo "* Update: $update";
    for initial in $initial_all
    do
	echo "** Size: $initial";
	range=$((2*$initial));

	for app in $app_all
	do
	    printf "*** Application: $app\n";
	    if [ "$app" = "ht" ];
	    then

		for load_factor in $load_factor_all
		do
		    echo "**** Load factor: $load_factor";
		    dat="$dat_folder/$app.i$initial.u$update.l$load_factor.dat";
		    ./scripts/scalability.sh "$cores" ./bin/$app -d$duration -i$initial -r$range -u$update | tee $dat;
		done;

	    else

		dat="$dat_folder/$app.i$initial.u$update.dat";
		./scripts/scalability.sh "$cores" ./bin/$app -d$duration -i$initial -r$range -u$update | tee $dat;
	    fi;
	done;
    done;
done;