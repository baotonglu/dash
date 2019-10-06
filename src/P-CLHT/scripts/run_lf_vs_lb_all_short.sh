#!/bin/bash

out_file="data/lf_vs_lb_all_short_"$(date | gawk '// {print $2"_"$3}')".dat";
echo "Output file: $out_file";
printf "" > $out_file;

initials="8 32 128 512 2048 8192";
updates="0 1 10 100";

lb_ll_gl="./bin/lb-ll_gl_ticket -x2";
if [ $(uname -n) = "diassrv8" ];
then
    lb_ll_gl="./bin/lb-ll_gl_hticket -x2";
elif [ $(uname -n) = "lpd48core" ];
then
    lb_ll_gl="./bin/lb-ll_gl_clh -x2";
fi;

for i in $initials
do
    for u in $updates
    do
	r=$((2*$i));	
	settings="-i$i -r$r -u$u";
	echo "## $settings" | tee -a $out_file;
	./scripts/scalability8.sh all ./bin/lf-ll "./bin/lb-ll_ticket -x2" "$lb_ll_gl" ./bin/lf-sl ./bin/lb-sl_ticket ./bin/lf-ht "./bin/lb-ht_ticket -x2" "./bin/lb-ht_gl_ticket -x2"  $settings -l4 | tee -a $out_file;
    done;
done;
