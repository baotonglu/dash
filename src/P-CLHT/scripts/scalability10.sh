#!/bin/bash

cores=$1;
shift;

source scripts/lock_exec;
source scripts/config;

prog1="$1";
shift;
prog2="$1";
shift;
prog3="$1";
shift;
prog4="$1";
shift;
prog5="$1";
shift;
prog6="$1";
shift;
prog7="$1";
shift;
prog8="$1";
shift;
prog9="$1";
shift;
prog10="$1";
shift;
params="$@";


printf "#     %-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s\n" "$prog1" "$prog2" "$prog3" "$prog4" "$prog5" "$prog6" "$prog7" "$prog8" "$prog9" "$prog10";
echo "#cor  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb  throughput %linea scalb";

prog=$prog1;

printf "%-6d" 1;
thr1a=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1a;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog2;

thr1b=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1b;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog3;

thr1c=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1c;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog4;

thr1d=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1d;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog5;

thr1e=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1e;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog6;

thr1f=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1f;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog7;

thr1g=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1g;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog8;

thr1h=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1h;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog9;

thr1i=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1i;
printf "%-7.2f" 100.00;
printf "%-7d" 1;

prog=$prog10;

thr1j=$($run_script ./$prog $params -n1 | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
printf "%-11d" $thr1j;
printf "%-7.2f" 100.00;
printf "%-6d\n" 1;


for c in $cores
do
    if [ $c -eq 1 ]
    then
	continue;
    fi;

    printf "%-6d" $c;

    prog=$prog1;
    thr1=$thr1a;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog2;
    thr1=$thr1b;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog3;
    thr1=$thr1c;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog4;
    thr1=$thr1d;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog5;
    thr1=$thr1e;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog6;
    thr1=$thr1f;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog7;
    thr1=$thr1g;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog8;
    thr1=$thr1h;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog9;
    thr1=$thr1i;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-7.2f" $scl;

    prog=$prog10;
    thr1=$thr1j;

    thr=$($run_script ./$prog $params -n$c | grep "#txs" | cut -d'(' -f2 | cut -d. -f1);
    printf "%-11d" $thr;
    scl=$(echo "$thr/$thr1" | bc -l);
    linear_p=$(echo "100*(1-(($c-$scl)/$c))" | bc -l);
    printf "%-7.2f" $linear_p;
    printf "%-6.2f\n" $scl;


done;

source scripts/unlock_exec;
