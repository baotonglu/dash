#!/bin/bash


app1=$1;
app2=$2;
shift; shift;
bu=$1;
shift; shift;
P=$@


printf "%-20s%-20s\n" $app1 $app2 $app3


# cores="2 6 36 48"
cores="32 36 48"

for n in $cores
do
    echo " -- $n threads";
    Pc="$bu $n $P"

    ./$app1 $Pc;
    ./$app2 $Pc;

    # P="-d1000 -n$n -l$acc -c1 -w1 -a0 -p0"; 
    # v1=$(./$app1 $P | awk '// { print $5 }');
    # v2=$(./$app2 $P | awk '// { print $5 }'); 
    # v3=0  # $(./$app3 $P | awk '// { print $5 }'); 

    # if [ $((v1)) -gt $((v2)) ];
    # then
    # 	if [ $((v1)) -gt $((v3)) ];
    # 	then
    # 	    max=$((v1));
    # 	else
    # 	    max=$((v3));
    # 	fi;
    # else
    # 	if [ $((v2)) -gt $((v3)) ];
    # 	then
    # 	    max=$((v2));
    # 	else
    # 	    max=$((v3));
    # 	fi;
    # fi;

    # printf "%-10d%-10.2f" $((v1)) $(echo "$v1/$max" | bc -l)
    # printf "%-10d%-10.2f" $((v2)) $(echo "$v2/$max" | bc -l)
    # printf "%-10d%10.2f\n" $((v3)) $(echo "$v3/$max" | bc -l)
    
done;

