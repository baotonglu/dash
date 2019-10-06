#!/bin/sh

if [ $# -gt 0 ];
then
    app_all="$@";
else
    app_all="ll sl ht";
fi

dat_folder="./data";
gp_folder="./gp";
out_folder="./plots"
gp_template="./scripts/lock-free.gp";

update_all="10 20";
initial_all="1000 10000 20000 50000";
load_factor_all="1 2 10 100";

for update in $update_all
do
    echo "* Update: $update";
    for initial in $initial_all
    do
	echo "** Size: $initial";

	for app in $app_all
	do
	    printf "*** Application: $app\n";
	    if [ "$app" = "ht" ];
	    then

		for load_factor in $load_factor_all
		do
		    echo "**** Load factor: $load_factor";

		    dat="$dat_folder/$app.i$initial.u$update.l$load_factor.dat";
		    gp="$gp_folder/$app.i$initial.u$update.l$load_factor.gp";
		    eps="$out_folder/$app.l$load_factor.u$update.i$initial.eps";
		    title="Lock-free $app / Size: $initial / Update: $update / Load factor: $load_factor";

		    cp $gp_template $gp
		    cat << EOF >> $gp
set title "$title"; 
set output "$eps";
plot \\
"$dat" using 1:(\$2) title  "Througput" ls 2 with lines, \\
"$dat" using 1:(\$4) axis x1y2 title  "Scalability" ls 1 with points 
EOF
		done;

	    else

		dat="$dat_folder/$app.i$initial.u$update.dat";
		gp="$gp_folder/$app.i$initial.u$update.gp";
		eps="$out_folder/$app.u$update.$initial.eps";
		title="Lock-free $app / Size: $initial / Update: $update";

		cp $gp_template $gp
		cat << EOF >> $gp
set title "$title"; 
set output "$eps";
plot \\
"$dat" using 1:(\$2) title  "Througput" ls 2 with lines, \\
"$dat" using 1:(\$4) axis x1y2 title  "Scalability" ls 1 with points 
EOF

	    fi;
	done;
    done;
done;

gnuplot $gp_folder/*.gp;

