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

update_all="0 20 50"
initial_all="1024 16384"
load_factor_all="2 8";

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

		    dat="$dat_folder/$app.i$initial.l$load_factor.u$update.dat";
		    gp="$gp_folder/$app.i$initial.l$load_factor.u$update.gp";
		    eps="$out_folder/$app.l$load_factor.i$initial.u$update.eps";
		    title="Lock-free $app / Load factor: $load_factor / Size: $initial / Update: $update";

		    cp $gp_template $gp
		    cat << EOF >> $gp
set title "$title"; 
set output "$eps";
plot \\
"$dat" using 1:(\$2) title  "lb - Througput" ls 2 with lines, \\
"$dat" using 1:(\$4) axis x1y2 title  "lb - Scalability" ls 1 with points, \\
"$dat" using 1:(\$5) title  "lz - Througput" ls 4 with lines, \\
"$dat" using 1:(\$7) axis x1y2 title  "lz - Scalability" ls 3 with points, \\
"$dat" using 1:(\$8) title  "lf - Througput" ls 6 with lines, \\
"$dat" using 1:(\$10) axis x1y2 title  "lf - Scalability" ls 5 with points 
EOF

		    gnuplot $gp;
		done;
	    else		# not hash table, aka not load_factor

		dat="$dat_folder/$app.i$initial.u$update.dat";
		gp="$gp_folder/$app.i$initial.u$update.gp";
		eps="$out_folder/$app.i$initial.u$update.eps";
		title="Lock-free $app / Size: $initial / Update: $update";

		cp $gp_template $gp


		if [ "$app" = "ll" ];
		then
		    cat << EOF >> $gp
set title "$title"; 
set output "$eps";
plot \\
"$dat" using 1:(\$2) title  "lb - Througput" ls 2 with lines, \\
"$dat" using 1:(\$4) axis x1y2 title  "lb - Scalability" ls 1 with points, \\
"$dat" using 1:(\$5) title  "lz - Througput" ls 4 with lines, \\
"$dat" using 1:(\$7) axis x1y2 title  "lz - Scalability" ls 3 with points, \\
"$dat" using 1:(\$8) title  "lf - Througput" ls 6 with lines, \\
"$dat" using 1:(\$10) axis x1y2 title  "lf - Scalability" ls 5 with points 
EOF

		else		# sl

		    cat << EOF >> $gp
set title "$title"; 
set output "$eps";
plot \\
"$dat" using 1:(\$2) title  "lb - Througput" ls 2 with lines, \\
"$dat" using 1:(\$4) axis x1y2 title  "lb - Scalability" ls 1 with points, \\
"$dat" using 1:(\$5) title  "lf - Througput" ls 6 with lines, \\
"$dat" using 1:(\$7) axis x1y2 title  "lf - Scalability" ls 5 with points 
EOF
		fi;

		gnuplot $gp;
	    fi;
	done;
    done;
done;


