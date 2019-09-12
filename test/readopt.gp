#set terminal postscript eps color "Times-Roman" 24
set term pdfcairo lw 2 font "times new roman, 20"
#set size ratio 0.67 
set output "read_opt.pdf"

set grid y
set title "2 stash buckets, 1 thread"
set style data histograms
set style fill solid 0.4 border
set style histogram clustered gap 2
#set xlabel "Operation"
set ylabel "Throughput(ops/nsec)"
set xtics font ",15,Bold"
#set ytics font ",30,Bold"
#set yrange [0:1.5]

set key top left
set key Left
set key reverse

plot './data/read_opt_2_1.csv' using 2:xticlabels(1) fs pattern 3 lt 1 lw 1 lc rgb 'black' title columnheader(2),\
 '' using 3:xticlabels(1) fs pattern 4 lt 1 lw 1 lc rgb 'black' title columnheader(3) 
#plot './data/read_opt.csv' using 2:xtic(1) linecolor 'black' title 'Static Insert Throughput(200 millions)'