set term pdfcairo lw 2 font "times new roman, 18"
#set size ratio 0.67 
set output "lock.pdf"

set grid y
set style fill pattern border -1
set xlabel "Number of threads" font ",28"
set ylabel "Throughput(Ops/usec)" font ",28"

set xtics font ",23,Bold"
set ytics font ",23,Bold"
set xrange [1:25]
#set yrange [0:10]
set xtics (1,2,4,8,16,24)
#set ytics (0,2,4,6,8,10)

set key top left
set key Left
set key reverse

set key font ",18"

plot "./data/pos_opt_lock.csv" using 1:2 with linespoints lw 2 pt 9 ps 1.3 title "Optimistic (pos search)",\
 "./data/neg_opt_lock.csv" using 1:2 with linespoints lw 2 pt 11 ps 1.3 title "Optimistic (neg search)",\
 "./data/pos_spin_lock.csv" using 1:2 with linespoints lw 2 pt 5 ps 1.3 title "Spin(pos search)",\
 "./data/neg_spin_lock.csv" using 1:2 with linespoints lw 2 pt 7 ps 1.3 title "Spin (neg search)"