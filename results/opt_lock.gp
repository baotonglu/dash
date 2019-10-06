#set terminal postscript eps color "Times-Roman" 24
set term pdfcairo lw 2 font "times new roman, 14"
#set size ratio 0.7 
set output "opt_lock.pdf"

set grid y
set style fill pattern border -1
#set title "Positive Search" font ",20"
set title "Negative Search" font ",20"
set xlabel "Number of threads" font ",20"
set ylabel "Throughput(ops/nsec)" font ",20"

#set xtics font ",30,Bold"
#set ytics font ",30,Bold"
set xrange [1:24]
#set yrange [0.4:1]
#set xtics (1,2,4,8,16,32,64,128)

set key top left
#set key Left
set key reverse

set key font ",14"

plot "./data/neg_spin_lock.csv" using 1:2 with linespoints pt 7 ps 1 title "Spin Lock",\
 "./data/neg_opt_lock.csv" using 1:2 with linespoints pt 6 ps 1 title "Optimistic Lock"
 