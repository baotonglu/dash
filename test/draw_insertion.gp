set term pdfcairo lw 2 font "times new roman, 14"
#set size ratio 0.67 
set output "insertion.pdf"

#set grid y
set style fill pattern border -1
set xlabel "Number of threads" font ",23"
set ylabel "Throughput(operations/nsec)" font ",23"

set xtics font ",23,Bold"
set ytics font ",23,Bold"
set xrange [1:25]
set yrange [0:10]
set xtics (1,2,4,8,16,24)
set ytics (0,2,4,6,8,10)

set key top left
set key Left
set key reverse

set key font ",14"

plot "./insertion.csv" using 1:2 with linespoints lw 3 pt 7 ps 1.3 title "CCEH(4)",\
 "./insertion.csv" using 1:3 with linespoints lw 3 pt 9 ps 1.3 title "Level",\
  "./insertion.csv" using 1:4 with linespoints lw 3 pt 11 ps 1.3 title "DASH\\_EX",\
   "./insertion.csv" using 1:5 with linespoints lw 3 pt 13 ps 1.3 title "DASH\\_LH"