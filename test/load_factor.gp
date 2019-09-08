#set terminal postscript eps color "Times-Roman" 24
set term pdfcairo lw 2 font "times new roman, 14"
#set size ratio 0.67 
set output "load_factor.pdf"

#set grid y
set style fill pattern border -1
#set title "Load Factor"
set xlabel "Sub-hash-table Size (KB)"
set ylabel "Maximum Load Factor"
#set xtics font ",30,Bold"
#set ytics font ",30,Bold"
set xrange [1:32]
set yrange [0.4:1]
set logscale x
set xtics (1,2,4,8,16,32)

set key bottom left
#set key Left
set key reverse

set key font ",10"

plot "./data/bucketized.csv" using 1:2 with linespoints pt 7 ps 1 title "Bucketized", "./data/probe.csv" using 1:2 with linespoints pt 6 ps 1 title "+Probing",
"./data/balanced.csv" using 1:2 with linespoints pt 7 ps 1 title "+ Balanced Insert", "./data/displacement.csv" using 1:2 with linespoints pt 7 ps 1 title "+ Displacement",