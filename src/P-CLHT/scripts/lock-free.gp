set xlabel "Number of threads"
set ylabel "Throughput (Ops/s)"
set y2label "Scalability"

#nomirror - no tics on the right and top
#scale - the size of the tics
set xtics 6 nomirror scale 2
set ytics auto nomirror scale 2
set y2tics auto nomirror scale 2

#remove top and right borders
set border 3 back
#add grid
set style line 12 lc rgb '#808080' lt 2 lw 1
set grid back ls 12

set xrange [0:]
set yrange [0:]

#the size of the graph
set size 2.0, 1.0

#some nice line colors
#lc: line color; lt: line type (1 - conitnuous); pt: point type
#ps: point size; lw: line width
set style line 1 lc rgb '#0060ad' lt 1 pt 2 ps 1.5 lw 5 
set style line 2 lc rgb '#dd181f' lt 1 pt 5 ps 2 lw 5
set style line 3 lc rgb '#8b1a0e' pt 1 ps 1.5 lt 1 lw 5
set style line 4 lc rgb '#5e9c36' pt 6 ps 2 lt 2 lw 5
set style line 5 lc rgb '#663399' lt 1 pt 3 ps 1.5 lw 5 
set style line 8 lc rgb '#299fff' lt 1 pt 8 ps 2 lw 5
set style line 9 lc rgb '#ff299f' lt 1 pt 9 ps 2 lw 5
set style line 6 lc rgb '#cc6600' lt 4 pt 4 ps 2 lw 5
set style line 7 lc rgb '#cccc00' lt 1 pt 7 ps 2 lw 5

#move the legend to a custom position (can also be moved to absolute coordinates)
set key outside

set terminal postscript color  "Helvetica" 24 eps enhanced

#set term tikz standalone color solid size 5in,3in
#set output "test.tex"

#for more details on latex, see http://www.gnuplotting.org/introduction/output-terminals/
#set term epslatex #size 3.5, 2.62 #color colortext
#size can also be in cm: set terminal epslatex size 8.89cm,6.65cm color colortext
#set output "test.eps"

# plot \
# "stress_rw_TTAS.txt" using 1:2 title "TTAS" ls 1 with linespoints, \
# "stress_rw_ARRAY.txt" using 1:2 title "Array" ls 2 with linespoints, \
# "stress_rw_MCS.txt" using 1:2 title "MCS" ls 3 with linespoints, \
# "stress_rw_HTICKET.txt" using 1:2 title "HTicket" ls 8 with linespoints, \
# "stress_rw_SPINLOCK.txt" using 1:2 title "Spinlock" ls 9 with linespoints, \
# "stress_rw_HCLH.txt" using 1:2 title "HCLH" ls 4 with linespoints, \
# "stress_rw_CLH.txt" using 1:2 title "CLH" ls 5 with linespoints, \
# "stress_rw_MUTEX.txt" using 1:2 title "Mutex" ls 7 with linespoints, \
# "stress_rw_TICKET.txt" using 1:2 title "Ticket" ls 6 with linespoints 

