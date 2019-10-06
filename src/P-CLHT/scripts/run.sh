#!/bin/bash
if [ `uname` = "SunOS" ]
then
for n in {1..64}
do
echo "./hclh " $1 " " $2 " " $4 " " ${n}
./latency_hclh $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_hclh.txt
./throughput_hclh $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_hclh.txt
done

for n in {1..64}
do
echo "./ttas " $1 " " $2 " " $4 " " ${n}
./latency_ttas $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_ttas.txt
./throughput_ttas $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_ttas.txt
done

for n in {1..64}
do
echo "./mcs " $1 " " $2 " " $4 " " ${n}
./latency_mcs $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_mcs.txt
./throughput_mcs $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_mcs.txt
done

for n in {1..64}
do
echo "./array " $1 " " $2 " " $4 " " ${n}
./latency_array $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_array.txt
./throughput_array $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_array.txt
done

for n in {1..64}
do
echo "./ticket " $1 " " $2 " " $4 " " ${n}
./latency_ticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_ticket.txt
./throughput_ticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_ticket.txt
done

for n in {1..64}
do
echo "./mutex " $1 " " $2 " " $4 " " ${n}
./latency_mutex $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_mutex.txt
./throughput_mutex $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_mutex.txt
done

for n in {1..64}
do
echo "./hticket " $1 " " $2 " " $4 " " ${n}
./latency_hticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_hticket.txt
./throughput_hticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_hticket.txt
done

else
# Other UNIX (Linux, etc.) specific stuff
for n in {1..48}
do
echo "./hclh " $1 " " $2 " " $4 " " ${n}
./latency_hclh $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_hclh.txt
./throughput_hclh $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_hclh.txt
done

for n in {1..48}
do
echo "./ttas " $1 " " $2 " " $4 " " ${n}
./latency_ttas $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_ttas.txt
./throughput_ttas $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_ttas.txt
done

for n in {1..48}
do
echo "./mcs " $1 " " $2 " " $4 " " ${n}
./latency_mcs $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_mcs.txt
./throughput_mcs $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_mcs.txt
done

for n in {1..48}
do
echo "./array " $1 " " $2 " " $4 " " ${n}
./latency_array $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_array.txt
./throughput_array $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_array.txt
done

for n in {1..48}
do
echo "./ticket " $1 " " $2 " " $4 " " ${n}
./latency_ticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_ticket.txt
./throughput_ticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_ticket.txt
done

for n in {1..48}
do
echo "./mutex " $1 " " $2 " " $4 " " ${n}
./latency_mutex $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_mutex.txt
./throughput_mutex $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_mutex.txt
done

for n in {1..48}
do
echo "./hticket " $1 " " $2 " " $4 " " ${n}
./latency_hticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/latency_hticket.txt
./throughput_hticket $1 ${n} $2 $3 $4 $5 $6 $7 >> results/throughput_hticket.txt
done
fi

echo "./sequential " $1 " " $2 " " 1
./sequential $1 1 $2 $3 $4 $5 $6 $7 $8 >> results/sequential.txt

if [ `uname` = "SunOS" ]
then
cd plot_sparc
else
cd plot
fi
./put_acq.sh
./put_rel.sh
./put_opt.sh
./put_tot.sh
./get_acq.sh
./get_rel.sh
./get_opt.sh
./get_tot.sh
./remove_acq.sh
./remove_rel.sh
./remove_opt.sh
./remove_tot.sh
./tput.sh
cp tput.eps "../throughput_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7".eps"
./speedup.sh
cp speedup.eps "../speedup_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7".eps"
cd ..

cd tex
if [ `uname` = "SunOS" ]
then
/home/guerel/local/texlive/2012/bin/sparc-solaris/latex latency_sparc.tex
/home/guerel/local/texlive/2012/bin/sparc-solaris/dvips latency_sparc.dvi
cp latency_sparc.ps "../latency_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7".ps"
cd ..
rm results/*.txt plot_sparc/*.eps tex/*.ps
else
latex latency.tex
dvipdf latency.dvi
cp latency.pdf "../latency_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7".pdf"
cd ..
rm results/*.txt plot/*.eps tex/*.pdf
fi
