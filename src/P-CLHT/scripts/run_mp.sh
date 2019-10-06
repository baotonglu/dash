#!/bin/bash

for n in {3..48..3}
do
echo "./mp " $1 " " $2 " " $4 " " ${n}
./throughput_mp $1 ${n} $2 $3 $4 $5 $6 $7 $8 >> results/throughput_mp.txt
done

echo "./sequential " $1 " " $2 " " $4 " " 1
./sequential $1 1 $2 $3 $4 $5 $6 $7 $8 >> results/sequential.txt

cd plot
./tput_mp.sh
cp tput_mp.eps "../throughput_mp_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8".eps"
./speedup_mp.sh
cp speedup_mp.eps "../speedup_mp_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8".eps"
cd ..

rm results/*.txt plot/*.eps