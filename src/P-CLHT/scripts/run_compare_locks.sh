#!/bin/bash

p="lb-ll -u10";
echo "## settings: $p";
./scripts/run_compare.sh $p
p="lb-ll -u10 -i128 -r256";
echo "## settings: $p";
./scripts/run_compare.sh $p

p="lb-ll -u10 -x2";
echo "## settings: $p";
./scripts/run_compare.sh $p
p="lb-ll -u10 -i128 -r256 -x2";
echo "## settings: $p";
./scripts/run_compare.sh $p

p="lb-sl  -u10";
echo "## settings: $p";
./scripts/run_compare.sh $p
p="lb-sl -u10 -i128 -r256";
echo "## settings: $p";
./scripts/run_compare.sh $p

p="lb-ht -u10 -l4";
echo "## settings: $p";
./scripts/run_compare.sh $p
p="lb-ht -u10 -l4 -i128 -r256";
echo "## settings: $p";
./scripts/run_compare.sh $p


p="lb-ht -u10 -l4 -x2"
echo "## settings: $p";
./scripts/run_compare.sh $p
p="lb-ht -u10 -l4 -i128 -r256 -x2";
echo "## settings: $p";
./scripts/run_compare.sh $p

