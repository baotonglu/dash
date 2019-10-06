#!/bin/bash

## global  --------------------------------------------------------------------------------------------------------
out_folder="./data";
cores="all";
duration="2000";

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
update=10;
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## settings --------------------------------------------------------------------------------------------------------
initial=1000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;

## settings --------------------------------------------------------------------------------------------------------
initial=10000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;

## settings --------------------------------------------------------------------------------------------------------
initial=20000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;

## settings --------------------------------------------------------------------------------------------------------
initial=50000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
update=20;
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## settings --------------------------------------------------------------------------------------------------------
initial=1000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;

## settings --------------------------------------------------------------------------------------------------------
initial=10000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;

## settings --------------------------------------------------------------------------------------------------------
initial=20000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;

## settings --------------------------------------------------------------------------------------------------------
initial=50000;
range=$((2*$initial));

load_factor=1;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;
load_factor=2;
out="$out_folder/ht.i$initial.u$update.l$load_factor.dat";
./scripts/scalability.sh "$cores" ./bin/lf-ht -d$duration -i$initial -r$range -u$update -l$load_factor | tee $out;


