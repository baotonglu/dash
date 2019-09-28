#!/bin/bash

if [ $(uname -n) = "maglite" ];
then
    MAKE="/home/trigonak/sw/make/make";
else
    MAKE=make
fi;

lock=mutex
$MAKE $lock GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_$lock lb
$MAKE $lock
./scripts/bins_add_suffix.sh $lock lb

lock=spin
$MAKE $lock GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_$lock lb
$MAKE $lock
./scripts/bins_add_suffix.sh $lock lb

lock=tas
$MAKE $lock GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_$lock lb
$MAKE $lock
./scripts/bins_add_suffix.sh $lock lb

lock=ticket
$MAKE $lock GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_$lock lb
$MAKE $lock
./scripts/bins_add_suffix.sh $lock lb


if [ $(uname -n) = "lpd48core" ] || [ $(uname -n) = "diassrv8" ];
then
    lock=hticket
    $MAKE $lock GRANULARITY=GLOBAL_LOCK;
    ./scripts/bins_add_suffix.sh gl_$lock lb

    lock=clh
    $MAKE $lock GRANULARITY=GLOBAL_LOCK;
    ./scripts/bins_add_suffix.sh gl_$lock lb
fi;
