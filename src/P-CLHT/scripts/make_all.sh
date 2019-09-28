#!/bin/bash

if [ $(uname -n) = "maglite" ];
then
    MAKE="/home/trigonak/sw/make/make";
else
    MAKE=make
fi;

$MAKE lockfree;
$MAKE ticket GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_ticket lb
$MAKE ticket
./scripts/bins_add_suffix.sh ticket lb
$MAKE hticket GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_hticket lb
$MAKE clh GRANULARITY=GLOBAL_LOCK;
./scripts/bins_add_suffix.sh gl_clh lb

