#!/bin/bash

./run_rep.sh $@ | gawk '// { one=$1; sum+=$2; i++ } END { print one, sum/i }';


