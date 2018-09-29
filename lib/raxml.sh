#!/bin/bash

set -e

FA=$1
N=$2
MODEL=$3
CPUS=$4
OUTDIR=$5
NAME=$6

# create deterministic seed from alignment
SEED=$(cat $FA | md5sum | grep -Po "\d" | tail -10 | tr -d "\n")

raxmlHPC-PTHREADS-AVX -T $CPUS -s $FA -p $SEED -x $SEED -f a -N $N -m $MODEL -w $PWD/$OUTDIR -n $NAME
