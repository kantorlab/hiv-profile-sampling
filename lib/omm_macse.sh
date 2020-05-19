#!/bin/bash
SIF=$PWD/scratch/omm_macse_sha256.096cd4607b78cd6aaf0d8af1e232e43824405320d155c5343f9c0a713595976c.sif
IN=$(basename $1)
OUT=$PWD/$2
cp $1 /tmp/$IN
cd /tmp
DIR=omm_macse_$(echo $1 | tr '/' '_')
rm -rf $DIR
singularity run $SIF --out_dir $DIR --out_file_prefix omm_macse --in_seq_file $IN --java_mem 4000m
cp $DIR/omm_macse_final_align_NT.aln $OUT
rm -f $IN
rm -rf $DIR
