#!/bin/bash
set -e
STARTDIR=$PWD
SIF=$PWD/scratch/omm_macse_sha256.096cd4607b78cd6aaf0d8af1e232e43824405320d155c5343f9c0a713595976c.sif
IN=$PWD/$1
OUT=$PWD/$2
SIZE=$(grep -v '^>' $IN | tr -d '\n' | wc -c)
if [ $SIZE -gt 0 ]
then
  DIR="omm_macse_${1//\//_}"
  mkdir -p $DIR
  cd $DIR
  cp $IN input.fa
  singularity run $SIF --out_dir omm_macse --out_file_prefix omm_macse --in_seq_file input.fa --java_mem 4000m
  cp omm_macse/omm_macse_final_align_NT.aln $OUT
  cd $STARTDIR
  rm -rf $DIR
else
  cp $IN $OUT
fi
