#!bin/sh
one_file=$1
outDir=$2
fastq-dump --gzip --split-3 -O ${outDir} ${one_file}