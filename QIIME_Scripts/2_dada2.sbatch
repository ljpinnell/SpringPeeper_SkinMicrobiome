#!/bin/bash
#SBATCH -J DADA2 -o dada2_2.out -t 04:00:00 -p medium --mem=200G --nodes=1 --ntasks-per-node=48
module load QIIME2/2023.2

qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza \
--p-n-threads 48 --o-table table_2.qza \
--p-trim-left 19 --p-trunc-len 194 \
--o-representative-sequences rep-seqs_2.qza --output-dir DADA2_2

# export dada2 stats
qiime tools export --input-path DADA2_2/denoising_stats.qza --output-path DADA2_2/stats
