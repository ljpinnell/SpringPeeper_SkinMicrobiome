#!/bin/bash
#SBATCH -t 00:05:00 -p short --mem=1G -n 1 --ntasks=1
module load QIIME2/2023.2

## make phyloseq dir
mkdir exported

## export to phyloseq dir
qiime tools export --input-path rep-seqs.qza --output-path exported/
qiime tools export --input-path taxonomy.qza --output-path exported/
qiime tools export --input-path rooted-tree.qza --output-path exported/
qiime tools export --input-path table.qza --output-path exported/

## still need to change headers and ammend biom table with taxonomy
