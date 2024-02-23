#!/bin/bash
#SBATCH -J ArtVandeley -o import_reads.out -t 00:10:00 -p short --mem=10G --nodes=1 --ntasks-per-node=1
module load QIIME2/2021.8

## importing all the demuxed reads
qiime tools import --type 'SampleData[SequencesWithQuality]' \
--input-path ../metadata/manifest.tsv \
--input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza

## visualize demuxed reads
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
