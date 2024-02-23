#!/bin/bash
#SBATCH -J classify -o classify.out -t 08:00:00 -p medium --mem=100G --nodes=1 --ntasks=24
module load QIIME2/2023.2

## classifying
qiime feature-classifier classify-sklearn --i-classifier /scratch/group/vero_research/QIIME_Classifiers/qiime2-2023.2/SILVA138.1/515f-806r/silva-138.1-ssu-nr99-515f-806r-classifier.qza \
--i-reads rep-seqs_2.qza --o-classification taxonomy.qza

## making the no-chloro-no-mito directory
mkdir no-chloro-no-mito

## copy the taxonomy file into it
cp taxonomy.qza no-chloro-no-mito/

## filtering ASVs
qiime taxa filter-table --i-table table_2.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-table no-chloro-no-mito/table.qza

## filtering rep-seqs
qiime taxa filter-seqs --i-sequences rep-seqs_2.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-sequences no-chloro-no-mito/rep-seqs.qza

## chahnge into no
cd no-chloro-no-mito

## tree
## making a midpoint-rooted tree from our rep-seqs
qiime alignment mafft --i-sequences rep-seqs.qza  --p-n-threads 'auto' \
--o-alignment aligned-rep-seqs.qza
qiime alignment mask --i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --p-n-threads 'auto' \
--o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza
