module load QIIME2/2023.2

# adding taxonomy to BIOM file
biom add-metadata -i feature-table.biom -o table-with-taxonomy.biom --observation-metadata-fp biom-taxonomy.tsv --sc-separated taxonomy
