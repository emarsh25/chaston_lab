#!/bin/sh

#  standard_qiime2.sh
#
#
#  Created by John Chaston on 7/25/18.
#  Edited by Evan Marshman on 5/13/21.
#

tablename=$1
treename=$2
name=$3
sdepth=$4
mapper=$5
echo $2
echo $3
echo $4
echo $5

mkdir $3

## just E, B, and ED overtime
qiime feature-table filter-samples \
--i-table table-$tablename"".qza \
--m-metadata-file $mapper \
--p-where "$6" \
--o-filtered-table table_$name"".qza

# make taxonomic assignments
#qiime feature-classifier classify-sklearn \
#--i-classifier gg-13-8-99-515-806-nb-classifier.qza \
#--i-reads filtered_rep-seqs.qza \
#--o-classification taxonomy.qza

# filter out wolbachia reads
#qiime taxa filter-table --i-table table_$name"".qza --i-taxonomy filtered_taxonomy.qza --p-exclude Wolbachia --o-filtered-table table_$name""_noW.qza

qiime feature-table summarize \
--i-table table_$name"".qza \
--o-visualization table_$name"".qzv \
--m-sample-metadata-file $mapper

## run the core metrics
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree-$treename"".qza \
--i-table table_$name"".qza \
--p-sampling-depth $sdepth \
--m-metadata-file $mapper \
--output-dir core-metrics-results-$name""

## run the alpha rarefaction
qiime diversity alpha-rarefaction \
--i-table table_$name"".qza \
--i-phylogeny rooted-tree-$treename"".qza \
--p-max-depth $sdepth \
--m-metadata-file $mapper \
--o-visualization alpha-rarefaction-$name"".qzv

## creates ace richness matrix
qiime diversity alpha \
  --i-table table_$name"".qza \
  --p-metric ace \
  --o-alpha-diversity core-metrics-results-$name""/ace_vector.qza

# convert the distance matrices and pcoas
cd core-metrics-results-$name""
qiime tools export \
--input-path weighted_unifrac_distance_matrix.qza \
--output-path weighted_unifrac_distance_matrix/
qiime tools export \
--input-path unweighted_unifrac_distance_matrix.qza \
--output-path unweighted_unifrac_distance_matrix/
qiime tools export \
--input-path unweighted_unifrac_pcoa_results.qza \
--output-path unweighted_unifrac_pcoa_results/
qiime tools export \
--input-path weighted_unifrac_pcoa_results.qza \
--output-path weighted_unifrac_pcoa_results/
qiime tools export \
--input-path bray_curtis_pcoa_results.qza \
--output-path bray_curtis_pcoa_results/
qiime tools export \
--input-path bray_curtis_distance_matrix.qza \
--output-path bray_curtis_distance_matrix/
qiime tools export \
--input-path ace_vector.qza \
--output-path ace_richness/
qiime tools export \
--input-path shannon_vector.qza \
--output-path shannon_richness/
qiime tools export \
--input-path rarefied_table.qza \
--output-path rarefied_table/
biom convert -i rarefied_table/feature-table.biom -o rarefied_table.txt --to-tsv
cd ..

mv *$3* $3/
