#!/bin/sh

#  qiime_labsScript.sh
#
#
#  Created by Evan Marshman on 7/25/18.
#

name=$1
lab=$2
depth=$3
mapper=$4
echo $2
echo $3
echo $4

mkdir $2


qiime feature-table filter-samples \
  --i-table table-$name"".qza \
  --m-metadata-file $mapper \
  --p-where $5 \
  --o-filtered-table table-$lab""-$name"".qza
  
qiime feature-table summarize \
  --i-table table-$lab""-$name"".qza \
  --o-visualization table-$lab""-$name"".qzv \
  --m-sample-metadata-file $mapper
  
qiime taxa barplot \
  --i-table table-$lab""-$name"".qza \
  --i-taxonomy taxonomy-$name"".qza \
  --m-metadata-file $mapper \
  --o-visualization taxa-bar-plots-$lab""-$name"".qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-$name"".qza \
  --i-table table-$lab""-$name"".qza \
  --p-sampling-depth $depth \
  --m-metadata-file $mapper \
  --output-dir core-metrics-results_$lab""-$name""
  
qiime taxa collapse \
  --i-table table-$lab""-$name"".qza \
  --i-taxonomy taxonomy-$name"".qza \
  --p-level 6 \
  --o-collapsed-table table-$lab""-$name""_l6.qza
  
qiime composition add-pseudocount \
  --i-table table-$lab""-$name""_l6.qza \
  --o-composition-table comp-table-$lab""-$name""_l6.qza
  
qiime composition ancom \
  --i-table comp-table-$lab""-$name""_l6.qza \
  --m-metadata-file $mapper \
  --m-metadata-column genotype_class \
  --o-visualization l6-ancom-genotype-$lab""-$name"".qzv
  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_$lab""-$name""/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $mapper \
  --m-metadata-column genotype_class \
  --o-visualization core-metrics-results_$lab""-$name""/unweighted-unifrac-genotype-significance.qzv \
  --p-pairwise
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_$lab""-$name""/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $mapper \
  --p-formula "genotype_class*sex" \
  --o-visualization permanova-sex-type-$lab""-$name""

mv *$2* $2/

