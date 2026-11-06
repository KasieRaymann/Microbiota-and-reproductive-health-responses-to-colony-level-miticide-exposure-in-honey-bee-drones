# QIIME2 Pipeline | Visualize alpha diversity graphs and generate stats

# Author:        Layla Jeries, with some ChatGPT assistance (chat log available)
# Date created:  2026-08-18
# Last modified: 2026-11-06

########################################## AI USAGE DISCLAIMER ###########################################
## This disclaimer has been included because some or all parts of this script were written with ChatGPT ##
## assistance. The author accepts full responsibility for the entirety of this script. This disclaimer, ##
## in-line comments, and all other documentation are written by the author and without the aid of       ##
## ChatGPT or a similar AI/chatbot.                                                                     ##
##########################################################################################################

### ABOUT ###
  # This script builds off the full_qiime2_pipeline.
  # This version was modified for highthroughput processing and for analyzing "guts" and "reproductive" samples separately.

### USAGE ###
  # ./alphas.sh [group] [feature_table_name_ONLY]
  # Do NOT put the extension when calling the feature table file!
  # Run script from main anlysis directory.
  
#!/bin/bash

group="$1"
features="$2"

# Observed features 
qiime diversity alpha-group-significance \
  --i-alpha-diversity "core_diversity_${group}/observed_features_vector.qza" \
  --m-metadata-file mappingv2.txt \
  --o-visualization "core_diversity_${group}/observed_features_vector.qzv"

# Faith's PD
qiime diversity alpha-group-significance \
  --i-alpha-diversity "core_diversity_${group}/faith_pd_vector.qza" \
  --m-metadata-file mappingv2.txt \
  --o-visualization "core_diversity_${group}/faith_pd-group-significance.qzv"

# Evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity "core_diversity_${group}/evenness_vector.qza" \
  --m-metadata-file mappingv2.txt \
  --o-visualization "core_diversity_${group}/evenness-group-significance.qzv"

# Shannon
qiime diversity alpha \
  --i-table "${features}.qza" \
  --p-metric shannon \
  --o-alpha-diversity "core_diversity_${group}/shannon_vector.qza"
  
qiime diversity alpha-group-significance \
  --i-alpha-diversity "core_diversity_${group}/shannon_vector.qza" \
  --m-metadata-file mappingv2.txt \
  --o-visualization "core_diversity_${group}/shannon-group-significance.qzv"

### END OF SCRIPT ###