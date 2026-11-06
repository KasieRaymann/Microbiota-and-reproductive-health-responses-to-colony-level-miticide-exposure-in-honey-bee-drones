# QIIME2 Pipeline | Perform stats on beta diversity results (PERMANOVA)

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
  # This script will run through PERMANOVAS for each beta diversity metric (Bray-Curtis, Weighted Unifrac, Unweighted Unfrac)
  # These will be run only for the metadata column/variable you call for.
  # This script will also create child directories under /core_diversity/ called /permanovas/ and then a child directory under that for the /metadata_column/ you called.
  # You can only run this script after having generated the [metric]_distance_matrix.qza files.

### USAGE ###
  # ./permanova.sh [group] [metadata_column] [metadata_column] [metadata_column] [...]
  # Run script from main analysis directory.

#!/bin/bash

group="$1"
for metadata_column in "${@:2}"
do
  mkdir -p "core_diversity_${group}/permanovas/${metadata_column}"     # Will create a directory if it does not exist; no return error if directory already exists

  for metric in bray_curtis weighted_unifrac unweighted_unifrac        # Can be changed to run through other metrics if desired.       
  do
    distance_matrix="core_diversity_${group}/${metric}_distance_matrix.qza"
    output="core_diversity_${group}/permanovas/${metadata_column}/${metric}_distance_matrix_${metadata_column}.qzv"   # Will output .qzvs into core_diversity/PERMANOVs/metadata_column_name/

    qiime diversity beta-group-significance \
      --i-distance-matrix "$distance_matrix" \
      --m-metadata-file mappingv2.txt \
      --m-metadata-column "$metadata_column" \
      --o-visualization "$output" \
      --p-pairwise
  done

done

### END OF SCRIPT ###