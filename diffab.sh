# QIIME2 Pipeline | Differential abundance

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
  # This script will run through all the differential abundance analyses for the metadata column you want and for what you want to compare to from within that column.
  # These will be run only for the comparison observations that you call, not every type within that metadata column.
  # This script will also create child directories under /core_diversity/ called /diff_abundance/ and then a child directory under that for the reference variable you called.
  # You can only run this script after having generated a collapsed taxon table; you must call its same level too.

### USAGE ###
  # ./diffab.sh [group] [taxon_level_number(1-7)] [metadata_column] [compare] [compare] [...]
  # You can list as many [compare] variables as you'd like, but make sure they match your metadata.
  # Run script from main anlysis directory.

#!/bin/bash

group="$1"

level="$2"
if [ "$level" -eq 1 ]; then
  taxon="domain"
elif [ "$level" -eq 2 ]; then
  taxon="phylum"
elif [ "$level" -eq 3 ]; then
  taxon="class"
elif [ "$level" -eq 4 ]; then
  taxon="order"
elif [ "$level" -eq 5 ]; then
  taxon="family"
elif [ "$level" -eq 6 ]; then
  taxon="genus"
elif [ "$level" -eq 7 ]; then
  taxon="species"
fi

metadata_column="$3"

mkdir -p "core_diversity_${group}/diff_abundance/${metadata_column}"

for compare in "${@:4}"
do

  qiime composition ancombc \
    --i-table "coll-table-l${level}${taxon}_${group}.qza" \
    --m-metadata-file "mappingv2.txt" \
    --p-formula "${metadata_column}" \
    --p-reference-levels "${metadata_column}::${compare}" \
    --o-differentials "core_diversity_${group}/diff_abundance/${metadata_column}/l${level}${taxon}_ancombc_${metadata_column}vs${compare}.qza"

qiime composition da-barplot \
    --i-data "core_diversity_${group}/diff_abundance/${metadata_column}/l${level}${taxon}_ancombc_${metadata_column}vs${compare}.qza" \
    --p-significance-threshold 0.05 \
    --p-level-delimiter ';' \
    --o-visualization "core_diversity_${group}/diff_abundance/${metadata_column}/l${level}${taxon}_ancombc_${metadata_column}vs${compare}.qzv"

done

### END OF SCRIPT ###