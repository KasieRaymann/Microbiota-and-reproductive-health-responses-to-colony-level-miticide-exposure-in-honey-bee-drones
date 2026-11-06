# QIIME2R | Alpha diversity signficnance extractions for ggsignif

# Author:        Layla Jeries, with some ChatGPT assistance
# Date created:  2026-09-11
# Last modified: 2026-11-06

########################################## AI USAGE DISCLAIMER ###########################################
## This disclaimer has been included because some or all parts of this script were written with ChatGPT ##
## assistance. The author accepts full responsibility for the entirety of this script. This disclaimer, ##
## in-line comments, and all other documentation are written by the author and without the aid of       ##
## ChatGPT or a similar AI/chatbot.                                                                     ##
##########################################################################################################

### ABOUT ###
    # This script will temporarily unzip the [alpha metric]-group-significance.qzv files to then search and copy over the kruskal-wallis-pairwise.csv
    # These .csvs have the q-value significance scores for plotting in R using ggsignif.

### USAGE ###
    # ./alphaextract.sh
    # Please call from WITHIN the core_diversity directory you want to analyze.
    # All the temporary unzipped directories get automatically deleted after finshing the loop (see line 40)

#!/bin/bash

for metric in shannon faith_pd
do 
    mkdir -p "sig_csvs_${metric}"     # Will create a directory if does not exist; no return error if directory already exists
done

for metric in shannon faith_pd
do 
    unzip -q "${metric}-group-significance.qzv" -d "${metric}_temp"
    
    find "${metric}_temp"/*/data/ -type f -name 'kruskal-wallis-pairwise-C_T.csv' -exec cp -v {} sig_csvs_${metric} \;
    find "${metric}_temp"/*/data/ -type f -name 'kruskal-wallis-pairwise-Best_Worst.csv' -exec cp -v {} sig_csvs_${metric} \;
    find "${metric}_temp"/*/data/ -type f -name 'kruskal-wallis-pairwise-CT_BestWorst.csv' -exec cp -v {} sig_csvs_${metric} \;

    rm -r "${metric}_temp"
done

### END OF SCRIPT ###