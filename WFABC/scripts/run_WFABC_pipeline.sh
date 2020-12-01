#!/bin/bash

# This script will run data through the WFABC pipeline

# change directory
cd ./WFABC

# make directories
mkdir -p ./results ./results/converted #./results/unconverted
# mkdir -p ./results/byYr ./results/byYr/converted ./results/byYr/unconverted

echo "Creating input files..."
Rscript ./scripts/1_createDF.R

echo "Running the WFABC analysis..."
./scripts/2_runWFABC.sh

echo "Converting the Posteriors..."
Rscript ./scripts/3_convertPosteriors.R

echo "Analyzing posteriors and creating plots..."
# Rscript ./scripts/4_analyze.R
# Rscript ./scripts/4_analyze_byMo.R
# Rscript ./scripts/4_analyze_underSelection.R
# mv Rplots.pdf ./results/byYr/unconverted/Rplots.pdf
# 
# Rscript ./scripts/4_analyze_underSelection_converted.R
# mv Rplots.pdf ./results/byYr/converted/Rplots.pdf
# 
# Rscript ./scripts/4_analyze_byMo_sel.R
# mv Rplots.pdf ./results/unconverted/Rplots.pdf
# 
Rscript ./scripts/4_analyze_byMo_sel_converted.R
mv Rplots.pdf ./results/converted/Rplots.pdf
