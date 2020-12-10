#!/bin/bash

# This script will run data through the WFABC pipeline

# change directory
cd ./WFABC

# make directories
mkdir -p ./results

echo "Creating input files..."
Rscript ./scripts/1_createDF.R

echo "Running the WFABC analysis..."
./scripts/2_runWFABC.sh

echo "Converting the Posteriors..."
Rscript ./scripts/3_convertPosteriors.R

echo "Analyzing posteriors and creating plots..."
Rscript ./scripts/4_analyze_byMo_sel_converted.R
mv Rplots.pdf ./results/Rplots.pdf

# change directory
cd ../