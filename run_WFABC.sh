#!/bin/bash

# This script will run data through the WFABC pipeline

# change directory
cd ~/Dropbox/GouldLab/Project_mosquito/Database
mkdir -p ./WFABC

### Create input files
Rscript ./R_Scripts/IQTmosq/run_WFABC_createDF.R

### Run the WFABC analysis
# change directory
cd ./WFABC

# set some params
data1="./multiple_loci.txt"
# data2="./V1016I_underSelection.txt"
# data3="./F1534C_underSelection.txt"
# data4="./V1016I_byMo.txt"
# data5="./F1534C_byMo.txt"
data6="./multiple_loci_byMo.txt"
data7="./multiple_loci_byMo_sel.txt"

twoNe=1000
min_s=-1.0
max_s=1.0
min_h=0
max_h=1

# step 1 - estimating s
~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data1
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data2
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data3
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data4
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data5
~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data6
~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 -nboots 0 $data7

# # step 2 - estimating h
~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $data1
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
#   -fixed_N $twoNe -min_s $min_s -max_s $max_s \
#   -min_h $min_h -max_h $max_h \
#   $data2
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
#   -fixed_N $twoNe -min_s $min_s -max_s $max_s \
#   -min_h $min_h -max_h $max_h \
#   $data3
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
#   -fixed_N $twoNe -min_s $min_s -max_s $max_s \
#   -min_h $min_h -max_h $max_h \
#   $data4
# ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
#   -fixed_N $twoNe -min_s $min_s -max_s $max_s \
#   -min_h $min_h -max_h $max_h \
#   $data5
~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $data6
~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $data7
  
# # Run R script to create plots
# Rscript ../R_Scripts/IQTmosq/run_WFABC_analyze.R


# Looping over multiple files --------------
# # set some params
# data=("./WFABC_multiple_loci.txt" "./WFABC_V1016I_underSelection.txt" "./WFABC_F1534C_underSelection.txt")
# twoNe=1000
# min_s=-1.0
# max_s=1.0
# min_h=0
# max_h=1
#
# # step 1
# for i in "${data[@]}"; do
#   ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_1 \
#     -nboots 0 \
#     $data
# done
#
# # step 2 - estimating h
# for i in "${data[@]}"; do
#   ~/../../Applications/WFABC_v1.1/binaries/Mac/wfabc_2 \
#     -fixed_N $twoNe -min_s $min_s -max_s $max_s \
#     -min_h $min_h -max_h $max_h \
#     $data
# done
