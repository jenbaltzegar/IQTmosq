#!/bin/bash

# This script will run data through the WFABC pipeline

# set some params
V1016I="./V1016I_byMo_sel.txt"
F1534C="./F1534C_byMo_sel.txt"

twoNe=1000
min_s=-1.0
max_s=0
min_h=0
max_h=1

# choose Mac, Linux, or Windows
myos=Linux
# path=paste0("../programs/WFABC_v1.1/binaries/", myos, "/")
bin1=../programs/WFABC_v1.1/binaries/$myos/wfabc_1
bin2=../programs/WFABC_v1.1/binaries/$myos/wfabc_2

# step 1 - estimating s
$bin1 -nboots 0 $V1016I
$bin1 -nboots 0 $F1534C

# # step 2 - estimating h
$bin2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $V1016I
$bin2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $F1534C

