#!/bin/bash

# This script will run data through the WFABC pipeline

# set some params
V1016I="./V1016I_byMo_sel.txt"
F1534C="./F1534C_byMo_sel.txt"

twoNe=1000
min_s=0.0
max_s=5.0
min_h=0
max_h=1

# step 1 - estimating s
../programs/WFABC_v1.1/binaries/Linux/wfabc_1 -nboots 0 $V1016I
../programs/WFABC_v1.1/binaries/Linux/wfabc_1 -nboots 0 $F1534C

# # step 2 - estimating h
../programs/WFABC_v1.1/binaries/Linux/wfabc_2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $V1016I
../programs/WFABC_v1.1/binaries/Linux/wfabc_2 \
  -fixed_N $twoNe -min_s $min_s -max_s $max_s \
  -min_h $min_h -max_h $max_h \
  $F1534C

