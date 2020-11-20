#!/bin/bash

# This script will run data through the WFABC pipeline

# # change directory
# cd ~/Dropbox/GouldLab/Project_mosquito/Database

# ### Run the WFABC analysis
# # change directory
# cd ./WFABC

################################
echo "Set some parameters"
################################
data1="./WFABC_popgensim_Fit0.11-.07-.08.txt"

twoNe=1000
min_s=-1.0
max_s=1.0
min_h=0
max_h=1

################################
echo "Estimating s"
################################

# step 1 - estimating s
~/src/WFABC_v1.1/binaries/Linux/wfabc_1 -Xmx1024m -nboots 0 $data1

# ################################
# echo "Estimating h"
# ################################
#
# # step 2 - estimating h
# ~/src/WFABC_v1.1/binaries/Linux/wfabc_2 \
#   -fixed_N $twoNe -min_s $min_s -max_s $max_s \
#   -min_h $min_h -max_h $max_h \
#   $data1
