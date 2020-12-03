# This file takes output from meltcurve analysis folders and concatenates it into six files
# Creates on file for each meltcurve locus + replicate combination
# Started: 6/7/17

#######
# If you receive the following error -- 
# "Error in file(file, "rt") : cannot open the connection' 
# -- make sure all analysis files are available and named properly in MeltCurve folder"
#######

# 0 ################################################################################
### load required library
library(gdata)


# Specify 4 folder names to loop over in setwd() loop
wd <- c("MeltCurve_1016_rep1", "MeltCurve_1016_rep2", "MeltCurve_1534_rep1", "MeltCurve_1534_rep2"
        , "MeltCurve_410_rep1", "MeltCurve_410_rep2")

# Begin setwd() loop
for (d in 1:length(wd)) {
  setwd(file.path("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Data_Mosquito/MeltCurve/MeltCurve_results/", wd[[d]]))
  num_files <- length(list.files())
  
  ### Initialize objects for rbind() loop
  # First part of file path
  str_First <- file.path("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Data_Mosquito/MeltCurve/MeltCurve_results/"
                    , wd[[d]], "MeltCurve_")
  # Initial file index 
  i <- 1
  # Last part of file path
  str_Last <- "/Analysis/PlateMap_Genotype_Zone.csv"
  
  # Create empty dataframe
  df <- data.frame()
  
  for (i in 1:num_files){
    # Paste parts of file path together and loop over i
    str_Tmp = paste(str_First, as.character(i), str_Last ,sep = "")
    # Read in csv files based on file path specified above
    df_New = read.csv(str_Tmp, header = T, sep = ",")
    # Concatenate csv files into one dataframe
    df = rbind(df, df_New)
  } # end file loop
  
  write.csv(df, file = paste("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/"
                           , as.character(wd[d]), ".csv", sep = "")
            , row.names = F)  
} # end setwd() loop
  