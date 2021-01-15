# IQTmosq
Code to analyze population data for Aedes aegypti from Iquitos, Peru

## Notes about WFABC Analysis  
** For this analysis the user will need to set their operating system in /WFABC/scripts/2_runWFABC.sh for the analysis to complete properly. Options are (case-sensitive): Linux, Mac, or Windows.
** The WFABC analysis contained within this repo utilizes a previously developed program WFABC v1.1 developed by Foll et al. (2014 and 2015).  
  + Link to paper: http://jjensenlab.org/wp-content/uploads/2016/02/Foll_Shim_Jensen_2014.pdf  
  + Link to program: http://jjensenlab.org/software  
  
=======
Code to analyze population data for Aedes aegypti from IQT

## Dependencies
* R (circa 3.5)
- see `setup.R` for required R libraries
- tinytex install
* Linux? (WFABC)

## Key files
* `make.R`: run analysis & initialize figures
* `build.R`: render figures usinguse output of `make.R` to re
