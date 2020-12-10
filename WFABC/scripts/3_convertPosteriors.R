#!/usr/bin/env Rscript

### Convert WFABC posteriors to match Michael Vella's scale

# Do the conversion
convertPosteriors.onelocus("V1016I_byMo_sel_posterior_s.txt", "V1016I_byMo_sel_posterior_h.txt")
convertPosteriors.onelocus("F1534C_byMo_sel_posterior_s.txt", "F1534C_byMo_sel_posterior_h.txt")
