# Load data for analyses

# Isolation log contains all samples isolated
# Note: make sure that IsolationLog_IQT-Mosq.csv file has the date column in the MM-DD-YYYY format
kdrData <- read.csv("./data/IsolationLog_IQT-Mosq.csv")

# masterdec.csv contains GPS and neighborhood info
masterdec <- read.csv("./data/masterdec.csv")

# These files contain genotype information
merged1016_rep1 <- read.csv("./data/MeltCurve_1016_rep1.csv")
merged1016_rep2 <- read.csv("./data/MeltCurve_1016_rep2.csv")
merged1534_rep1 <- read.csv("./data/MeltCurve_1534_rep1.csv")
merged1534_rep2 <- read.csv("./data/MeltCurve_1534_rep2.csv")
merged410_rep1 <- read.csv("./data/MeltCurve_410_rep1.csv")
merged410_rep2 <- read.csv("./data/MeltCurve_410_rep2.csv")

# exptZone contains zone information for 2013 & 2014 experiments
exptZone <- read.csv("./data/Location_Zone.csv")

