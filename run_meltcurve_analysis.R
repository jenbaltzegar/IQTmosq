# This script will apply meltcurve functions to kdrData and create objects required for downstream analysis
# requires run_prepData.R to create kdrData object

### Parse data for allele frequency analysis -----
# Parse data by year and/or month to use for allele frequency analysis

# Select based on year -----
dat.mosq <- sqldf("Select haplotype, strftime('%Y', newDate) as year from kdrData where newDate between '2000-01-01' and '2017-12-31'")
# For haplotypes at all years - requires function_mc.haps.R -----
mc.haps.yr <- (
    dat.mosq
    ##! note: year column is now first (instead of last)
    %>% group_by(year)
    %>% group_modify(~mc.haps(.x))
)

# Select based on month from each year -----
# Select based on month from all 2000 data
mc.dat <- sqldf("Select newDate, V1016I_converted, F1534C_converted,  strftime('%Y-%m', newDate) as YrMon from kdrData where newDate between '2000-01-01' and '2017-12-31'") 
## make sure all yr.mon combos are present
yr.mon <- with(
    expand.grid(mon=1:12, yr=2000:2017), 
    data.frame(YrMon = sprintf('%d-%02d', yr, mon), month=mon, year=yr)
)
mc.dat <- (
    merge(yr.mon, mc.dat, all=T, by='YrMon')
    %>% group_by(YrMon)
)    

# For 1016 at all month+year combos -----
# mc.1016.byMo - For 1016 locus by month & year -----
mc.1016.byMo <- group_modify(mc.dat, ~mc.1016(.x))
# write df to csv
#write.csv(mc.1016.byMo, file = "./data/mc.1016.byMo.csv", row.names = FALSE)
# For 1534 at  all month+year combos -----
# mc.1534.byMo - For 1534 locus by month & year -----
mc.1534.byMo <- group_modify(mc.dat, ~mc.1534(.x))
# write df to csv
#write.csv(mc.1534.byMo, file = "./data/mc.1534.byMo.csv", row.names = FALSE)

##? TODO: check observation counts for 2014, e.g. March treatment
##  aa = subset(kdrData, project_code == 'treatment' & year(newDate)==2014 & month(newDate)==3)
dat.zone <- sqldf("
    Select V1016I_converted, F1534C_converted,  newDate, project_code as zone,
    strftime('%Y', newDate) as year,  strftime('%m', newDate) as month 
    from kdrData
    where project_code is not null
    and V1016I_converted is not null 
    and not V1016I_converted = 'error'
    and F1534C_converted is not null
    and not F1534C_converted = 'error'
    and (
        newDate between '2013-01-01' and '2013-10-31' 
        or newDate between '2014-01-01' and '2014-10-31'
    ) 
    and not (
        newDate > '2014-01-01' and project_code = 'other'
    )
")

## 
.id.x <- with(dat.zone, 
    year==2013 & zone != 'treatment'
)
## indentify 2013 "expanded buffer"
dat.zone <- within(dat.zone, {
    zone[.id.x] <- 'buffer'
    ## create logical column to mark them
    expanded <- F
    expanded[.id.x] <- T
})

# For 1016 - by zone, 2013 & 2014 -----
mc.1016.zone <- (
    dat.zone
    ##! note: year column is now first (instead of last)
    %>% group_by(year, month, zone)
    %>% group_modify(~mc.1016(.x))
)

#oct2014t <- sqldf("Select * from trt where newDate between '2014-10-01' and '2014-10-31'")

#oct2013b.x <- buff.x[which(buff.x$newDate >= "2013-10-01" & buff.x$newDate <= "2013-10-31" & !is.na(buff.x$project_code)),]
# Select based on 2014 month from buffer zone
#jan2014b <- sqldf("Select * from buff where newDate between '2014-01-01' and '2014-01-31'")

#mc.1016.t14 <- cbind(month, dftrt)
#mc.1016.b13.expanded <- cbind(month, dfbuff13.x)
#mc.1016.b14 <- cbind(month, dfbuff)

# View dataframes -----
# mc.haps.yr
# mc.1016.t13
# mc.1016.t14
# mc.1016.b13.expanded
# mc.1016.b14
# mc.1016.byMo
# mc.1534.byMo
