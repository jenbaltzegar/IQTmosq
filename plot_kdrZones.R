## depends on setup.R and load_data.R

## process data: 
## minimum observations per month
.min.obs <- 5
## helper function
mk.date <- function(yr, mon, day=1) as.Date(paste(yr, mon, day, sep='-'))

mc.1016.sub <-  (
    ## remove rows with too few observations
    subset(mc.1016.zone, n > .min.obs)
    ## add columns
    %.>% within(., {
        date <- mk.date(year, month)
        ## use factor for glm (ANOVA, not slope)
        fmonth <- factor(month)
        n.allele <- 2*n
    })
)

## separate model for each, results in list
kdr.mod <- dlply(mc.1016.sub, 'year', function(dat) {
    ## logistic regression of allele freqency
    mod <- glm(data=dat,
        freqR ~ fmonth*zone,
        family=binomial, weights=n.allele
    )
    ## compute contrasts
    emm <- emmeans(mod, spec=~zone|fmonth, type='response')
    contr <- as.data.frame(
        contrast(emm, method='trt.vs.ctrl', type='response')
    )
    list(mod=mod, emm=emm, contr=contr)
})

## pull out contrasts, combine years
kdr.mod.contr <- (
    ldply( kdr.mod, function(yr) yr$contr, .id='year')
    %.>% within(., {
        date <- mk.date(year, as.character(fmonth))
    })
)

## pull out glm CI
kdr.mod.ci <- (
    ldply( kdr.mod, function(yr) as.data.frame(yr$emm), .id='year')
    %.>% within(., {
        date <- mk.date(year, as.character(fmonth))
    })
)
mc.1016.sub <- merge(mc.1016.sub, kdr.mod.ci)  
        
kdrZones <- (
  ggplot(mc.1016.sub, 
    aes(x=date, y=freqR, color=zone, ymin=asymp.LCL, ymax=asymp.UCL, label=n)
  ) +
  facet_grid( ~ year, scale='free_x', space='free_x') +
  labs(x = "Month", y = "Frequency") +
  scale_color_manual(
    name = "Treatment", 
    values = c("treatment" = "red", "buffer" = "blue"), 
    labels = c("Spray Zone", "Buffer Zone"), 
    breaks=c("treatment","buffer")
  ) +
  # Add yellow background to represent spray periods
  # Spray date range: 4/29/13 (4.96) - 6/3/13 (6.1)
  geom_rect(
    aes(xmin=xmin, xmax=xmax),
    inherit.aes=F,
    data=data.frame(
        xmin=as.Date(c('2013-04-15', '2014-01-25', '2014-04-15')),
        xmax=as.Date(c('2013-05-15', '2014-02-05', '2014-05-15')),
        year=c(2013, 2014, 2014)
    ),
    ymin = -Inf, ymax = Inf, 
    fill="yellow", alpha = 0.4
  ) +
  geom_text(
    aes(x=date, label=label),
    inherit.aes=F,
    y=1.0, size = 5, 
    data=data.frame( 
        date=as.Date(c('2013-05-01', '2014-05-01', '2014-02-01')),
        label= c("Experimental", "Experimental", "City Wide"),
        year=c(2013, 2014, 2014)
    )
  ) +
  ## show p-values
  geom_text(
    ## all
    #data=kdr.mod.contr, 
    ## only small-ish p-values
    data=subset(kdr.mod.contr, p.value<0.1), 
    inherit.aes=F,
    nudge_x = -3,
    nudge_y = .05+ c(1,0)*0.08,
    mapping=aes(x=date, y=0.1, label=sprintf('p=%0.2g', p.value))
    ## override: just asterisk
    #label='*'
  )+
  #Add data
  geom_line() +
  ## width in units of days
  geom_point(size = 3) +
  geom_errorbar(width=3,  size = 0.7, alpha=0.8) +
  geom_text( 
    ## treatment above buffer
    aes(y=(0.0 + 0.07*(zone=='treatment'))),
    fontface=2
  ) +
  ## add "N=" in each panel
  geom_text(
    fontface=2, 
    inherit.aes=F,
    mapping=aes(x=x, y=y, label=label),
    data=data.frame(
        y=0.035, 
        x=as.Date(c('2013-03-20', '2013-12-15')), 
        label='N=', year=c('2013','2014')
    )
  )+
  ylim(c(0.0,1)) +
  #xlim(c(1,11)) +
  ## label by month (%B is full month)
  scale_x_date(
    date_breaks = '1 month', date_labels='%b',
    expand=expansion(add=c(10,15))
  ) +
  my_theme +
  theme(
    legend.position='top',
    legend.direction='horizontal',
  ) +
  guides(color=guide_legend(ncol=2, override.aes=list(size=3, linetype=0)))
)

# View plot
#kdrZones13
 
# Write plot to png
# ggsave(filename = paste0("figures/kdrZones/kdrZones_2013/kdrZones13_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')
#ggsave(filename = paste0("./data/kdrZones13_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')

# 
# # Write plot to pdf
# pdf(file = paste("figures/kdrZones/kdrZones_", Sys.Date(), ".pdf", sep = ""), 11, 8.5)
# print(kdrZones)
# dev.off()
