#===============================================
#           Test Method Settings - EPOS
#===============================================

library(dplyr)
library(tidyverse)
library(geohippos)
library(dartR)
library(tictoc)


#Read in dataframe 

eposdf <- read.csv("ss100df.csv")
# 
 eposdf$filename <- str_replace(eposdf$filename, "c:/temp", "/data/scratch/isobel/vcf")
# 
# eposdf

#some checks on the input file
# nLoc(gls)
# nInd(gls)
# pop(gls)<- rep("A", nInd(gls))
# sfs <-gl.sfs(gls)

#check os to find correct binaries
os <- tolower(Sys.info()['sysname'])
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate


##############Epos testing################
settings <- expand_grid(minbin = 1:2, greedy = c("", " -E 2"), .name_repair = "minimal")
settings
test.epos <- as_tibble(crossing(eposdf, settings))
test.epos



###convert all to gls and run epos
library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")
test.epos$gls <- mclapply(1:nrow(test.epos), function(x) {
  out <- gl.read.vcf(test.epos$filename[[x]])
  return(out)

}, mc.cores=20)


#==================Run EPOS for all rows================================

test.epos$eposout <- mclapply(1:nrow(test.epos), function(x) {
  out <- gl.epos(test.epos$gls[[x]], epos.path = paste0("./binaries/epos/",os),l = L, u=mu, boot=20, minbinsize = test.epos$minbin[[x]], other.options = test.epos$greedy[[x]])
  return(out)
}, mc.cores = 20)


#================Extract loci for each run====================================
library(dartR)

for (i in 1:nrow(test.epos)) {
  test.epos$loci[[i]] <- nLoc(test.epos$gls[[i]])
}



#extract loci data
loci.labels <- curdata %>% group_by(pop_init, crash_prop) %>%
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$pop_init, loci.labels$crash_prop), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL)

gp <- ggplot(data = curdata, aes(x = X.Time, y = Median, colour = as.factor(pop_init):as.factor(crash_prop))) +
  geom_line() +
  geom_vline(xintercept = curdata$ts[1], colour = "blue") +
  geom_vline(xintercept = (curdata$ts[1] - curdata$tl[1]), colour = "black") +
  ggtitle(paste0("model: " , curdata$model[1] ,
                   ", start: " , curdata$ts[1] , "ybp" ,
                   ", length: " , curdata$tl[1] , "yrs" )) +
  theme_minimal() +
  xlim(0, 250) +
  ylim(0, 600) +
  labs(colour = "pop init: crash %") +
  facet_grid(greedy ~ minbin)

grid.arrange(gp, tb,heights = c(10, 1), widths = c(10, 1))
table(curdata$loci, curdata$pop_init + curdata$crash_prop)

n <- test.epos %>% filter(model == "bottle", ss == 30)
n$loci
