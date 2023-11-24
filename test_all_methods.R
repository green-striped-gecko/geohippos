#################################################################
#                     RUN ALL - TESTING                         #
#################################################################


####load necessary packages
library(dplyr)
library(tidyverse)
library(geohippos)
library(dartR)
library(parallel)
library(tictoc)

###prep data for all methods
tr_03 <- read.csv("ss100df.csv")
tr_03

runcode <- "tr_03"

df  <- tr_03
df$filename <- str_replace(df$filename, "c:/temp", "/data/scratch/isobel/vcf")

##generate gls files for all
setwd("~/R/geohippos")
df$gls <- mclapply(1:nrow(df), function(x) {
  out <- gl.read.vcf(df$filename[[x]])
  return(out)
}, mc.cores=20)



##generate loci numbers for all
for (i in 1:nrow(df)) {
  df$loci[[i]] <- nLoc(df$gls[[i]])
}



#Separate gls and loci columns
dfgls <- df %>% select(runnumb, gls, loci)


####specify settings for each method:
df.cross <- df[1:10]



#check os to find correct binaries
os <- tolower(Sys.info()['sysname'])
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate



#===================================#
#           Epos testing            #
#===================================#

#Define settings to test
ep.settings <- expand_grid(minbin = 1:2, greedy = c("", " -E 2", " -E 5"), .name_repair = "minimal")

#Generate all combinations of settings
test.epos <- as_tibble(crossing(df.cross, ep.settings))

#Read gls and loci information back in
test.epos <- left_join(test.epos, dfgls, by = "runnumb")

#Run epos (select for background job)
library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")
test.epos$eposout  <- mclapply(1:nrow(test.epos), function(x) {
  out <- gl.epos(test.epos$gls[[x]], epos.path = paste0("./binaries/epos/",os),l = L, u=mu, boot=20, minbinsize = test.epos$minbin[[x]], other.options = test.epos$greedy[[x]])
  return(out)
}, mc.cores = 20)

colnames(test.epos[ncol(test.epos)]) <- paste0("eposout_", runcode)



#===================================#
#      StairwayPlot testing         #
#===================================#

#Define settings to test
sw.settings <- expand_grid(minbin = c(1:3), breakpoints = c(4, 5, 8), .name_repair = "minimal")


#Generate all combinations of settings
test.sw <- as_tibble(crossing(df.cross, sw.settings))

#Read gls and loci information back in
test.sw <- left_join(test.sw, dfgls, by = "runnumb")

all.test.sw <- test.sw

test.sw <- test.sw[1:5,]

test.sw$stairway <- mclapply(1:nrow(test.sw), function(x) {
  out <- gl.stairway2(test.sw$gls[[x]], verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize = test.sw$minbin[[x]], cleanup = T, nrand = test.sw$breakpoints[[x]])
  saveRDS(out, file =  file.path("/data/scratch/isobel/results",paste0(runcode, test.sw$runnumb[x],".rds")))
  return(out)
}, mc.cores = 20)
