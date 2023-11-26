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
library(future)

###prep data for all methods
tr_03 <- read.csv("ss100df.csv")
tr_03

runcode <- "tr_03"

df  <- tr_03
#df$filename <- str_replace(df$filename, "c:/temp", "/data/scratch/isobel/vcf")
df$filename <- str_replace(df$filename, "runnum:", "runnum")

##generate gls files for all
setwd("~/R/geohippos")

##MCLAPPLY gls file generation####
df$gls <- mclapply(1:nrow(df), function(x) {
  out <- gl.read.vcf(df$filename[[x]])
  return(out)
}, mc.cores=20)

##generate gls files - non-parallel######
for (i in 1:nrow(df)) {
  df$gls[[i]] <- gl.read.vcf(df$filename[[i]])
}

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

################################################################################
#============             SETUP SETTING COMBINATIONS              =============#
################################################################################

#====== EPOS =========#
#Define settings to test
ep.settings <- expand_grid(minbin = 1:2, greedy = c("", " -E 2", " -E 5"), .name_repair = "minimal")

#Generate all combinations of settings
test.epos <- as_tibble(crossing(df.cross, ep.settings))

#Read gls and loci information back in
test.epos <- left_join(test.epos, dfgls, by = "runnumb")

#====== STAIRWAY =======#
#Define settings to test
sw.settings <- expand_grid(minbin = c(1:2), breakpoints = c(4, 6), .name_repair = "minimal")

#Generate all combinations of settings
test.sw <- as_tibble(crossing(df.cross, sw.settings))

#Read gls and loci information back in
test.sw <- left_join(test.sw, dfgls, by = "runnumb")

#===== GONE ==========#
test.g <- df


################################################################################
#============             RUN METHODS WITH MCLAPPLY               =============#
################################################################################


#===================================#
#           Epos testing            #
#===================================#
#Run epos (select for background job)
library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")
test.epos$eposout  <- mclapply(1:nrow(test.epos), function(x) {
  out <- gl.epos(test.epos$gls[[x]], epos.path = paste0("./binaries/epos/",os),l = L, u=mu, boot=20, minbinsize = test.epos$minbin[[x]], other.options = test.epos$greedy[[x]])
  saveRDS(out, file =  file.path("/data/scratch/isobel/results/", paste0("Epos_", runcode, test.epos$runnumb[x],".rds")))
  return(out)
}, mc.cores = 10)

colnames(test.epos[ncol(test.epos)]) <- paste0("eposout_", runcode)

#===================================#
#      StairwayPlot testing         #
#===================================#

all.test.sw <- test.sw

test.sw <- test.sw[1:5,]

test.sw$stairway <- mclapply(1:nrow(test.sw), function(x) {
  out <- gl.stairway2(test.sw$gls[[x]], verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize = test.sw$minbin[[x]], cleanup = T, nrand = test.sw$breakpoints[[x]])
  saveRDS(out, file =  file.path("/data/scratch/isobel/results/", paste0("Stairway_", runcode, test.sw$runnumb[x],".rds")))
  return(out)
}, mc.cores = 10)

#===================================#
#           GONE testing            #
#===================================#

library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")

####CURRENT SETTINGS#####
##THREADS = 5
##MAF = 0.0

test.g$GONE_00 <- mclapply(1:nrow(test.g), function(x) {
  out <- gl.gone(test.g$gls[[x]],gone.path = paste0("./binaries/gone/",os), )
  saveRDS(out, file = file.path('/data/scratch/isobel/results/', paste0("GONE_", runcode, test.g$runnumb[x], ".rds")))
  return(out)
}, mc.cores = 10)


################################################################################
#==============       RUN METHODS WITH PARSAPPLY                  =============#
################################################################################

library(doParallel)
cl <- makeCluster(6) 
ret <- clusterEvalQ(cl, library("geohippos"))
ret <- clusterEvalQ(cl, library("future"))
ret <- clusterExport(cl, varlist = c("test.sw","L","mu")) 


#============= EPOS ================#
epos_sims <- parSapply(cl, 1:nrow(test.epos), function(x) {
  
  #load data from fire.  
  res_epos <- NULL
  #epos  
  res_epos <- gl.epos(test.epos$gls[[x]], epos.path = paste0("./binaries/epos/",os),l = L, u=mu, boot=20, minbinsize = test.epos$minbin[[x]], other.options = test.epos$greedy[[x]])

  #save all three in one output
  out <- list(epos=res_epos)
  saveRDS(out, file =  file.path("C:/Users/Isobel/Desktop/Honours/RDS_files/", paste0("Epos_", runcode, test.epos$runnumb[x],".rds"))) 
  return(TRUE)
})

stopCluster(cl)



stair_sims <- parSapply(cl, 1:nrow(test.sw), function(x) {
  
  #load data from fire.  
  res_stair <- NULL

  res_stair <- gl.stairway2(test.sw$gls[[x]], verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=1, L=L, minbinsize = test.sw$minbin[[x]], cleanup = T, nrand = test.sw$breakpoints[[x]])
  
  #save all three in one output
  out <- res_stair
  #saveRDS(out, file =  file.path("~/RDS_files/", paste0("Stairway_", runcode, test.sw$runnumb[x],".rds")))
  return(TRUE)
})

stopCluster(cl)

gone_sims <- parSapply(cl, 1:nrow(test.g), function(x) {
  
  #load data from fire.  
  res_gone <- NULL
  
  res_gone <- gl.gone(test.g$gls[[x]],gone.path = paste0("./binaries/gone/",os), )
  
  #save all three in one output
  out <- list(gone=res_gone)
  saveRDS(out, file =  file.path("C:/Users/Isobel/Desktop/Honours/RDS_files/", paste0("GONE_", runcode, test.sw$runnumb[x],".rds")))
  return(TRUE)
})

stopCluster(cl)


