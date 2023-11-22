#===============================================
#           Test Method Settings - Stairways
#===============================================

library(dplyr)
library(tidyverse)
library(geohippos)
library(tictoc)

#Read in dataframe 
df <- as.data.frame(read.csv("./maindf.csv", stringsAsFactors = F))

swdf <- as_tibble(df)

##change directory for filenames##
swdf$filename <- str_replace(swdf$filename, "c:/temp", "~/R/vcf.testing")

#check os to find correct binaries
os <- tolower(Sys.info()['sysname'])
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate


#===================Set parameter variations to test ===========================
settings <- expand_grid(minbin = c(1:2), breakpoints = c(4))
settings
test.sw <- crossing(swdf, settings)

reserve <- test.sw %>% filter(minbin == 2)
test.sw <- test.sw %>% filter(minbin == 1)
test.sw <- reserve
#===================convert to gls for all runs ==============================

test.sw$gls <- mclapply(1:nrow(test.sw), function(x) {
  out <- gl.read.vcf(test.sw$filename[x])
  return(out)
  
}, mc.cores=20)

for (i in 1:nrow(test.sw)) {
  test.sw$gls[[i]] <- gl.read.vcf(test.sw$filename[[i]])
}


#==================Run StairwayPlot for all rows================================
test.sw$stairway <- mclapply(1:nrow(test.sw), function(x) {

  out <- gl.stairway2(test.sw$gls[[x]], verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize = test.sw$minbin[[x]], cleanup = T, nrand = test.sw$breakpoints[[x]])
  # newname <- str_replace(test.sw$filename[[x]], "vcf$", "csv")
  # outname <- str_replace(newname, "^~/R/vcf.testing/", "/home/R/geohippos/filesout/")
  # write_csv(out, file = outname)
  return(out)
    
}, mc.cores = 20)

stairout <- test.sw$stairway


#==================Extract loci for each run===================================

for (i in 1:nrow(test.sw)) {
  test.sw$loci[i] <- nLoc(test.sw$gls[[i]])
}
test.sw$loci

#==================Extract Stairway Data into usable dataframe=================
extract_inf <- function(df1, df2) {
  if(nrow(df1) != length(df2)) {
    print("Error: Output data length not equal to input data length")
    return()
  }
  res_stair <- data.frame();
  for (i in nrow(df1)) {
    new_data <- cbind(df1[i, c(2:9, 11:12, 15)], df2[i]);
    res_stair <- rbind(res_stair, new_data);
  }
  return(res_stair)
}

res <- extract_inf(test.sw, stairout)
res

mod = "decline"

curdata <- res %>% filter(model == mod)

#=================Extract loci data

write.csv(x = res, file = "~/stairwayrundata.csv", col.names = T)

res
class(res)
write_csv(x = res, file = "~/stairwayplotdata.csv", col_names = T)

