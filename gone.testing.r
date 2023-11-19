#===============================================
#           Test Method Settings - GONE
#===============================================

library(dplyr)
library(tidyverse)
library(geohippos)
library(dartR)
library(tictoc)

#Read in dataframe 
df <- as.data.frame(read.csv("./maindf.csv", stringsAsFactors = F))

gdf <- as_tibble(df)

##change directory for filenames##
gdf$filename <- str_replace(gdf$filename, "c:/temp", "~/R/geohippos/vcf.testing")


#check os to find correct binaries
os <- tolower(Sys.info()['sysname'])
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate


#===================Set parameter variations to test ===========================

#test.g <- crossing(gdf, settings)
test.g <- as_tibble(gdf)

#===================convert to gls for all runs ==============================

test.g$gls <- mclapply(1:nrow(test.g), function(x) {
  out <- gl.read.vcf(test.g$filename[[x]])
  return(out)
  
}, mc.cores=20)


#==================Run GONE for all rows================================
test.g$GONE_01 <- mclapply(1:nrow(test.g), function(x) {
  out <- gl.gone(test.g$gls[[x]],gone.path = paste0("./binaries/gone/",os), )
  return(out)
}, mc.cores = 20)


#==================Extract loci for each run===================================

for (i in 1:nrow(test.g)) {
  test.g$loci[[i]] <- nLoc(test.g$gls[[i]])
}

gone.maf.0.02 <- test.g


#=================Convert data to accessible format for plotting=============
extract_epos <- function(df) {
  res_epos <- data.frame()
  for (i in 1:nrow(df)) {
    new_data <- cbind( df[i,c(1:10, 14)], df$eposout1[[i]]);
    res_epos <- rbind(res_epos, new_data);
  }
  return(res_epos);
}
