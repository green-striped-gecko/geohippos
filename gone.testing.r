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
# 
# extract_gone <- function(df) {
#   res_epos <- data.frame()
#   for (i in 1:nrow(df)) {
#     new_data <- cbind( df[i,c(2:8, 13, 14)], df$GONE_01[[i]]);
#     res_epos <- rbind(res_epos, new_data);
#   }
#   return(res_epos);
# }
# 
# res <- extract_gone(all.gone)
# 
# gp <- ggplot(data = res, aes(x = Generation, y = Geometric_mean, colour = as.factor(pop_init):as.factor(crash_prop))) +
#   geom_line() +
#   geom_vline(xintercept = res$ts[1], colour = "blue") +
#   geom_vline(xintercept = (res$ts[1] - res$tl[1]), colour = "black") +
#   ggtitle(paste0( "Sample size", res$ss[1],
#                  ", start: " , res$ts[1] , "ybp" ,
#                  ", length: " , res$tl[1] , "yrs" )) +
#   theme_minimal() +
#   xlim(0, 1500) +
#   ylim(0, 600) +
#   labs(colour = "pop init: crash %") +
#   facet_wrap(model ~ maf, nrow = 3) 
# gp
all.gone
