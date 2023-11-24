#===============================================
#           Test Method Settings - GONE
#===============================================

library(dplyr)
library(tidyverse)
library(geohippos)
library(dartR)
library(tictoc)
library(parallel)

#Read in dataframe 
#df <- as.data.frame(read.csv("./maindf.csv", stringsAsFactors = F))

gdf <- as_tibble(df)

##change directory for filenames##
gdf$filename <- str_replace(gdf$filename, "c:/temp", "/data/scratch/isobel/vcf")


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
library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")


test.g$GONE_10 <- mclapply(1:nrow(test.g), function(x) {
  out <- gl.gone(test.g$gls[[x]],gone.path = paste0("./binaries/gone/",os), )
  saveRDS(out, file = file.path('/data/scratch/isobel/results/tr-02-', paste0(test.g$runnumb[x], ".rds")))
  return(out)
}, mc.cores = 20)


#==================Extract loci for each run===================================

locs <- lapply(test.g$gls, function(x) {
  out <- nLoc(x)
  return (out)
})

test.g$loci <- unlist(locs)

#=================Convert data to accessible format for plotting=============
# extract_epos <- function(df) {
#   res_epos <- data.frame()
#   for (i in 1:nrow(df)) {
#     new_data <- cbind( df[i,c(1:10, 14)], df$eposout1[[i]]);
#     res_epos <- rbind(res_epos, new_data);
#   }
#   return(res_epos);
# }

outdf <- df_extract_output(test.g, 13, c(1:9, 11, 12))
fname <- "gonerun22nov1"
write_csv(x = outdf, file = paste0("/data/scratch/isobel/results/", fname, ".csv"), col_names = T)

curdata <- stairway001.testdf %>% filter(model == "decline")

test.g

test.g.long <- test.g %>%  gather(key = "maf", value = "Ne", -rep, -pop_init, -tl, -ts, -ss, -model, -runnumb, -filename, -gls, -loci)

res.g <- df_extract_output(test.g.long, 12, 1:8)


sapply(1:nrow(test.g.long), function(x) {
  
  df <- data.frame(test.g.long$Ne[[x]])
  
  dummy <- test.g.long[x,c(1:8,10:11)]
  df <- rbind(df,dummy)
  return(df)
})


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