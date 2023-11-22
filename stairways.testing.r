#===============================================
#           Test Method Settings - Stairways
#===============================================

library(dplyr)
library(tidyverse)
library(geohippos)
library(dartR)
library(tictoc)

#Read in dataframe 
df <- as.data.frame(read.csv("./big samplesize df.csv", stringsAsFactors = F))

swdf <- as_tibble(df)

##change directory for filenames##
swdf$filename <- str_replace(swdf$filename, "c:/temp", "/data/scratch/isobel/vcf")

#check os to find correct binaries
os <- tolower(Sys.info()['sysname'])
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate


#===================Set parameter variations to test ===========================
settings <- expand_grid(minbin = c(1:3), breakpoints = c(4, 5, 8))
settings
test.sw <- crossing(swdf, settings)

#===================convert to gls for all runs ==============================

test.sw$gls <- mclapply(1:nrow(test.sw), function(x) {
  out <- gl.read.vcf(test.sw$filename[x])
  return(out)
  
}, mc.cores=20)



#==================Run StairwayPlot for all rows================================
test.sw$stairway <- mclapply(1:nrow(test.sw), function(x) {

  out <- gl.stairway2(test.sw$gls[[x]], verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize = test.sw$minbin[[x]], cleanup = T, nrand = test.sw$breakpoints[[x]])

  return(out)
    
}, mc.cores = 20)


#==================Extract loci for each run===================================
locs <- lapply(test.sw.all$gls, function(x) {
  out <- nLoc(x)
  return (out)
})

test.sw.all$loci <- unlist(locs)


#=================Extract loci data
# 
# write.csv(x = res, file = "~/stairwayrundata.csv", col.names = T)
# 
# res
# class(res)
# write_csv(x = res, file = "~/stairwayplotdata.csv", col_names = T)
# 

###convert existing outputs to long form

####Dataset Stairway001.test
#pop_init: 500, 200
#crash_prop: 0.5, 0.1
#tl: 100
#ts: 200
#ss: 20
#minbinsize: 1,2
#breakpoints: 4,8,10
outdf <- df_extract_output(test.sw.all, 14, 2:12)
fname <- "stairwaytest001"
write_csv(x = outdf, file = paste0("/data/scratch/isobel/results/", fname, ".csv"), col_names = T)

curdata <- stairway001.testdf %>% filter(model == "decline")


gp <- ggplot(data = curdata, aes(x = year, y = Ne_median, colour = as.factor(pop_init), linetype = as.factor(crash_prop))) +
  geom_line() +
  geom_vline(xintercept = curdata$ts[1], colour = "blue") +
  geom_vline(xintercept = (curdata$ts[1] - curdata$tl[1]), colour = "black") +
  ggtitle(paste0("Stairway Plot -- model: " , curdata$model[1] ,
                 ", start: " , curdata$ts[1] , "ybp" ,
                 ", length: " , curdata$tl[1] , "yrs" ,
                 ", sample size = ", curdata$ss[1])) +
  theme_minimal() +
  xlim(0, 1000) +
  ylim(0, 600) +
  labs(colour = "pop init: crash %") +
  facet_grid(breakpoints~minbin, labeller = label_both)
gp

loci.labels <- test.sw.all %>% group_by(pop_init, crash_prop) %>%
  summarise(loc = mean(loci))
loc.long <- data.frame(key = paste(loci.labels$pop_init, loci.labels$crash_prop), loci = loci.labels$loc)
loc.table <- t(loc.long)
tb <- tableGrob(loc.table, rows = NULL)

grid.arrange(gp,tb, heights = c(10, 2))
