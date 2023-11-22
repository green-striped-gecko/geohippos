#===============================================
#           Test Method Settings
#===============================================
library(dplyr)
library(data.table)
library(tidyr)
library(grid)
#Read in dataframe 

eposdf <- df


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


##############Set up Epos testing parameters ################

settings <- expand_grid(minbin = 1:2, greedy = c(""), .name_repair = "minimal")
settings


no.gls <- eposdf %>% select(-gls)
test.epos <- crossing(no.gls, settings)
test.epos

#======================Extract gls files======================================
library(geohippos)
setwd("C:/Users/Isobel/Desktop/Honours/geohippos")



get_gls <- function(df) {
  for (i in 1:nrow(df)) {
  df$gls[[i]] <- gl.read.vcf(df$filename[[i]]);
  print(paste("Genlight ", i, " of ", nrow(df), "extracted."))
  }
  print("All genlights extracted!")
  return(df)
}

#Generate gls files
eposdf <- get_gls(eposdf)

gls.files <- eposdf %>% select(runnumb, gls)

#===============Join GLS files to main dataset=================================
all.epos <- left_join(test.epos, gls.files, join_by(runnumb))

##======================Setup parallel=======================================
library(parallel)
cl <- makeCluster(4)
clusterCall(cl, function(x) x^2, 10)
clusterEvalQ(cl, library(geohippos))
clusterEvalQ(cl, library(dartR))
clusterEvalQ(cl, library(adegenet))
clusterExport(cl, "small.eposdf", envir = .GlobalEnv)
clusterExport(cl, c("os", "L", "mu"), envir = .GlobalEnv)


#=============RUN EPOS IN PARALLEL==================================


epos.output <- parSapply(cl, 1:nrow(small.eposdf), function (x) {
  out = gl.epos(
    small.eposdf$gls[[x]], 
    epos.path = paste0("./binaries/epos/",os),
    l = L, u=mu, 
    boot=30, 
    minbinsize = 1)
  return(out)} )


test.epos

#================Extract loci for each run====================================

for (i in 1:nrow(test.epos)) {
  test.epos$loci[i] <- nLoc(test.epos$gls[[i]])
}

for (i in 1:nrow(res2211)) {
  res2211$loci[i] <- as.numeric(res2211$loci[i][[1]])
}


smallres <- extract_epos(small.eposdf, epos.output)

res2211 <- df_extract_output(test.epos, out_index = 13, keep = c(1:11,14))

bin = 1 #(1:4)
othop = "" #c("", " -E 2", " -E 5", " -E 10", " -E 20")
p0 <- 200
cp <- 0.1
mod = "decline"

curdata <- trialres %>% filter(model == mod)

cursims <- siminds %>% filter(model == mod)

#extract loci data
loci.labels <- res2211 %>% group_by(model, ss) %>% 
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$model, loci.labels$ss), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL)

gp <- ggplot(data = res2211, aes(x = X.Time, y = Median, linetype = as.factor(greedy), colour = as.factor(minbin), group = interaction(greedy, minbin))) +
  geom_path() +
  #geom_line(data = cursims, aes(x = X.Time, y = nind, colour = as.factor(pop_init):as.factor(crash_prop))) +
  geom_vline(xintercept = res2211$ts[1], colour = "blue") +
  geom_vline(xintercept = (res2211$ts[1] - res2211$tl[1]), colour = "black") +
  ggtitle(paste0("Initial pop size = ", res2211$pop_init[1] ,", start: " , res2211$ts[1] , "ybp" ,", length: " , res2211$tl[1] , "yrs" )) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(0, 1000) +
  ylim(0, 600) +
  facet_grid(ss~model) 

gp

grid.arrange(gp, tb, ncol = 2)
gp

grid.arrange(tb )
trialres
