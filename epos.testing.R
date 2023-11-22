#===============================================
#           Test Method Settings
#===============================================
library(dplyr)
library(data.table)
library(tidyr)
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
  test.epos$loci[[i]] <- nLoc(test.epos$gls[[i]])
}

#==============add simulated ind filename to dataframe===============================

modelindex <- data.frame(model = c("decline", "expansion", "bottle", "stable"), n = 1:4)


test.epos <- test.epos %>% 
  mutate(simindfile = paste0("C:/Users/Isobel/Desktop/Honours/geohippos/sim_ind_folder/sim_ind_", 
                             match(model, modelindex$model), 
                             pop_init, crash_prop, ts, tl,".txt"))
test.epos$simindfile

##read in simind data
for (i in 1:nrow(test.epos)) {
  test.epos$siminddata[[i]] <- read.delim2(test.epos$simindfile[[i]], header = T, sep = ",")
}

#==============Extract data into usable dataframe=============================#

extract_epos <- function(df1, df2) {
  if(nrow(df1) != ncol(df2)) {
    print("Error: Dataframe rows not equal to column entries!")
    return()
  };
  res_epos <- data.frame();
  for (i in 1:nrow(df1)) {
    new_data <- cbind(df1[i,c(1:11)], as.data.frame(df2[,i]));
    res_epos <- rbind(res_epos, new_data);
  }
  return(res_epos);
}

df_extract_epos <- function(df1) {
  res_epos <- data.frame();
  for (i in 1:nrow(df1)) {
    new_data <- cbind(df1$eposout1[[i]], df1[i,c(1:7, 9:10)]);
    res_epos <- rbind(res_epos, new_data);
  }
  return(res_epos);
}

extract_sim <- function(df) {
  sim_data <- data.frame()
  for (i in 1:nrow(df)) {
    sim_data <- rbind(sim_data, df$siminddata[[i]]);
  }
  return(sim_data);
}

 

#=================Extract data from parallel epos output=========================
for (i in 1:ncol(epos.output)) {
  test.epos$eposout[[i]] <- as.data.frame(epos.output[,i])
}


siminds <- extract_sim(test.epos)
siminds$X.Time <- 2000 - siminds$cycle
siminds$model <- modelindex$model[siminds$model]

siminds$crash_prop <- as.double(siminds$crash_prop)

res500 <- extract_epos(test.epos, epos.output)
res500

smallres <- extract_epos(small.eposdf, epos.output)



trialres <- df_extract_epos(trial)


left_join(x = res500, y = siminds, by = c("pop_init", "crash_prop", "model", "ts", "tl"), relationship = "many-to-many")



bin = 1 #(1:4)
othop = "" #c("", " -E 2", " -E 5", " -E 10", " -E 20")
p0 <- 200
cp <- 0.1
mod = "decline"

curdata <- trialres %>% filter(model == mod)

cursims <- siminds %>% filter(model == mod)

#extract loci data
loci.labels <- curdata %>% group_by(pop_init, crash_prop) %>% 
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$pop_init, loci.labels$crash_prop), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL)

gp <- ggplot(data = curdata, aes(x = X.Time, y = Median, colour = as.factor(pop_init):as.factor(crash_prop))) +
  geom_line() +
  #geom_line(data = cursims, aes(x = X.Time, y = nind, colour = as.factor(pop_init):as.factor(crash_prop))) +
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
gp


trialres
