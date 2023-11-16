#===============================================
#           Test Method Settings
#===============================================
library(dplyr)

#Read in dataframe 

#eposdf <- df


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


####Run parallel
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel()

##############Epos testing################
# settings <- expand_grid(minbin = 1:4, greedy = c("", " -E 2", " -E 5", " -E 10", " -E 20"), .name_repair = "minimal")
# settings
# test.epos <- crossing(eposdf, settings)
# test.epos
# 
library(geohippos)
setwd("C:/Users/Isobel/Desktop/Honours/geohippos")
###convert all to gls and run epos
test.epos.all <- test.epos
test.epos <- test.epos.all
test.epos$eposout <- as.data.frame()
for (i in 1:nrow(test.epos)) {
  #test.epos$gls[[i]] <- gl.read.vcf(test.epos$filename[[i]]);
  test.epos$eposout1[[i]] <- gl.epos(test.epos$gls[[i]], epos.path = paste0("./binaries/epos/",os),l = L, u=mu, boot=50, minbinsize = test.epos$minbin[[i]])
}

#================Extract loci for each run====================================

for (i in 1:nrow(test.epos)) {
  test.epos$loci[[i]] <- nLoc(test.epos$gls[[i]])
}




for (i in 1:nrow(settings)) {
  currset <- settings[i,];
  curdf <- test.epos %>% filter(minbin == currset$minbin, greedy == currset$greedy);
  #ggplot(data = curdf, )
  print(head(curdf));
}

test.epos %>% filter(minbin == 1)
 
trial <- test.epos[1:5,]
trial

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

extract_epos <- function(df) {
  res_epos <- data.frame()
  for (i in 1:nrow(df)) {
    new_data <- cbind( df[i,c(1:10, 14)], df$eposout1[[i]]);
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

siminds <- extract_sim(test.epos)
siminds$X.Time <- 2000 - siminds$cycle
siminds$model <- modelindex$model[siminds$model]

siminds$crash_prop <- as.double(siminds$crash_prop)

ress <- extract_epos(test.epos)
ress

left_join(x = ress, y = siminds, by = c("pop_init", "crash_prop", "model", "ts", "tl"), relationship = "many-to-many")



bin = 1 #(1:4)
othop = "" #c("", " -E 2", " -E 5", " -E 10", " -E 20")
p0 <- 200
cp <- 0.1
mod = "stable"

curdata <- ress %>% filter(model == mod)

cursims <- siminds %>% filter(model == mod)

#extract loci data
loci.labels <- curdata %>% group_by(pop_init, crash_prop) %>% 
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$pop_init, loci.labels$crash_prop), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL)

gp <- ggplot(data = curdata, aes(x = X.Time, y = Median, colour = as.factor(pop_init):as.factor(crash_prop))) +
  geom_line() +
  geom_line(data = cursims, aes(x = X.Time, y = nind, colour = as.factor(pop_init):as.factor(crash_prop))) +
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
