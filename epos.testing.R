#===============================================
#           Test Method Settings
#===============================================

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

test.epos$eposout1
i = 1
data <- test.epos$eposout1[[i]]
name <- test.epos$filename[[i]]
data <- data %>% gather(key = "stat", value = "val", -X.Time)

ggplot(data, aes(x = X.Time, y = val, colour = stat)) +
  geom_line() +
  ggtitle(name)


for (i in 1:nrow(settings)) {
  
  print(settings[i,]);
}

test.epos %>% filter(minbin == 1)
 
