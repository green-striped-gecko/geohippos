#############run slim gui outputs ########

guidf <- expand.grid(rep =1:5,
                  pop_init = c(500, 200, 100), 
                  crash_prop = c(0.5, 0.1), 
                  tl = c(100, 50, 25),
                  ts = c(200, 100, 50),
                  ss = c(15, 20, 30),
                  model = c("decline", "expansion", "bottle", "stable"))



library(dartR)
library(devtools)
#install_git("https://github.com/green-striped-gecko/geohippos")
library(geohippos)

library(tictoc)

###output directory##
outdir <- "c:/temp/"

guidf <- expand.grid(rep =1,
                     pop_init = c(200), 
                     crash_prop = c(0.1), 
                     tl = c(100),
                     ts = c(100),
                     ss = c(15),
                     model = c("decline", "expan", "bottle", "stable"))


guidf <- guidf %>% 
  mutate(filename = paste0(outdir, model, "_p0_", pop_init , "_cp_" , crash_prop , "_ts_" , ts , "_tl_" , tl , ".vcf"))
guidf

#####create gls file for each vcf####
guidf$gls <- lapply(guidf$filename, function(x) geohippos::gl.read.vcf(x))
guidf$sls <- lapply(guidf$gls, function(x) gl.sfs(x))

guidf$sls[1]

ggplot(data = guidf) +
  geom_col(aes(x = sfs, stat = "identity"))

results <- tibble(guidf)
results

tic()

gls <- geohippos::gl.read.vcf(guidf$filename[3])

gls$chromosome <- factor(ceiling(gls$position/1e8)) #slim simulation
table(gls@chromosome)


#some checks on the input file
nLoc(gls)
nInd(gls)
pop(gls)<- rep("A", nInd(gls))
sfs <-gl.sfs(gls)


#check os to find correct binaries
os <- tolower(Sys.info()['sysname']) 
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate


> setwd("C:/Users/Isobel/Desktop/Honours/geohippos")

guidf
for (i in 1:nrow(guidf)) {
  guidf$eposgreedy[[i]] <- gl.epos(guidf$gls[[i]], epos.path = paste0("./binaries/epos/",os),
                        l = L, u=mu, boot=50);
  guidf$epos10[[i]] <- gl.epos(guidf$gls[[i]], epos.path = paste0("./binaries/epos/",os),
                               l = L, u=mu, boot=50, other.options = " -E 10");
  guidf$epos20[[i]] <- gl.epos(guidf$gls[[i]], epos.path = paste0("./binaries/epos/",os),
                               l = L, u=mu, boot=50, other.options = " -E 20");
  guidf$epos30[[i]] <- gl.epos(guidf$gls[[i]], epos.path = paste0("./binaries/epos/",os),
                               l = L, u=mu, boot=50, other.options = " -E 30");
  guidf$swrand0[[i]] <- gl.stairway2(guidf$gls[[i]], simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nrand = NULL, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T);
  guidf$swrand1[[i]] <- gl.stairway2(guidf$gls[[i]], simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nrand = 1, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T);
  guidf$swrand2[[i]] <- gl.stairway2(guidf$gls[[i]], simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nrand = 2, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T);
}



###############
# Neestimator #
###############
system.time(
  Ne_ldnest <- gl.LDNe(gls,neest.path = paste0("./binaries/NEestimator/",os), singleton.rm = F, critical = c(0))
)

###############
#     Epos    #
###############
system.time(
  Ne_epos <- gl.epos(gls, epos.path = paste0("./binaries/epos/",os), l = L, u=mu, boot=50)
  
)


plot(Median ~ (X.Time), data=Ne_epos, type="l", lwd=2, xlim = c(0,400))
points(LowerQ ~ (X.Time), data=Ne_epos, type="l", col="blue", lty=2)
points(UpperQ ~ (X.Time), data=Ne_epos, type="l", col="orange", lty=2)

###############
#  STAIRWAYS  #
###############

##breakpoint variations
#nrand = 1: 5 breakpoints, nrand = 2: 8 breakpoints


gls <- guidf$gls[[1]]


system.time(
  Ne_sw <- gl.stairway2(gls, simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nrand = 2, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T)
)
plot(Ne_sw$year, Ne_sw$Ne_median, type="l", lwd=2, xlab="year", ylab="Ne")
points(Ne_sw$year, Ne_sw$Ne_12.5., type="l", lty=2, col="blue")
points(Ne_sw$year, Ne_sw$Ne_87.5., type="l", lty=2, col="blue")
points(Ne_sw$year, Ne_sw$Ne_2.5., type="l", lty=2, col="red")
points(Ne_sw$year, Ne_sw$Ne_97.5., type="l", lty=2, col="red")



###############
#     GONE    #
###############
system.time(
  Ne_gone <- geohippos::gl.gone(gls,gone.path = paste0("./binaries/gone/windows/",os)) #runs parallel via InputParamters
)



gone.path = paste0("./binaries/gone/",os)
file.exists(file.path(gone.path, progs))

gl.gone(outfile = )
plot(Ne_gone$Generation, Ne_gone$Geometric_mean, type="l")



ress <- list()


dummy <-data.frame(nr =1:nrow(Ne_epos))
dummy$year <- Ne_epos$X.Time
dummy$Ne<- Ne_epos$Median
dummy$method <- "epos"
ress[[1]]<- dummy

dummy <-data.frame(nr =1:nrow(Ne_sw))
dummy$year <- Ne_sw$year
dummy$Ne <- Ne_sw$Ne_median
dummy$method <- "sw"
ress[[2]]<- dummy


dummy <-data.frame(nr =1:nrow(Ne_gone))
dummy$year <- Ne_gone$Generation
dummy$Ne <- Ne_gone$Geometric_mean
dummy$method <- "gone"
ress[[3]]<- dummy

dummy <-data.frame(nr =1:nrow(Ne_snep))
dummy$year <- Ne_snep$GenAgo
dummy$Ne <- Ne_snep$Ne
dummy$method <- "snep"
ress[[4]]<- dummy


dummy <-data.frame(nr =1:nrow(Ne_LinkNe))
dummy$year <- 1/(2*Ne_LinkNe$MEAN_C)
dummy$Ne <- Ne_LinkNe$NE
dummy$method <- "LinkNe"
ress[[5]]<- dummy

##sim 
dummy <-data.frame(nr =1:3)

dummy$year <- c(0 , 50, 550)
dummy$Ne <- c(302, 500 ,500)
##dummy$year <- c(0,25,25,500)
#dummy$Ne <- c(200,200,100,100)
#fox
#dummy$year <- c(0,10,30,30,50,50,500)
#dummy$Ne <- c(600,600,600,200,200,10000,10000)


dummy$method="sim"
ress[[5]]<- dummy
res <- do.call(rbind, ress)
library(ggplot2)


ggplot(res, aes(x=year, y=Ne, color=method))+geom_line()+facet_wrap(. ~  method, scales = "free")
ggplot(res, aes(x=year, y=Ne, color=method, group=method))+geom_line()+xlim(c(2,500))+ylim(c(0,1000))


toc()



i <- 0
pop <- 500
for (i in 1:50) {
  pop <- pop*0.99
  print(pop)
}

