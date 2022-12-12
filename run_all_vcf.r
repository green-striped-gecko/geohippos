

#library(dartR)
#library(devtools)
#install_github("https://github.com/green-striped-gecko/geohippos")
library(geohippos)

library(tictoc)
tic()
#gls <- geohippos::gl.read.vcf("./inst/extdata/slim_5c_100.vcf", verbose=0)
gls <- geohippos::gl.read.vcf("./vcf/slim_200_20-50y_100-5y.vcf")
#gls <- geohippos::gl.read.vcf("d:/downloads/slim_50_100inlast10.vcf")
#split chromosomes...

gls$chromosome <- factor(ceiling(gls$position/1e8)) #slim simulation
table(gls@chromosome)

#some checks on the input file
nLoc(gls)
nInd(gls)
pop(gls)<- rep("A", nInd(gls))
sfs <-gl.sfs(gls)


#gls <- gls[,gls@chromosome==1 | gls@chromosome==4]
#gls <- gls[,gls@chromosome==4]


#check os to find correct binaries
os <- tolower(Sys.info()['sysname']) 
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate
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
#plot(Median ~ (X.Time), data=Ne_epos, type="l", lwd=2)
#points(LowerQ ~ (X.Time), data=Ne_epos, type="l", col="blue", lty=2)
#points(UpperQ ~ (X.Time), data=Ne_epos, type="l", col="orange", lty=2)

###############
#  STAIRWAYS  #
###############
system.time(
  Ne_sw <- gl.stairway2(gls,simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T)
)
#plot(Ne_sw$year, Ne_sw$Ne_median, type="l", lwd=2, xlab="year", ylab="Ne")
#points(Ne_sw$year, Ne_sw$Ne_12.5., type="l", lty=2, col="blue")
#points(Ne_sw$year, Ne_sw$Ne_87.5., type="l", lty=2, col="blue")
#points(Ne_sw$year, Ne_sw$Ne_2.5., type="l", lty=2, col="red")
#points(Ne_sw$year, Ne_sw$Ne_97.5., type="l", lty=2, col="red")



###############
#     GONE    #
###############
system.time(
  Ne_gone <- gl.gone(gls,gone.path = paste0("./binaries/gone/",os)) #runs parallel via InputParamters
)
#plot(Ne_gone$Generation, Ne_gone$Geometric_mean, type="l")


###############
#   SNEP      #
###############
system.time(
  Ne_snep <- gl.snep(gls, snep.path = paste0("./binaries/snep/",os), n.cores=30)
)
#plot(Ne ~ GenAgo, data=Ne_snep, type="l")



###############
#   LinkNe    #
###############
system.time(
Ne_LinkNe <- gl.LinkNe(gls, outfile = "trun", LinkNe.path = paste0("./binaries/linkne/",os), perl = FALSE)
)
#plot(1/(2*Ne_LinkNe$MEAN_C), Ne_LinkNe$NE, type="l")

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
dummy <-data.frame(nr =1:4)

dummy$year <- c(0,5,   5, 50 , 50, 500)
dummy$Ne <- c(100,100,20, 20, 200 ,200)
##dummy$year <- c(0,25,25,500)
#dummy$Ne <- c(200,200,100,100)
#fox
#dummy$year <- c(0,10,30,30,50,50,500)
#dummy$Ne <- c(600,600,600,200,200,10000,10000)


dummy$method="sim"
ress[[6]]<- dummy
res <- do.call(rbind, ress)
library(ggplot2)


ggplot(res, aes(x=year, y=Ne, color=method))+geom_line()+facet_wrap(. ~  method, scales = "free")
ggplot(res, aes(x=year, y=Ne, color=method, group=method))+geom_line()+xlim(c(2,500))+ylim(c(0,1000))


toc()









