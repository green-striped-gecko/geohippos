####Run inference methods####


library(dartR)
library(devtools)
#install_git("https://github.com/green-striped-gecko/geohippos")
library(geohippos)

library(tictoc)
tic()
#gls <- geohippos::gl.read.vcf("./inst/extdata/slim_5c_100.vcf", verbose=0)
gls <- geohippos::gl.read.vcf("C:/temp/slim_300_0.98_multitest_decline.vcf")
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
genome_len =  5e8
mut_rate = 1e-8
pop.init = c(200, 300, 400, 500)
dec.rate = c(0.99, 0.98, 0.97)


run_epos <- function(simset, nboots) {
  ##create empty dataframe for results
  epos_results <- as.list(simset);
  i = 1
  for (i in 1:length(simset)) {
    
  
    L <- genome_len; #total length of chromosome (for sfs methods)
    mu <- mut_rate;  #mutation rate
    gls <- geohippos::gl.read.vcf(simset[i,]$filename);
    #split chromosomes...
    
    gls$chromosome <- factor(ceiling(gls$position/1e8)); #slim simulation
    table(gls@chromosome);
    #some checks on the input file
    if (is.na(nLoc(gls))) {
      print("Warning: nLoc = NA")
    };
    if (is.na(nInd(gls))) {
      print("Warning: nInd = NA")
    };
    pop(gls)<- rep("A", nInd(gls))
    #sfs <-gl.sfs(gls);
    
    #epos_results[i,]$sfs <- sfs;
    epos_results[i,]$nboots <- nboots;
    
  system.time(
    Ne_epos <- gl.epos(gls, epos.path = paste0("./binaries/epos/",os), l = L, u=mu, boot= nboots))
    if (i == 1){
      epos_results <- Ne_epos;
      epos_results
  };
  return(epos_results)
}

test <- run_epos(simset, 10)

plot(Median ~ (X.Time), data=Ne_epos, type="l", lwd=2, xlim = c(0,100))
points(LowerQ ~ (X.Time), data=Ne_epos, type="l", col="blue", lty=2)
points(UpperQ ~ (X.Time), data=Ne_epos, type="l", col="orange", lty=2)

###############
#  STAIRWAYS  #
###############
system.time(
  Ne_sw <- gl.stairway2(gls,simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T)
)
plot(Ne_sw$year, Ne_sw$Ne_median, type="l", lwd=2, xlab="year", ylab="Ne")
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
plot(Ne_gone$Generation, Ne_gone$Geometric_mean, type="l")


###############
#   SNEP      #
###############
system.time(
  Ne_snep <- gl.snep(gls, snep.path = paste0("./binaries/snep/",os), n.cores=3)
)
plot(Ne ~ GenAgo, data=Ne_snep, type="l")



###############
#   LinkNe    #
###############
system.time(
  Ne_LinkNe <- gl.LinkNe(gls, outfile = "trun", LinkNe.path = paste0("./binaries/linkne/",os), perl = FALSE)
)
plot(1/(2*Ne_LinkNe$MEAN_C), Ne_LinkNe$NE, type="l")

ress <- list()


