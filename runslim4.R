
#better by isobel
install.packages("reticulate")

library(slimr)
library(dartR)
gl.set.verbosity(0)


slim_script(
  slim_block(initialize(),
             {
               defineConstant("alpha2", 1e-8);
               defineConstant("L", 5e8);
               #defineConstant("simlen", 5000);
               #how many outputs (simlen=only once at the end)
               initializeSLiMOptions(nucleotideBased=T);
               initializeAncestralNucleotides(randomNucleotides(L));
               initializeMutationTypeNuc("m1", 0.5, "f",0.0 );
               m1.mutationStackPolicy = "l";
               m1.convertToSubstitution=T;
               defineConstant("mm2", matrix(c(0.0,alpha2,0.0,0.0, alpha2,0.0,0.0,0.0, 0.0,0.0,0.0,alpha2, 0.0,0.0,alpha2,0.0 ),nrow=4, ncol=4));
               #defineConstant("mm2", mmJukesCantor((1e-8)/3));
               
               initializeGenomicElementType("g1", m1, 1.0, mm2);
               initializeGenomicElement(g1, 0, L-1);
               
               #initializeRecombinationRate(rates = 1e-8);
               initializeRecombinationRate(rates = c(1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8), ends=c(1e8-1,1e8,2e8-1,2e8,3e8-1,3e8,4e8-1,4e8,5e8-1));     
             }),
  slim_block(1, early(),
             {
               sim.addSubpop("p1", 500);
             }),
  slim_block(2000, early(),
             {
               newsize <- asInteger(p1.individualCount * 0.99);
               p1.setSubpopulationSize(newsize);
             }),

  
  slim_block(2050,late(),
             {
               nn = "c:/temp/slim_500_slowdecline2.vcf";
               p1.outputVCFSample(sampleSize=20, replace=F,  outputMultiallelics=F,filePath=nn,  simplifyNucleotides=T);
               sim.simulationFinished();
             })
) -> script_1

script_1

#run slim from within slimr

slim_run(script_1)


library(geohippos)

#load back into R 

gls <- geohippos::gl.read.vcf("c:/temp/slim_500_slowdecline2.vcf", verbose=0)

#to split into chromosomes...
gls$chromosome <- factor(ceiling(gls$position/1e8)) #slim simulation
table(gls@chromosome)
nLoc(gls)
nInd(gls)
pop(gls)<- rep("A", nInd(gls))
sfs <-gl.sfs(gls)


os <- tolower(Sys.info()['sysname']) 
if (os=="darwin") os <- "mac"

L <- 10e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate

###############
#     Epos    #
###############
system.time(
  Ne_epos <- gl.epos(gls, epos.path = paste0("./binaries/epos/",os), l = L, u=mu, boot=5, minbinsize =1)
  
)
plot(Median ~ (X.Time), data=Ne_epos, type="l", lwd=2, ylim = c(0, 400), xlim = c(0, 100))
points(LowerQ ~ (X.Time), data=Ne_epos, type="l", col="blue", lty=2)
points(UpperQ ~ (X.Time), data=Ne_epos, type="l", col="orange", lty=2)

###############
#  STAIRWAYS  #
###############
system.time(
  Ne_sw <- gl.stairway2(gls, simfolder = "stairwaytest", verbose = T,stairway.path="C:/Users/isobe/Desktop/Honours/geohippos/binaries/stairways", mu = mu, gentime = 1, run=TRUE, nreps = 2, parallel=1, L=L, minbinsize =1, cleanup = T))

plot(Ne_sw$year, Ne_sw$Ne_median, type="l", lwd=2, xlab="year", ylab="Ne")
points(Ne_sw$year, Ne_sw$Ne_12.5., type="l", lty=2, col="blue")
points(Ne_sw$year, Ne_sw$Ne_87.5., type="l", lty=2, col="blue")
points(Ne_sw$year, Ne_sw$Ne_2.5., type="l", lty=2, col="green")
points(Ne_sw$year, Ne_sw$Ne_97.5., type="l", lty=2, col="red")

###############
# Neestimator #
###############
system.time(
  Ne_ldnest <- gl.LDNe(gls,neest.path = paste0("./binaries/NEestimator/",os), singleton.rm = F, critical = c(0))
)

gl.gone(gls, outfile = "test", gone.path = "C:/Users/isobe/Desktop/Honours/geohippos/binaries/gone")



