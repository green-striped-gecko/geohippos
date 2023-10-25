library(slimr)
library(future)
library(dartR.base)
gl.set.verbosity(0)


df <- data.frame(ss=seq(10,50,10))


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
               initializeRecombinationRate(rates = c(1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8), ends=c(1e8-1,1e8,2e8-1,2e8,3e8-1, 3e8,4e8-1, 4e8,5e8-1));     
             }),
  slim_block(1,
             {
               sim.addSubpop("p1", 100);
             }),

  
  slim_block(1000,late(),
             {
             #p1.setSubpopulationSize(200+(sim.cycle-2100));
             p1.setSubpopulationSize(10);
             }),

  slim_block(1050,1100,late(),
             {
               r = log(10)/50
               p1.setSubpopulationSize(asInteger(10*exp(r*(sim.cycle-1050))));
               #p1.setSubpopulationSize(200+(sim.cycle-2100));
               p1.setSubpopulationSize(100);
             }),
  
  
  
    
  slim_block(1100,late(),
             {
               nn = paste("d:/downloads/slim_100-10-100_bottle",slimr_template("ss"),".vcf");
               p1.outputVCFSample(sampleSize=slimr_template("ss"), replace=F,  outputMultiallelics=F,filePath=nn,  simplifyNucleotides=T);
               sim.simulationFinished();
             })
) -> script_1

scripts <- slim_script_render(script_1, df)

#run slim from within slimr
plan(multisession, workers=5)
all <- slim_run(scripts,parallel = T)


library(geohippos)
#load back into R 

gls <- geohippos::gl.read.vcf("d:/downloads/slim_100-10-100_bottle 50 .vcf", verbose=0)

#to split into chromosomes...
gls$chromosome <- factor(ceiling(gls$position/1e8)) #slim simulation
table(gls@chromosome)
nLoc(gls)
nInd(gls)
pop(gls)<- rep("A", nInd(gls))
sfs <-gl.sfs(gls)


os <- tolower(Sys.info()['sysname']) 
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate

###############
#     Epos    #
###############
system.time(
  Ne_epos <- gl.epos(gls, epos.path = paste0("./binaries/epos/",os), l = L, u=mu, boot=10, minbinsize =2)
  
)
plot(Median ~ (X.Time), data=Ne_epos, type="l", lwd=2)
points(LowerQ ~ (X.Time), data=Ne_epos, type="l", col="blue", lty=2)
points(UpperQ ~ (X.Time), data=Ne_epos, type="l", col="orange", lty=2)

###############
#  STAIRWAYS  #
###############
#system.time(
  Ne_sw <- gl.stairway2(gls,simfolder = "stairwaytest", verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize =1, cleanup = T)
#)
plot(Ne_sw$year, Ne_sw$Ne_median, type="l", lwd=2, xlab="year", ylab="Ne")
points(Ne_sw$year, Ne_sw$Ne_12.5., type="l", lty=2, col="blue")
points(Ne_sw$year, Ne_sw$Ne_87.5., type="l", lty=2, col="blue")
points(Ne_sw$year, Ne_sw$Ne_2.5., type="l", lty=2, col="red")
points(Ne_sw$year, Ne_sw$Ne_97.5., type="l", lty=2, col="red")





