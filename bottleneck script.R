#======================BOTTLENECK MODEL SCRIPT================================#

library(dartR)
library(devtools)
library(slimr)
library(future)
library(geohippos)
library(parallel)
library(plyr)
library(doParallel)
library(tidyr)
library(dplyr)
library(readr)

##If slim_run doesn't work
Sys.setenv(SLIM_HOME='/home/isobel/slim/bin/slim/bin')
#Sys.setenv(SLIM_HOME = "/data/scratch/isobel/win-library/4.1/slimr")
slim_setup()

#======================SIMULATION SCRIPT======================================
slim_script(
  slim_block(initialize(),
             {
               
               defineConstant("alpha2", 1e-8);
               defineConstant("L", 5e8);
               defineConstant("model", slimr_template("model"));
               defineConstant("pop_init", slimr_template("pop_init"));
               defineConstant("ts", slimr_template("ts"));
               defineConstant("tl", slimr_template("tl"));
               defineConstant("cp", slimr_template("crash_prop"));
               defineConstant("ss", slimr_template("ss"));
               defineConstant("filename", slimr_template("filename"));
               
               
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
  slim_block(1, early(),
             {
               ##Define initial population size
                 sim.addSubpop("p1", pop_init);

               # 
               # logfilename = paste0("sim_ind_folder/sim_ind_" + slimr_template("model") + pop_init + cp + ts + tl + ss +".txt");
               # log = community.createLogFile(logfilename, logInterval=10);
               # log.addCycle();
               # log.addCustomColumn("nind", "p1.individualCount;");
               # log.addCustomColumn("model", "model;");
               # log.addCustomColumn("pop_init", "pop_init;");
               # log.addCustomColumn("crash_prop", "cp;");
               # log.addCustomColumn("ts", "ts;");
               # log.addCustomColumn("tl", "tl;");
               # log.addCustomColumn("ss", "ss;");
               # 
             }),
  
  slim_block(1000,1200, late(),
             
             {
               startyr = 1200 - ts;
               endyr = 1200 - ts + tl;
               ###Abort simulation if trajectory length is not compatible with trajectory start time
               if(tl > ts) {
                 sim.simulationFinished();
                 ###create error message###
               }
               ##separate by trajectory start - delay until start time
               if ((sim.cycle + 1) > startyr) {

                   yr = sim.cycle - startyr;
                   
                   ###population crash
                   if (sim.cycle < (startyr + 0.1*tl)) {
                     x1 = pop_init;
                     xn = cp * pop_init;
                     tl1 = 0.1 * tl;
                     r1 = log(x1 - xn)/tl1;
                     p1.setSubpopulationSize(asInteger(round(
                       (x1 - xn)*exp((-r1) * yr) + xn)
                     ));
                   }
                   
                   # ###Simulate crash - 10% of trajectory length
                   # if (yr < (0.1 * tl) + 1) {
                   #   #crash length = cl
                   #   cl = 0.1*tl;
                   #   r = 4*log(10)/(cl);
                   #   p1.setSubpopulationSize(asInteger(round(
                   #     ((popin - k)/(1 + ((popin/k) - 1)*exp(r*(yr - 0.5*cl))))+k)
                   #   )); 
                   # }
                   
                   ##Recovery after crash
                   if (sim.cycle > (startyr + (0.8*tl))) {
                     x1 = p1.individualCount;
                     xn = pop_init;
                     tl2 = 0.2 * tl;
                     r2 = log(xn/x1)/tl2;
                     p1.setSubpopulationSize(asInteger(round(
                       x1*exp(r2*yr))
                     ));
                     if (p1.individualCount > pop_init) {
                       p1.setSubpopulationSize(asInteger(pop_init));
                     }
                   }
               }
               
             }),
  
  slim_block(1200,late(),
             {
               nn = filename;
               p1.outputVCFSample(sampleSize=ss, replace=F,  outputMultiallelics=F,filePath=nn,  simplifyNucleotides=T);
               sim.simulationFinished();   
             })
) -> script_bottle

###select folder for vcf files to go###
outdir <- "/data/scratch/isobel/vcf/"

###create dataframe of all test combinations
df <- expand.grid(rep =1,
                  pop_init = c(1000, 500), 
                  crash_prop = c(0.5), 
                  tl = c(100), #trajectory length (from ts)
                  ts = c(200), #trajectory start (ybp)
                  crl = c(10), #crash length - time for population to drop to crashed size
                  rl = c(20), #recovery length, time for population to recover from crashed size
                  lag = c(50), #lag time between pop crash and recovery
                  ss = c(100, 75))  


df$runnumb <- paste0("Run_",sprintf("%05d",1:nrow(df)))

#==============Run all df combinations through SLiM=============================
testscript <- slim_script_render(script_1, template = df)
testscript
`
###run multiple scripts
plan(multisession, workers = 7)
slim_run(testscript[[1]], parallel = T)

