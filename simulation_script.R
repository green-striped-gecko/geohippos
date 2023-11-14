
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
               if (model != 2) {
                 sim.addSubpop("p1", pop_init);
               }
               if (model == 2) {
                 exp_pop_init = cp * pop_init;
                 sim.addSubpop("p1", asInteger(round(exp_pop_init)));
               }
               
               logfilename = paste0("sim_ind_folder/sim_ind_" + slimr_template("model") + pop_init + cp + ts + tl + ss +".txt");
               log = community.createLogFile(logfilename, logInterval=10);
               log.addCycle();
               log.addCustomColumn("nind", "p1.individualCount;");
               log.addCustomColumn("model", "model;");
               log.addCustomColumn("pop_init", "pop_init;");
               log.addCustomColumn("crash_prop", "cp;");
               log.addCustomColumn("ts", "ts;");
               log.addCustomColumn("tl", "tl;");
               log.addCustomColumn("ss", "ss;");
               
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
                
                ######Decline model#######
                if (model == 1) {
                   
                   ###run until end of trajectory length, relative to start time
                     x1 = pop_init;
                     xn = cp * pop_init;
                     yr = sim.cycle - startyr;
                     r = log(x1 - xn)/tl;
                     p1.setSubpopulationSize(asInteger(round(
                       (x1 - xn)*exp((-r) * yr) + xn)
                       )); 
                  
                  }
                
                ######Expansion model######
                if (model == 2) {
                  
                  ###run until end of trajectory length, relative to start time
                  if (sim.cycle < (endyr + 1)) {
                      x1 = cp * pop_init;
                      xn =  pop_init;
                      yr = sim.cycle - startyr;
                      r = log(xn/x1)/tl;
                      
                      p1.setSubpopulationSize(asInteger(round(
                        x1*exp(r*yr))
                      ));
                  }
                    if (sim.cycle > endyr) {
                        p1.setSubpopulationSize(asInteger(pop_init));
                  }
                }
                
                ######Bottleneck model########
                if (model == 3) {
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
                  if (sim.cycle > (startyr + 0.1*tl)) {
                    x1 = p1.individualCount;
                    xn = pop_init;
                    tl2 = 0.9 * tl;
                    r2 = log(xn/x1)/tl2;
                    p1.setSubpopulationSize(asInteger(round(
                      x1*exp(r2*yr))
                    ));
                    if (p1.individualCount > pop_init) {
                      p1.setSubpopulationSize(asInteger(pop_init));
                    }
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
) -> script_1

###select folder for vcf files to go###
outdir <- "c:/temp/"

###create dataframe of all test combinations
df <- expand.grid(rep =1,
                  pop_init = c(500, 200), 
                  crash_prop = c(0.1, 0.5), 
                  tl = c(100),
                  ts = c(200),
                  ss = c(20),
                  model = c("decline", "expansion", "bottle", "stable"))  

df$filename <- paste0(outdir, "slim_rep", df$rep, "_pop", df$pop_init, "_tl", df$tl, "_crash", df$crash_prop, "_ts-ybp", df$ts, "_model", df$model, "_ss", df$ss, "test13-11.vcf")
df$runnumb <- paste0("Run_",sprintf("%05d",1:10))

#==============Filter out dataframes where ss > final pop=================================
#check for any rows 
df %>% filter(model == "decline") %>% filter(as.numeric(ss) > round((pop_init * crash_prop)))
#get index numbers
ind.remove <- which(df$model == "decline" & (df$ss > (df$pop_init * df$crash_prop)))
#remove from df
df <- df[-ind.remove,]
#==============Runs with higher sample size than final pop are removed===============




#==============Run all df combinations through SLiM==================================

testscript <- slim_script_render(script_1, template = df)
testscript

###run multiple scripts
plan(multisession, workers = 7)
slim_run(testscript, parallel = T)

#==============================================================================



#==============add simulated ind filename to dataframe===============================

modelindex <- data.frame(model = c("decline", "expansion", "bottle", "stable"), n = 1:4)

df <- df %>% 
  mutate(simindfile = paste0("C:/Users/Isobel/Desktop/Honours/geohippos/sim_ind_", 
                             which(modelindex$model == model), 
                                pop_init, crash_prop, ts, tl, ss, ".txt"))
df$simindfile

##read in simind data
for (i in 1:nrow(df)) {
  df$siminddata[[i]] <- read.delim2(df$simindfile[[i]], header = T, sep = ",")
}
#====================================================================================


#=====Reassign df=============
dfshorttest <- as_tibble(df)


#===================convert all to gls and run epos==================================
for (i in 1:nrow(df)) {
  df$gls[[i]] <- geohippos::gl.read.vcf(df$filename[i]);
  df$epos[[i]] <- gl.epos(df$gls[[i]], epos.path = paste0("./binaries/epos/",os),
                          l = L, u=mu, boot=50, )
}
#===============================================================================

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


####Run parallel
library(doParallel)
cl <- makeCluster(4)
registerDoParallel()

##############Epos testing################
settings <- expand_grid(minbin = 1:4, greedy = c(" -E 2", " -E 5", " -E 10", " -E 20"))
settings
test.epos <- crossing(df, settings)
test.epos



###convert all to gls and run epos
for (i in 1:nrow(test.epos[1:4,])) {
  # test.epos$gls[[i]] <- gl.read.vcf(test.epos$filename[i]);
  test.epos$eposout[[i]] <- gl.epos(test.epos$gls[[i]], epos.path = paste0("./binaries/epos/",os),
                          l = L, u=mu, boot=50, minbinsize = test.epos$minbin[[i]], other.options = test.epos$greedy[[i]])
}

ggplot(data = test.epos, aes(x))

df$filename[[1]]
gls <- gl.read.vcf(df$filename[[3]])
###############
#     Epos    #
###############
system.time(
  Ne_epos <- gl.epos(gls, epos.path = paste0("./binaries/epos/",os), l = L, u=mu, boot=50, other.options = " -E 10", minbinsize = 1)
  
)
i = 3
data = Ne_epos
filename = runtest$filename

plot(Median ~ (X.Time), data=Ne_epos, type="l", lwd=2, xlim = c(0,500), ylim = c(0, 1000), main = df$filename[[i]])
points(LowerQ ~ (X.Time), data=data, type="l", col="blue", lty=2)
points(UpperQ ~ (X.Time), data=data, type="l", col="orange", lty=2)





plot(Median ~ (X.Time), data=df$epos[[i]], type="l", lwd=2, xlim = c(0,400), ylim = c(0, 200), main = df$filename[[i]])
points(LowerQ ~ (X.Time), data=df$epos[[i]], type="l", col="blue", lty=2)
points(UpperQ ~ (X.Time), data=df$epos[[i]], type="l", col="orange", lty=2)

epos.list <- expand.grid(gls = df$gls, boots = c(25, 50, 75, 100))

for (i in 1:nrow(epos.list)) {
  epos.list$Ne[[i]] <- gl.epos(epos.list$gls[[i]], epos.path = paste0("./binaries/epos/",os), l = L, u=mu, boot=epos.list$boots[[i]])
}

epos.list$Ne

