####Decline function writing

library(slimr)
library(dartR)
library(tidyr)
gl.set.verbosity(0)

##If slim_run doesn't work
# Sys.setenv(SLIM_HOME = "C:/Users/Isobel/Documents/R/win-library/4.1/slimr")
# slim_setup()

####### SLiM Run Function Writing###########

###define parameters########
genome_len =  5e8
mut_rate = 1e8
pop.init = c(200, 300, 400, 500)
dec.rate = c(0.99, 0.98, 0.97)

simset <- expand.grid(genome_len = genome_len , mut_rate = mut_rate, pop.init = pop.init, dec.rate = dec.rate)

slim_decline <- function(genome_len, mut_rate, pop.init, dec.rate) {
  
  #create combinations of all variables
  simset <- expand.grid(genome_len = genome_len , mut_rate = mut_rate, pop.init = pop.init, dec.rate = dec.rate);
  
  #create scripts for every combination
  decline_set <- slim_script_render(decline_script, template = simset);
  
  #run slim for every script
  for (i in 1:length(decline_set)) {
    slim_run(decline_set[[i]]);
    print(decline_set[[i]]);
  }
  

}


##Define slim script for decline
slim_script(
  slim_block(initialize(),
             {
               
               defineConstant("alpha2", slimr_template("mut_rate", 1e-8, unquote_strings = T));
               defineConstant("L", slimr_template("genome_len", 5e8, unquote_strings = T));
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
               initializeGenomicElement(g1, 0, L - 1);
               
               #initializeRecombinationRate(rates = 1e-8);
               initializeRecombinationRate(rates =c(1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8), 
                                           ends=c(1e8-1,1e8,2e8-1,2e8,3e8-1,3e8,4e8-1,4e8,5e8-1));     
             }),
  slim_block(1, early(),
             {
               #defineConstant("pop.init", slimr_inline(pop.init));
               sim.addSubpop("p1", slimr_template("pop.init", 200, unquote_strings = T));
             }),
  slim_block(2000,2050, early(),
             {
               newsize <- asInteger(p1.individualCount * slimr_template("dec.rate", 0.99, unquote_strings = T));
               p1.setSubpopulationSize(newsize);
             }),
  
  #  slim_block(2050,2100,late(),
  #             {
  #             p1.setSubpopulationSize(100+(sim.cycle-2050)*4);
  #             #  p1.setSubpopulationSize(200);
  #               
  #             }),
  
  slim_block(2050,late(),
             {
               nn = paste0("c:/temp/slim_", slimr_template("pop.init", 200), "_", slimr_template("dec.rate", 0.99), "_multitest_decline.vcf");
               slimr_output(p1.outputVCFSample(sampleSize=40, replace=F,  outputMultiallelics=F, filePath=nn,  simplifyNucleotides=T), name = "nn");
               sim.simulationFinished();
             })
) -> decline_script
