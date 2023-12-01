#################################################################
#                     RUN ALL - TESTING                         #
#################################################################


####load necessary packages
library(dplyr)
library(tidyverse)
library(geohippos)
library(dartR)
library(parallel)
library(tictoc)

###prep data for all methods
tr_03 <- read.csv("ss100df.csv")
tr_03

runcode <- "tr_04"

df  <- tr_02
df$filename <- str_replace(df$filename, "c:/temp", "~/R/geohippos/vcf.testing")

##generate gls files for all
setwd("~/R/geohippos")
df$gls <- mclapply(1:nrow(df), function(x) {
  out <- gl.read.vcf(df$filename[[x]])
  return(out)
}, mc.cores=20)

df <- filteringdf

df.f <- df %>% filter(filtering == "prefiltered")
df.uf <- df %>% filter(filtering == "unfiltered")
##generate loci numbers for all
for (i in 1:nrow(df)) {
  df$loci[[i]] <- nLoc(df$gls[[i]])
}

df$loci <- unlist(df$loci)

sfs.list <- list()
for (i in 1:nrow(df.f)) {
  sfs.list[i] <- gl.sfs(df.f$gls[[i]])
}

w#Separate gls and loci columns
dfgls <- df %>% select(runnumb, gls, loci)


####specify settings for each method:
df.cross <- df[1:9]



#check os to find correct binaries
os <- tolower(Sys.info()['sysname'])
if (os=="darwin") os <- "mac"

L <- 5e8 #total length of chromosome (for sfs methods)
mu <- 1e-8  #mutation rate



#===================================#
#           Epos testing            #
#===================================#

#Define settings to test
ep.settings <- expand_grid(minbin = 1, greedy = " -E 2", .name_repair = "minimal")

#Generate all combinations of settings
test.epos <- as_tibble(crossing(df.cross, ep.settings))

test.epos <- df.f




#Read gls and loci information back in
test.epos <- left_join(test.epos, dfgls, by = "runnumb")

test.epos$id <- paste0("eposid",sprintf("%05d",1:nrow(test.epos)))

#Run epos (select for background job)
library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")
test.epos$eposout  <- mclapply(1:nrow(test.epos), function(x) {
  out <- gl.epos(test.epos$gls[[x]], epos.path = paste0("./binaries/epos/",os),l = L, u=mu, boot=20, other.options = " -E 2", minbinsize = 1)
  saveRDS(out, file =  file.path("/data/scratch/isobel/results/", paste0("epos_filttest2_", runcode, test.epos$id[x],".rds")))
  return(out)
}, mc.cores = 10)





#===================================#
#      StairwayPlot testing         #
#===================================#

#Define settings to test
sw.settings <- expand_grid(minbin = c(1:2), breakpoints = c(4, 6), .name_repair = "minimal")


#Generate all combinations of settings
test.sw <- as_tibble(crossing(df.cross, sw.settings))

test.sw <- df.f
#Read gls and loci information back in
test.sw <- left_join(test.sw, dfgls, by = "runnumb")

test.sw$id <- paste0("swid",sprintf("%05d",1:nrow(test.sw)))


library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")

test.sw$stairway <- mclapply(1:nrow(test.sw), function(x) {
  out <- gl.stairway2(test.sw$gls[[x]], verbose = T,stairway.path="./binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 30, parallel=5, L=L, minbinsize = 1, cleanup = T, nrand = 5)
  saveRDS(out, file =  file.path("/data/scratch/isobel/results/", paste0("Stairway_filttest2", runcode, test.sw$id[x],".rds")))
  return(out)
}, mc.cores = 10)



#===================================#
#           GONE testing            #
#===================================#

test.g <- df

library(geohippos)
library(dartR)
library(parallel)
setwd("~/R/geohippos")

####CURRENT SETTINGS#####
##THREADS = 5
##MAF = 0.0

test.g$GONE_10 <- mclapply(1:nrow(test.g), function(x) {
  out <- gl.gone(test.g$gls[[x]],gone.path = paste0("./binaries/gone/",os))
  saveRDS(out, file = file.path('/data/scratch/isobel/results/', paste0("GONE_10", runcode, test.g$runnumb[x], ".rds")))
  return(out)
}, mc.cores = 10)


#===============Stairway plotting prep=========================================#
tr_03_SW_out <- df_extract_output(test.sw, 16, c(2:12, 14:15))
tr_03_SW_out$loci <- unlist(tr_03_SW_out$loci)


#===============GONE plotting prep======================================
#gather outputs
test.g <- test.g %>% gather(key = "maf", value = "output", -c(1:12))

which(test.g$maf == "gone_00")
test.g$maf[test.g$maf == "gone10"] <- 0.10

test.g$loci <- unlist(test.g$loci)

tr_03_gone_out <- df_extract_output(test.g, 14, c(2:9, 12:13))
#unlist loci
tr_03_gone_out$loci <- unlist(tr_03_gone_out$loci)
#alter generation to X.time variable
tr_03_gone_out$X.Time <- 1200 - tr_03_gone_out$Generation

#extract relative generation
max_gens <- tr_03_gone_out %>% group_by(pop_init, model, ss, maf) %>% 
  summarize(gen_max = max(Generation))
tr_03_gone_out <- left_join(tr_03_gone_out, max_gens, by = c('pop_init', "model", "ss", "maf"))
tr_03_gone_out <- tr_03_gone_out %>% mutate(gen_rev = gen_max - Generation)

#===============PLOT OUTPUTS===========================#

####EPOS#######
tr_03_epos <- df_extract_output(test.epos, 16, c(2:12, 14:15))
tr_03_epos$loci <- unlist(tr_03_epos$loci)
tr_03_epos$greedy[tr_03_epos$greedy == ""] <- 0
tr_03_epos$greedy[tr_03_epos$greedy == " -E 2"] <- 2

curdata <- tr_03_epos %>% filter(model == "decline")
#extract loci data
loci.labels <- curdata %>% group_by(pop_init, ss, model) %>%
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$pop_init, loci.labels$model, loci.labels$ss), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL, theme=ttheme_minimal(base_size = 10))

decline <- ggplot(data = curdata, aes(x = X.Time, y = Median, colour = as.factor(pop_init), linetype = as.factor(ss))) +
  geom_line() +
  geom_vline(xintercept = curdata$ts[1], colour = "blue") +
  geom_vline(xintercept = (curdata$ts[1] - curdata$tl[1]), colour = "black") +
  ggtitle(paste0("Pop_init: 1000, 500, Sample size: 100, 75, ",
                 ", start: " , curdata$ts[1] , "ybp" ,
                 ", length: " , curdata$tl[1] , "yrs" )) +
  theme_minimal() +
  xlim(0, 2000) +
  ylim(0, 1200) +
  facet_grid(greedy~minbin, labeller = label_both)

decline$labels$title <- "Decline model"
expan$labels$title <- "Expansion model"
bottle$labels$title <- "Bottleneck model"
stable$labels$title <- "stable model"

grid.arrange(gp,tb, ncol = 2, widths = c(20, 4))

grid.arrange(decline, bottle, expan, stable, top = "Crash: 0.5, tl: 100, ts: 200")


####GONE#####
curdata <- tr_03_gone_out
loci.labels <- curdata %>% group_by(pop_init, model, ss) %>%
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$pop_init, loci.labels$model, loci.labels$ss), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL, theme=ttheme_minimal(base_size = 10))

gp <- ggplot(data = curdata, aes(x = gen_rev, y = Geometric_mean, colour = as.factor(pop_init), linetype = as.factor(ss))) +
  geom_line() +
  geom_vline(xintercept = curdata$ts[1], colour = "blue") +
  geom_vline(xintercept = (curdata$ts[1] - curdata$tl[1]), colour = "black") +
  ggtitle(paste0("Pop_init = 1000, 500, sample_size = 100, 75, ",
                 ", start: " , curdata$ts[1] , "ybp" ,
                 ", length: " , curdata$tl[1] , "yrs" )) +
  theme_minimal() +
  xlim(0, 1200) +
  ylim(0, 10000) +
  facet_grid(model~maf)

grid.arrange(gp,tb, ncol = 2, widths = c(20, 4))
table(curdata$loci, curdata$pop_init + curdata$crash_prop)

####STAIRWAY PLOT####
curdata <- tr_03_SW_out %>% filter(model == "stable")


loci.labels <- curdata %>% group_by(pop_init, ss) %>%
  summarise(loc = mean(loci))
loc.table <- data.frame(key = paste(loci.labels$pop_init, loci.labels$ss), loci = loci.labels$loc)
tb <- tableGrob(loc.table, rows = NULL, theme=ttheme_minimal(base_size = 10))

gp <- ggplot(data = curdata, aes(x = year, y = Ne_median, colour = as.factor(breakpoints), linetype = as.factor(minbin))) +
  geom_line() +
  geom_vline(xintercept = curdata$ts[1], colour = "blue") +
  geom_vline(xintercept = (curdata$ts[1] - curdata$tl[1]), colour = "black") +
  ggtitle(paste0("Model: ", curdata$model[1], " Pop_init: 1000, 500, sample size: 100, 75, ",
                 ", start: " , curdata$ts[1] , "ybp" ,
                 ", length: " , curdata$tl[1] , "yrs" )) +
  theme_minimal() +
  xlim(0, 1200) +
  ylim(0, 2000) +
  facet_grid(pop_init+ss ~ minbin)

grid.arrange(gp,tb, ncol = 2, widths = c(20, 4))

