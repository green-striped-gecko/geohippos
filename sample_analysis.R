library(geohippos)
library(data.table)
library(segmented)
install.packages("cmna", lib = "/data/scratch/isobel")
library(cmna, lib.loc = "/data/scratch/isobel")
###Extract standard dataset for analysis####
fn <- list.files("/data/scratch/isobel/results_gadi/", pattern = "*.rds", full.names = TRUE)
fn

##generate random sample to work on
fn_ss <- sample(fn, 100)

## Load data for sample ##
dat_list = sapply(fn_ss, function (x) data.table(readRDS(x)))

## Convert list names to run numbers ##
names(dat_list) <- unlist(lapply(names(dat_list), function(x) substr(x, 36, 48)))

read.delim2(test.epos$simindfile[[i]], header = T, sep = ",")

## Load data for sample ##
dat_list = sapply(fn_ss, function (x) data.table(readRDS(x)))

## Convert list names to run numbers ##
names(dat_list) <- unlist(lapply(names(dat_list), function(x) substr(x, 36, 48)))
## Generate list of potential files from folder
sim_list <- list.files("/data/scratch/isobel/FOLDER/", pattern = "*.txt", full.names = TRUE)

##Generate list of simulated data
out <- list()
sim_files <- lapply(names(dat_list), function(x) {
  names(x) <- x
  x <- read.delim2(grep(paste0("Run2_", substr(x, 5, 9)), sim_list, value = T), header = T, sep = ",")})
sim_files
names(sim_files) <- names(dat_list)
dat_list

for (i in 1:length(dat_list)) {
  dat_list[[i]][[4]] <- sim_files
}

sim_files_list <- lapply(sim_files, function(X) {
  as.list(X)
})
####Subset data list for GONE outputs#########################
x = dat_list[[5]][[3]]$Generation
y = dat_list[[5]][[3]]$Geometric_mean


####Subset data list for Epos #############################
x = dat_list[[1]][[1]]$X.Time
y = dat_list[[1]][[1]]$Median

####Subset data list for Stairway Plot
x = dat_list[[4]][[2]]$year
y = dat_list[[4]][[2]]$Ne_median


###add simulated data to list
mapply(c, dat_list, sim_files, SIMPLIFY = FALSE)

glm(y ~ x)

model <- segmented(glm(y ~ x))
plot(y ~ x) + plot.segmented(model, conf.level = 0.95, shade = TRUE, col.shade = "blue", add = T)


 lines.segmented(model, )

out <- list()
seg.plots <- lapply(dat_list[1:10], function(X) {
  epos.x <- X[[1]]$X.Time
  epos.y <- X[[1]]$Median
  epos.glm <- glm(epos.y~epos.x)
  epos.model <- segmented(glm(X[[1]]$Median ~ X[[1]]$X.Time))
  
  stair.x <- X[[2]]$year
  stair.y <- X[[2]]$Ne_median
  gone.x <- X[[3]]$Generation
  gone.y <- X[[3]]$Geometric_mean
  out$epos <- data.frame()
  
})


###testing automated glm inputs
ex <- unlist(dat_list$Run_05860_020[[1]][1])
ey <- unlist(dat_list$Run_05860_020[[1]][3])
glm(ey~ex)
test <- lapply(dat_list, function(X) {
  #Generate glms and segms for each method
  e.x <- unlist(X[[1]][1])
  e.y <- unlist(X[[1]][3])
  epos.glm <- glm(e.y ~ e.x)
  epos.seg <- segmented(epos.glm)
  
  s.x <- unlist(X[[2]][6])
  s.y <- unlist(X[[2]][7])
  stair.glm <- glm(s.y ~ s.x)
  stair.seg <- segmented(stair.glm)
  
  g.x <- unlist(X[[3]][1])
  g.y <- unlist(X[[3]][2])
  gone.glm <- glm(g.y ~ g.x)
  gone.seg <- segmented(gone.glm)
})

run_seg <- function(X) {
  #Generate glms and segms for each method
  e.x <- unlist(X[[1]][1])
  e.y <- unlist(X[[1]][3])
  if(length(e.y) < 50){
    print("Warning: Epos output smaller than expected")
  }
  epos.glm <- glm(e.y ~ e.x)
  epos.seg <- segmented(epos.glm)
  epos.plot <- plot(e.y ~ e.x) + plot.segmented(epos.seg, conf.level = 0.95, shade = TRUE, col.shade = "blue", add = T)
  
  s.x <- unlist(X[[2]][6])
  s.y <- unlist(X[[2]][7])
  stair.glm <- glm(s.y ~ s.x)
  stair.seg <- segmented(stair.glm)
  stair.plot <- plot(s.y ~ s.x) + plot.segmented(stair.seg, conf.level = 0.95, shade = TRUE, col.shade = "blue", add = T)
  
  g.x <- unlist(X[[3]][1])
  g.y <- unlist(X[[3]][2])
  gone.glm <- glm(g.y ~ g.x)
  gone.seg <- segmented(gone.glm)
  gone.plot <- plot(g.y ~ g.x) + plot.segmented(gone.seg, conf.level = 0.95, shade = TRUE, col.shade = "blue", add = T)
  
  models <- list(epos.glm, epos.seg, stair.glm, stair.seg, gone.glm, gone.seg)
  plots <- list(epos.plot, stair.plot, gone.plot)
  
  return(list(models, plots))
}
run_seg(dat_list$Run_05284_005)


for (i in 1:length(dat_list[1:10])) {
  run_seg(dat_list[i][[1]])
}

###### trim dataset to appropriate range
trimmed <- list()
for (i in 1:length(dat_list[1:10])) {
  curr_run <- dat_list[i][[1]]
  
  ###Trim epos data
  epos <- curr_run[[1]]
  ##trim to less than 250
  epos <- epos[epos$X.Time < 250, ]
  ##check for gap before cutoff and if more data outside of range
  if ((tail(epos, n=1)$X.Time < 225) & (nrow(epos) < nrow(curr_run[[1]]))) {
    #add next available datapoint
    newrow <- curr_run[[1]][nrow(epos) + 1, ]
    epos <- rbind(epos, newrow)
  }
  
  ###Trim Stairway data
  stair <- curr_run[[2]]
  ##trim to less than 250
  stair <- stair[stair$year < 250, ]
  ##check for gap before cutoff and if more data outside of range
  if ((tail(stair, n=1)$year < 225) & (nrow(stair) < nrow(curr_run[[2]]))) {
    newrow <- curr_run[[2]][nrow(stair) + 1, ]
    stair <- rbind(stair, newrow)
  }
  
  ###Trim gone data
  gone <- curr_run[[3]]
  ##trim to less than 250
  gone <- gone[gone$Generation < 250, ]
  ##check for gap before cutoff and if more data outside of range
  if ((tail(gone, n=1)$Generation < 225) & (nrow(gone) < nrow(curr_run[[3]]))) {
    newrow <- curr_run[[3]][nrow(gone) + 1, ]
    gone <- rbind(gone, newrow)
  }
  
  ###Add trimmed dataframes to output list
  new_run <- list(epos, stair, gone)
  
}
