##############################################################################
#                   Analysis helper functions
##############################################################################


###  Segmented Code

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


for (i in 1:length(dat_list[1:10])) {
  run_seg(dat_list[i][[1]])
}

##extract selected subset to trim


#==============================================================================#
#              Trim dataset to analysis range - single core
#==============================================================================#

trim_single <- function(data_list) {
  trimmed <- list()
  for (i in 1:length(data_list)) {
    print(i)
    curr_run <- data_list[i][[1]]
    ###test more efficient method
    new_run <- curr_run
    
    ### trim dataframes
    new_run[[1]] <- new_run[[1]][new_run[[1]]$X.Time <= 250 #& new_run[[1]]$X.Time >= 1
                                 , ]
    new_run[[2]] <- new_run[[2]][new_run[[2]]$year <= 250 #& new_run[[2]]$year >= 1
                                 , ]
    new_run[[3]] <- new_run[[3]][new_run[[3]]$Generation <= 250  #& new_run[[3]]$Generation >= 1
                                 , ]
    new_run[[4]] <- new_run[[4]][new_run[[4]]$year <= (250 + 1)  #& new_run[[4]]$year >= 1
                                 , ]
    
    ### check dataframes for large gaps - if present, add next available data point

    if (nrow(new_run[[1]]) > 0 & (tail(new_run[[1]]$X.Time, n=1) < 225) & (nrow(new_run[[1]]) < nrow(curr_run[[1]]))) {
      new_run[[1]] <- rbind(new_run[[1]], curr_run[[1]][nrow(new_run[[1]]) + 1, ])
    }
    if (nrow(new_run[[2]]) > 0 & (tail(new_run[[2]]$year, n=1) < 225) & (nrow(new_run[[2]]) < nrow(curr_run[[2]]))) {
      new_run[[2]] <- rbind(new_run[[2]], curr_run[[2]][nrow(new_run[[2]]) + 1, ])
    }
    if (nrow(new_run[[3]]) > 0 & (tail(new_run[[3]]$Generation, n=1) < 225) & (nrow(new_run[[3]]) < nrow(curr_run[[3]]))) {
      new_run[[3]] <- rbind(new_run[[3]], curr_run[[3]][nrow(new_run[[3]]) + 1, ])
    }
    
    ###transfer into output with run name
    trimmed[i] <- data_list[i]
    names(trimmed[i]) <- names(data_list[i])
    trimmed[i][[1]] <- new_run
  }
  names(trimmed) <- names(data_list)
  return(trimmed)
}
#==============================================================================#
#==============================================================================#

#==============================================================================#
#              Trim dataset to analysis range - parallel
#==============================================================================#

trim <- function(data_list, ncores = 10) {
  trimmed <- mclapply(1:length(data_list), function(i) {
    print(i)
    curr_run <- data_list[i][[1]]
    new_run <- curr_run
    ### trim dataframes
    new_run[[1]] <- new_run[[1]][new_run[[1]]$X.Time <= 250 #& new_run[[1]]$X.Time >= 1
                                 , ]
    new_run[[2]] <- new_run[[2]][new_run[[2]]$year <= 250 #& new_run[[2]]$year >= 1
                                 , ]
    new_run[[3]] <- new_run[[3]][new_run[[3]]$Generation <= 250  #& new_run[[3]]$Generation >= 1
                                 , ]
    new_run[[4]] <- new_run[[4]][new_run[[4]]$year <= (250 + 1)  #& new_run[[4]]$year >= 1
                                 , ]
    
    ### check dataframes for large gaps - if present, add next available data point
    
    if (nrow(new_run[[1]]) > 0 & (tail(new_run[[1]]$X.Time, n=1) < 225) & (nrow(new_run[[1]]) < nrow(curr_run[[1]]))) {
      new_run[[1]] <- rbind(new_run[[1]], curr_run[[1]][nrow(new_run[[1]]) + 1, ])
    }
    if (nrow(new_run[[1]][new_run[[1]]$X.Time > 250, ])) {
      new_run[[1]][new_run[[1]]$X.Time > 250, ]$X.Time <- 250
    }
    new_run[[1]][new_run[[1]]$X.Time > 250, ]
    if (nrow(new_run[[2]]) > 0 & (tail(new_run[[2]]$year, n=1) < 225) & (nrow(new_run[[2]]) < nrow(curr_run[[2]]))) {
      new_run[[2]] <- rbind(new_run[[2]], curr_run[[2]][nrow(new_run[[2]]) + 1, ])
    }
    if (nrow(new_run[[2]][new_run[[2]]$year > 250, ])) {
      new_run[[2]][new_run[[2]]$year > 250, ]$year <- 250
    }
    if (nrow(new_run[[3]]) > 0 & (tail(new_run[[3]]$Generation, n=1) < 225) & (nrow(new_run[[3]]) < nrow(curr_run[[3]]))) {
      new_run[[3]] <- rbind(new_run[[3]], curr_run[[3]][nrow(new_run[[3]]) + 1, ])
    }
    if (nrow(new_run[[3]][new_run[[3]]$Generation > 250, ])) {
      new_run[[3]][new_run[[3]]$Generation > 250, ]$Generation <- 250
    }
    
    return(new_run)
  }, mc.cores = ncores)
  names(trimmed) <- names(data_list)
  return(trimmed)
}
#==============================================================================#
#==============================================================================#


#==============================================================================#
#                              Metadata function  
#==============================================================================#
                               
metadata <- function(data_list) {
  meta <- data.frame(run = names(data_list));
  meta <- meta %>% mutate(loci = as.numeric(substr(run, 11, 13)));
  for (i in 1:nrow(meta)) {
    model <- sim_files[meta$run[i]][[1]]$model[3]
    pop_init <- sim_files[meta$run[i]][[1]]$pop_init[3]
    ts <- sim_files[meta$run[i]][[1]]$ts[3]
    tl <- sim_files[meta$run[i]][[1]]$tl[3]
    ss <- sim_files[meta$run[i]][[1]]$ss[3]
    cp <- sim_files[meta$run[i]][[1]]$crash_prop[3]
    meta$model[i] <- model
    meta$pop_init[i] <- pop_init
    meta$ts[i] <- ts
    meta$tl[i] <- tl
    meta$ss[i] <- ss
    meta$cp[i] <- cp
  }
  return(meta)
}
#==============================================================================#
#==============================================================================#



#==============================================================================#
#                   Convert to long dataframe (single core)
#==============================================================================#
##convert data to one long dataframe
list_to_long <- function(data_list) {
  dat_long <- data.frame()
  for (i in 1:length(data_list)) {
    epos <- data.frame(method = "epos", year = data_list[[i]][[1]]$X.Time, Ne = data_list[[i]][[1]]$Median)
    stair <- data.frame(method = "stairway", year = data_list[[i]][[2]]$year, Ne = data_list[[i]][[2]]$Ne_median)
    gone <- data.frame(method = "gone", year = data_list[[i]][[3]]$Generation, Ne = data_list[[i]][[3]]$Geometric_mean)
    sim <- data.frame(method = "sim", year = data_list[[i]][[4]]$year, Ne = data_list[[i]][[4]]$nind)

    newrun <- rbind(epos, stair, gone, sim)
    newrun$model <- data_list[[i]][[4]]$model[1]
    newrun$pop_init <- data_list[[i]][[4]]$pop_init[1]
    newrun$tl <- data_list[[i]][[4]]$tl[1]
    newrun$ts <- data_list[[i]][[4]]$ts[1]
    newrun$ss <- data_list[[i]][[4]]$ss[1]
    newrun$cp <- data_list[[i]][[4]]$crash_prop[1]
    newrun$runnum <- names(data_list[i])
    newrun$loci <- as.numeric(substr(names(data_list[i]), 11, 13))
    
    dat_long <- rbind(dat_long, newrun)
  }
  return(dat_long)
}
#==============================================================================#
#==============================================================================#


#==============================================================================#
#                   Convert to long dataframe - Parallel
#==============================================================================#

par_data_long <- function(data_list, ncores = 10) {
  
  long_list <- mclapply(1:length(data_list), function(i) {
    
    epos <- data.frame(method = "epos", year = data_list[[i]][[1]]$X.Time, Ne = data_list[[i]][[1]]$Median)
    stair <- data.frame(method = "stairway", year = data_list[[i]][[2]]$year, Ne = data_list[[i]][[2]]$Ne_median)
    gone <- data.frame(method = "gone", year = data_list[[i]][[3]]$Generation, Ne = data_list[[i]][[3]]$Geometric_mean)
    sim <- data.frame(method = "sim", year = data_list[[i]][[4]]$year, Ne = data_list[[i]][[4]]$nind)
    
    newrun <- rbind(epos, stair, gone, sim)
    newrun$model <- data_list[[i]][[4]]$model[1]
    newrun$pop_init <- data_list[[i]][[4]]$pop_init[1]
    newrun$tl <- data_list[[i]][[4]]$tl[1]
    newrun$ts <- data_list[[i]][[4]]$ts[1]
    newrun$ss <- data_list[[i]][[4]]$ss[1]
    newrun$cp <- data_list[[i]][[4]]$crash_prop[1]
    newrun$runnum <- names(data_list[i])
    newrun$loci <- as.numeric(substr(names(data_list[i]), 11, 13))
    return(newrun)
  }, mc.cores = ncores)
  names(long_list) <- names(data_list)
  out <- rbindlist(long_list)
  return(data.frame(out))
}
#==============================================================================#
#==============================================================================#


#==============================================================================#
#                  Extract full loci data from vcf files
#==============================================================================#

library(geohippos)
library(dartR)
library(parallel)
fn_loc <- list.files("/data/scratch/isobel/vcf/", pattern = "*.vcf", full.names = TRUE)

#for (i in 1:length(fn)) {

glsnloc <- function(i)
{
  gls <- geohippos::gl.read.vcf(fn_loc[i], verbose = 0)
  #pop(gls)<- rep("A", nInd(gls))
  #gls$chromosome <- factor(ceiling(gls$position/1e8)) #slim simulation
  #gls@chromosome <- factor("1", )
  return(nLoc(gls))
}

ss <- mclapply(1:length(fn_loc), function(x) glsnloc(x), mc.cores  =25)

ss_loc <- lapply(1:length(fn_loc), function(x) glsnloc(x))

names(ss_all) <- substr(fn_loc, 28, 36)

data.frame(unlist(ss_loc))
