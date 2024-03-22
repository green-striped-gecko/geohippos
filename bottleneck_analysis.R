#=============================================================================#
#                         Bottleneck sample analysis
#=============================================================================#

library(geohippos)
library(data.table)
library(viridis)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
###Extract standard dataset for analysis####

###set results folder
res_folder1 <- "/data/scratch/isobel/bottle_res1/"
res_folder2 <- "/data/scratch/isobel/bottle_res2/"
res_folder3 <- "/data/scratch/isobel/bottle_res3/"


###set logfile folder
log_folder <- "/data/scratch/isobel/bottle_log/"

fn <- list.files(res_folder1, pattern = "*.rds", full.names = TRUE)
fn2 <- list.files(res_folder2, pattern = "*.rds", full.names = TRUE)
fn3 <- list.files(res_folder3, pattern = "*.rds", full.names = TRUE)

## Load data for sample ##
dat_list = sapply(fn3, function (x) data.table(readRDS(x)))

## Convert list names to run numbers ##
names(dat_list) <- unlist(lapply(names(dat_list), function(x) substr(x, 42, 54)))

# ## Generate list of potential files from folder
# sim_list <- list.files(log_folder, pattern = "*.txt", full.names = TRUE)
# 
# ##Generate list of simulated data
# sim_files <- lapply(sim_list, function(x) read.delim2(x, header = T, sep = ","))
# 
# #extract model and pop_init for runname
# names(sim_files) <- lapply(sim_files, function(x) {
#   name_list <- paste0(x$model[1], x$pop_init[1])
#   return(name_list)
# })
# 
# ##import run_list dataframe from sim script and change runnum name
# df <- df %>% rename(runnum = runnumb)
# 
# 
# #add loci data to dataframe
# df <-  left_join(df, loci_df[, c(1:2, 15)], by = c("model", "pop_init"), multiple = "first")
# df <- df %>% rename(loci_total = loci, runnum = runnumb)

###Add simulated data to main data list
for (i in 1:length(dat_list)) {
  id <- grep(substr(names(dat_list)[i], 1, 9), df$runnum)
  code <- paste0(df$model[id], df$pop_init[id])
  
  dat_list[i][[1]][4] <- sim_files[code]
  dat_list[i][[1]][[4]]$year <- 2301 - dat_list[i][[1]][[4]]$cycle
}

data_trimmed <- trim(dat_list, ncores = 2, trim_pt = 350)

data_long <- par_bottle_long(data_trimmed, ncores = 5)

#add runnum column without loci
data_long$runnum2 <- substr(data_long$runnum, 1, 9)

##add sample size column by runnum
dl3 <- left_join(data_long, df[, c(12:13, 15)], by = c("runnum2" = "runnum"))

dl1$set <- "set1"
dl2$set <- "set2"
dl3$set <- "set3"

data <- rbind(dl1, dl2, dl3)

#calculate conversion factor
data$c_factor <- data$loci_total/(data$loci * 1000)

#recalculate year and Ne
data <- data %>% 
  mutate("c_year" = if_else(method %in% c("epos", "stairway"), year * c_factor, year),
         "c_Ne" = if_else(method %in% c("epos", "stairway"), Ne * c_factor, Ne))




###add runcode column to data
data$runcode <- paste0(data$model, data$pop_init)
df$runcode <- paste0(df$model, df$pop_init)
data$run_set <- paste0(data$runnum, data$set)
#=============================================================================#
#                   Plot all combinations 
#=============================================================================#

methods <- c("epos", "stairway", "gone")

combinations <- expand.grid("runcode" = unique(df$runcode), "method" = methods)

# Loop over each combination
for (i in 1:nrow(combinations)) {
  combo_data <- subset(data, runcode == combinations$runcode[i] & method == combinations$method[i])
  sim_data <- subset(data, runcode == combinations$runcode[i] & method == "sim")
  
  # Check if plot_data is empty
  if (nrow(combo_data) == 0) {
    # cat("Skipping empty plot for combination:", combinations[i, "pop_init"], combinations[i, "cp"], 
    #     combinations[i, "tl"], combinations[i, "ts"], "\n")
    next  # Skip to the next combination
  }
  
  # Plot year against Ne with facets for each combination of loci and ss
  p <- ggplot(combo_data, aes(x = c_year, y = c_Ne, group = run_set, colour = as.character("Ne Method"))) +
    geom_line(linewidth = 0.2) +
    geom_line(data = sim_data, aes(x = year, y = Ne, colour = as.character("Sim")), linewidth = 0.2) +
    facet_grid(loci ~ ss, scales = "free", labeller = label_both)+
    ylim(c(0, 2*max(sim_data$c_Ne))) +
    xlim(c(0, 500)) +
    labs(subtitle = paste("Year vs Ne: ", combinations$runcode[i], "_", combinations$method[i]),
         x = "Year",
         y = "Ne") +
    theme_minimal() +
    scale_colour_manual(values = c("Ne Method" = "blue", "Sim" = "red"), name = "Method")
  
  # Save the plot
  filename <- paste0("/data/scratch/isobel/bottle_plots/plot_corrected", combinations$runcode[i], "_", combinations$method[i], ".png")
  tryCatch(ggsave(filename, plot = p, width = 30, height = 15, units = "cm"),
           error = function(e) {
             cat("Error occurred while saving plot:", filename, "\n")
             cat("Error message:", conditionMessage(e), "\n")
           })
}

write_rds(data, file = "/data/scratch/isobel/bottle_data.rds")



interp_list <- split(data, list(data$run_set, data$method))

####calculate error
library(parallel)
interp_out_test <-mclapply(1:length(interp_list[1:10]), function(i) {
  x <- interp_list[[i]]
  df <- data.frame(approx(x$c_year, x$c_Ne, xout = 1:350))
  colnames(df) <- c("year", "Ne")
  df[, 3:9] <- x[2, c(1, 4:7, 9, 11)]
  return(df);
},
mc.cores = 5)

names(interp_out) <- names(interp_list)

###for each run/set - extract simulation dataset for comparison
grep(paste0(substr(names(interp_out)[2000], 1, 18), "sim"), names(interp_out))

interp_out[3887][[1]]



####Calculate error for bottleneck plots
bottle_rmse <- lapply(1:length(interp_out), function(i) {
  run_set <- substr(names(interp_out[i]), 1, 18)
  method <- substr(names(interp_out[i]), 19, 30)
  curr_data <- interp_out[i][[1]]
  sim_data <- interp_out[grep(paste0(run_set, "sim"), names(interp_out))][[1]]
  n <- nrow(curr_data)
  all_data <- left_join(curr_data, sim_data[, 1:2], by = "year")
  all_data$error_sq <- ((all_data$Ne.y - all_data$Ne.x)^2)/n
  rmse <- all_data %>% summarise(RMSE = sqrt(sum(error_sq)))
  out <- data.frame("method" = method, rmse, curr_data[1, 3:11])
  return(out)
})

##combine output into one fataframe
df_b_rmse <- do.call("rbind", bottle_rmse)

df_b_rmse$run_code <- paste0(df_b_rmse$model, df_b_rmse$pop_init)

write_csv(df_b_rmse, file = "/data/scratch/isobel/bottle_rmse_df.csv")

method_list <- c("epos", "stairway", "gone")

model_list <- list("croc", "frog", "seal", "whale")

lapply(model_list, function(x) {
  p <- ggplot(data = df_b_rmse %>% filter(model == x, method != "sim") , aes(x = as.factor(ss), y = log(RMSE), colour = as.factor(run_code))) +
    geom_boxplot() +
    theme_minimal() +
    facet_grid(method ~ loci, labeller = "label_both", scales = "free") +
    labs(
      x = "Sample Size",
      y = "Root mean square error") +
    scale_colour_viridis(discrete = TRUE, option = "A", end = 0.8)
  return(p)
})


###################################################################
# Whale8000 not in final bottle_rmse dataframe
# Whale8000 runnumbers: Run_00241 - Run_00280
ggplot(data = df_b_rmse %>% filter(method != "sim") , aes(x = as.factor(ss), y = log(RMSE), colour = as.factor(run_code))) +
  geom_boxplot() +
  theme_minimal() +
  facet_grid(method ~ loci, labeller = "label_both", scales = "free") +
  labs(
    x = "Loci",
    y = "Root mean square error") +
  scale_colour_viridis(discrete = TRUE, option = "A", end = 0.8)
#ylim(c(0, 25000))
  

##plot by individual parameters
ggplot(data = df_b_rmse %>% filter(method != "sim", model == "croc") , aes(x = as.factor(ss), y = log(RMSE), colour = as.factor(model))) +
  geom_boxplot() +
  theme_minimal() +
  facet_grid(method ~ loci, labeller = "label_both", scales = "free") +
  labs(
    x = "Loci",
    y = "Root mean square error") +
  scale_colour_viridis(discrete = TRUE, option = "A", end = 0.8)






