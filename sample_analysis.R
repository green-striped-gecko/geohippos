library(geohippos)
library(data.table)
library(viridis)
library(dplyr)
###Extract standard dataset for analysis####

###set results folder
res_folder1 <- "/data/scratch/isobel/results_gadi/"
res_folder2 <- "/data/scratch/isobel/results1"
res_folder3 <- "/data/scratch/isobel/results2/"

###set logfile folder
log_folder <- "/data/scratch/isobel/log_std_set2/"

fn <- list.files(res_folder1, pattern = "*.rds", full.names = TRUE)
fn2 <- list.files(res_folder2, pattern = "*.rds", full.names = TRUE)

##generate random sample to work on
fn_ss <- sample(fn, 1000)

## Load data for sample ##
dat_list = sapply(fn, function (x) data.table(readRDS(x)))

## Convert list names to run numbers ##
names(dat_list) <- unlist(lapply(names(dat_list), function(x) substr(x, 36, 48)))

## Generate list of potential files from folder
sim_list <- list.files(log_folder, pattern = "*.txt", full.names = TRUE)

##Generate list of simulated data
sim_files <- mclapply(names(dat_list), function(x) {
  x <- read.delim2(grep(paste0("Run2_", substr(x, 5, 9)), sim_list, value = T), header = T, sep = ",")}, 
  mc.cores = 20)

sim_files
names(sim_files) <- names(dat_list)
dat_list


###Add simulated data to main data list
for (i in 1:length(dat_list)) {
  dat_list[i][[1]][4] <- sim_files[i]
  dat_list[i][[1]][[4]]$year <- 2201 - dat_list[i][[1]][[4]]$cycle
}


##generate subset data_list using selected metadata subset
meta_sub <- metadata %>% filter(pop_init == 1000, ts == 200, tl == 100, cp == 0.5, loci == 50)
dat_sub <- list()
for (i in 1:nrow(meta_sub)) {
  dat_sub[i] <- dat_list[meta_sub$run[i]]
}
names(dat_sub) <- meta_sub$run

###subset trimmed dataset to only necessary runs
data_trimmed <- trim(dat_list, ncores = 20)

##convert data to one long dataframe
data_long <- par_data_long(data_trimmed, ncores = 20)

method_list <- split(data_long, data_long$method)
sim_data <- method_list[3]$sim
method_list <- method_list[c(1,2,4)]

# Assuming combinations is your dataframe

#Get combinations from std_simulation set
combinations <- unique(df[, c(1:4, 6)])

combinations <- combinations %>%
  mutate(model = recode(model, "decline" = 1, "expansion" = 2, "stable" = 3))


library(ggplot2)

# Assuming dat_long is your dataset
# Convert pop_init, tl, ts, crash_plot, and method to factors to ensure correct ordering

dat_long$method <- factor(dat_long$method)


# Loop over each combination
for (i in 1:nrow(combinations)) {
  plot_data <- subset(dat_long, pop_init == combinations[i, "pop_init"] & crash_plot == combinations[i, "crash_plot"] & 
                        tl == combinations[i, "tl"] & ts == combinations[i, "ts"])
  
  # Check if plot_data is empty
  if (nrow(plot_data) == 0) {
    cat("Skipping empty plot for combination:", combinations[i, "pop_init"], combinations[i, "crash_plot"], 
        combinations[i, "tl"], combinations[i, "ts"], "\n")
    next  # Skip to the next combination
  }
  
  # Get unique combinations of loci and ss
  loci_ss_combinations <- expand.grid(loci = levels(dat_long$loci), ss = levels(dat_long$ss))
  
  # Loop over each method
  for (method in levels(dat_long$method)) {
    # Plot year against Ne with facets for each combination of loci and ss
    p <- ggplot(plot_data, aes(x = year, y = Ne)) +
      geom_point() +
      geom_line(data = subset(plot_data, method == method), aes(group = method), color = "black") +
      facet_grid(loci ~ ss, scales = "free", space = "free") +
      labs(title = paste("Year vs Ne (pop_init =", combinations[i, "pop_init"], ", crash_plot =", combinations[i, "crash_plot"],
                         ", tl =", combinations[i, "tl"], ", ts =", combinations[i, "ts"], ", method =", method, ")"),
           x = "Year",
           y = "Ne") +
      theme_minimal()
    
    # Save the plot
    filename <- paste0("plot_", combinations[i, "pop_init"], "_", combinations[i, "crash_plot"], "_", 
                       combinations[i, "tl"], "_", combinations[i, "ts"], "_", method, ".png")
    tryCatch(ggsave(filename, plot = p),
             error = function(e) {
               cat("Error occurred while saving plot:", filename, "\n")
               cat("Error message:", conditionMessage(e), "\n")
             })
  }
}
