#=============================================================================#
#                         Failed runs
#=============================================================================#

#Gather all not run filenames
notrun1 <- read.csv(file = "/data/scratch/isobel/notrun1/home/589/bmg589/gls/groups.csv")
notrun2 <- read.csv(file = "/data/scratch/isobel/notrun2/groups.csv")

std_nr <- rbind(notrun1, notrun2)


#Gather all bottle not run files
bottle1 <- read.csv(file = "/data/scratch/isobel/bottle_notrun1/groups.csv")
bottle2 <- read.csv(file = "/data/scratch/isobel/bottle_notrun2/groups.csv")
bottle3 <- read.csv(file = "/data/scratch/isobel/bottle_notrun3/groups.csv")

bottle_nr <- rbind(bottle1, bottle2, bottle3)


##extract run_code and loci
std_nr <- std_nr  %>% mutate(run_code = 
                     if_else(grepl("Run_", fn) == TRUE, substr(fn, 5, 9), substr(fn, 6, 10)),
                     loci = if_else(grepl("Run_", fn) == TRUE, substr(fn, 11, 13), substr(fn, 12,14))) 

##add run_code column to std_sim dataframe
df$run_code <- substr(df$runnumb, 6, 10)

#join dataframes

std_nr <- left_join(std_nr, df, by = "run_code")

##save dataframe
write_csv(std_nr, file = "/data/scratch/isobel/std_notrun.csv")


#####Analyse interp_list for short runs

lapply(interp_list[1:300], function(X) {
  if (max(X$c_year) > 2000) {
    return(NULL)
  }
  else(return(max(X$c_year)) )
  
})

sim_data_all <- new_data %>% filter(method == "sim")

max_years <- new_data %>% group_by(runnum, method) %>% summarise(max_yr = max(c_year))

low_max <- max_years %>% filter(max_yr < 225 & max_yr > 150)

###extract plot_data for short inference runs
short_runs <- left_join(low_max, new_data, by = c("runnum", "method"))
short_runs$cp <- as.numeric(short_runs$cp)

##extract unique cases

short_alldata <- left_join(low_max, new_data[, c(1, 4:11)], by = c("runnum", "method"), multiple = "first")


short_alldata2 <- short_alldata %>% ungroup() %>% select(-c(ss, loci, runnum, max_yr))

short_alldata2 <- unique(short_alldata2)

combo <- short_alldata2
short_runs %>% filter(model == short_alldata2$model[3] & pop_init == short_alldata2$pop_init[3])


############################################################################
#                      Plot all combos for short runs

for (i in 1:nrow(combo)) {
  combo_data <- short_runs %>% filter(pop_init == combo$pop_init[i] &
                                        cp == combo$cp[i] &
                                        tl == combo$tl[i] &
                                        ts == combo$ts[i] &
                                        model == combo$model[i])
  
  sim_data <- sim_data_all %>% filter(pop_init == combo$pop_init[i] &
                                        cp == combo$cp[i] &
                                        tl == combo$tl[i] &
                                        ts == combo$ts[i] &
                                        model == combo$model[i])
  if (nrow(combo_data) == 0) {
    print(paste("Skipping empty plot for combination: ", i))
    next 
  }
  plot_data <- combo_data %>% filter(method == combo$method[i])
  p <- ggplot(plot_data, aes(x = c_year, y = c_Ne, group = runnum, colour = as.character("Ne Method"))) +
    geom_line(linewidth = 0.2) +
    geom_line(data = sim_data, aes(x = year, y = Ne, colour = as.character("Sim")), linewidth = 0.2) +
    facet_grid(loci ~ ss, scales = "free", labeller = label_both)+
    ylim(c(0, max(sim_data$Ne) + 500)) +
    xlim(c(0, 250)) +
    labs(subtitle = paste("Year vs Ne (pop_init =", combo$pop_init[i], ", cp =", combo$cp[i],
                          ", tl =", combo$tl[i], ", ts =", combo$ts[i], "model = ", combo$model[i], "method =", combo$method[i], ")"),
         x = "Year",
         y = "Ne") +
    theme_minimal() +
    scale_colour_manual(values = c("Ne Method" = "blue", "Sim" = "red"), name = "Method")
  
  filename <- paste0("/data/scratch/isobel/notrun_plots/plot_", combo$pop_init[i], "_", combo$cp[i],
                     "_tl_", combo$tl[i], "_ts_", combo$ts[i], "_model_", combo$model[i], "_method_", combo$method[i], ".png")
  ggsave(filename, plot = p)
  
}



short_fails <- max_years %>% filter(max_yr <= 150)

short_fails <- left_join(short_fails, new_data[, c(1, 4:11)], by = c("runnum", "method"), multiple = "first")

