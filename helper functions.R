#######################################
#         HELPER FUNCTIONS            #
#######################################


#Extract data from outputs#

extract_output <- function(df1, df2) {
  if(nrow(df1) != length(df2)) {
    print("Error: Dataframe rows not equal to column entries!")
    return()
  };
  res <- data.frame();
  for (i in 1:nrow(df1)) {
    new_data <- cbind(df1[i,c(1:11)], as.data.frame(df2[[i]]));
    res <- rbind(res, new_data);
  }
  return(res);
}


df_extract_output <- function(df1, out_index, keep) {
  x <- out_index;
  res <- data.frame();
  for (i in 1:nrow(df1)) {
    new_data <- cbind(df1[i,x][[1]], df1[i, keep]);
    res <- rbind(res, new_data);
  }
  return(res);
}



#==============add simulated ind filename to dataframe===============================

modelindex <- data.frame(model = c("decline", "expansion", "bottle", "stable"), n = 1:4)


test.epos <- test.epos %>% 
  mutate(simindfile = paste0("C:/Users/Isobel/Desktop/Honours/geohippos/sim_ind_folder/sim_ind_", 
                             match(model, modelindex$model), 
                             pop_init, crash_prop, ts, tl,".txt"))
test.epos$simindfile

##read in simind data
for (i in 1:nrow(test.epos)) {
  test.epos$siminddata[[i]] <- read.delim2(test.epos$simindfile[[i]], header = T, sep = ",")
}

#==============Extract sim data into usable dataframe=============================#


extract_sim <- function(df) {
  sim_data <- data.frame()
  for (i in 1:nrow(df)) {
    sim_data <- rbind(sim_data, df$siminddata[[i]]);
  }
  return(sim_data);
}

