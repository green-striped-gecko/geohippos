#######################################
#         HELPER FUNCTIONS            #
#######################################


#Extract data from outputs#

extract_output <- function(df1, df2) {
  if(nrow(df1) != ncol(df2)) {
    print("Error: Dataframe rows not equal to column entries!")
    return()
  };
  res <- data.frame();
  for (i in 1:nrow(df1)) {
    new_data <- cbind(df1[i,c(1:11)], as.data.frame(df2[,i]));
    res <- rbind(res, new_data);
  }
  return(res);
}


df_extract_output <- function(df1, x) {
  res <- data.frame();
  for (i in 1:nrow(df1)) {
    new_data <- cbind(df1$x[[i]], df1[i,c(1:7, 9:10)]);
    res <- rbind(res, new_data);
  }
  return(res);
}