#' Wrapper function to run GONE
#' 
#' @export 


gl.gone <- function(x,outfile="dummy" ,gone.path)  #need to be able to change parameters....
{
  # check OS
  os <- Sys.info()['sysname'] 
  progs <- c("script_Gone.sh","PROGRAMMES","INPUT_PARAMETERS_FILE")
  fex <- file.exists(file.path(gone.path, progs))
  if (all(fex)) {
    file.copy(file.path(gone.path, progs),
              to = tempdir(),
              overwrite = TRUE, recursive = TRUE)
  } else{
    cat(
      error(
        "  Cannot find",
        progs[!fex],
        "in the specified folder given by gone.path:",
        gone.path,
        "\n"
      )
    )
    stop()
  }
  gone.path <- tempdir()
  
  
  
  
op <- getwd()
setwd(gone.path)
if (os=="Linux") {
  system("chmod 777 script_Gone.sh")
  system("chmod 777 ./PROGRAMMES/*")
}

#TODO: create INPUT PARAMETER FILE
gl2plink(x, outfile = outfile, outpath =gone.path, phen_value=-9, verbose=0)

#SET PATH of INPUT PARAMETER FILE in script_GONE.sh
shell(paste0("script_Gone.sh ",outfile))

setwd(op)
res<-  read.csv(file.path(gone.path,paste0("Output_Ne_",outfile)), sep="\t",skip = 1, header = T)
return(res)
}
  
  