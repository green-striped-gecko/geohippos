#' Wrapper function to run GONE
#' 
#' @export 


gl.gone <- function(x,outfile="dummy" ,gone.path, cleanup=TRUE)  #need to be able to change parameters....
{
  # check OS
  os <- Sys.info()['sysname'] 
  #create individual tempdir
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  progs <- c("script_Gone.sh","PROGRAMMES")
  fex <- file.exists(file.path(gone.path, progs))
  if (all(fex)) {
    file.copy(file.path(gone.path, progs),
              to = tempd,
              overwrite = TRUE, recursive = TRUE)
  } else{
    cat("  Cannot find",
        progs[!fex],
        "in the specified folder given by gone.path:",
        gone.path,
        "\n")
    stop()
  }
  gone.path <- tempd

op <- getwd()
setwd(gone.path)
if (os=="Linux") {
  system("chmod 777 script_Gone.sh")
  system("chmod 777 ./PROGRAMMES/*")
}

#TODO: create INPUT PARAMETER FILE
gl2plink(x, outfile = outfile, outpath =gone.path, phen_value=-9, verbose=0)

#SET PATH of INPUT PARAMETER FILE in script_GONE.sh
system(paste0("./script_Gone.sh ",outfile))


res<-  read.csv(paste0("Output_Ne_",outfile), sep="\t",skip = 1, header = T)
setwd(op)

if (cleanup) unlink(tempd, recursive = T)
return(res)
}
  

  