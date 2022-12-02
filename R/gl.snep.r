#' Wrapper function to run SNEP
#' 
#' @export
gl.snep <- function(x, outfile="dummy",snep.path, n.cores=1, args="") {
  
  # check OS
  os <- Sys.info()['sysname'] 
  if (os=="Windows") progs <- "SNep1.11.exe"
  if (os=="Linux") progs <- "SNeP1.1"
  if (os=="Darwin") stop("Currently there is no mac binary for Snep")
  fex <- file.exists(file.path(snep.path, progs))
  if (all(fex)) {
    file.copy(file.path(snep.path, progs),
              to = tempdir(),
              overwrite = TRUE)
  } else{
    cat(
      error(
        "  Cannot find",
        progs[!fex],
        "in the specified folder given by snep.path:",
        snep.path,
        "\n"
      )
    )
    stop()
  }
  snep.path <- tempdir()
  gl2plink(x, outfile = outfile, outpath =snep.path, verbose=0)
  op <- getwd()
  setwd(snep.path)
  snpcmd <- progs
  if (os=="Linux") {
    snpcmd <- paste0("./",snpcmd)
    system(paste0("chmod 555 ",progs))
  }
  system(paste0(snpcmd," -ped ",outfile,".ped -t ", n.cores," ", args))
  
  res <- read.csv(paste0(outfile,".NeAll"), sep="\t", header=T)
  setwd(op)
  return(res)
}
