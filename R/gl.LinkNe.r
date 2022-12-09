#' Runs LinkNe on a genlight/dartR obejct
#' 
#' This is an implementation/wrapper to run LinkNe \link[Hollenbeck et al. 2016]{http://www.nature.com/hdy/journal/v117/n4/abs/hdy201630a.html}.
#' @param x genlight/dartR object 
#' @param outfile name of the output file
#' @param L length of a single chromosome (used to create the map file)
#' @param LinkNe.path path to LinkNe (could be a binary (e.g. LinkNe.exe) or the perl version (LinkNe.pl))
#' @param perl switch to run the perl version instead of the binary. To run the perl version your system must have perl installed and accessible from any folder (e.g. in the path under windows)
#' @param other.args additional arguments passed to LinkNe (\link[LinkNe settings]{https://github.com/chollenbeck/LinkNe} )
#' @verbose verbosity (though we cannot suppress console output during LinkNe runs)
#' @export
#' @author Custodian: Bernd Gruber


gl.LinkNe <- function(x, outfile="LinkNe", L=1e8, LinkNe.path,  perl=FALSE, temp.path=tempdir(), other.args="", verbose=1)
{
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # check if temp.path exists, create or stop it
  if (!dir.exists(temp.path)) {
    cc <- dir.create(temp.path)
  }
  # check OS
  os <- Sys.info()['sysname'] 
  progs <- c("LinkNe")
  if (os=="Windows") progs <- paste0(progs,".exe") 
  if (perl) progs <- "LinkNe.pl" 
  fex <- file.exists(file.path(LinkNe.path, progs))
  if (all(fex)) {
    file.copy(file.path(LinkNe.path, progs),
              to = temp.path,
              overwrite = TRUE)
    if (os=="Linux" | os=="Darwin") system(paste("chmod 777 ", file.path(temp.path, progs)))
  } else{
    cat(
      cat(
        "  Cannot find",
        progs[!fex],
        "in the specified folder given by LinkNe.path:",
        LinkNe.path,
        "\n"
      )
    )
    stop()
  }
  old.path <- getwd()  
  setwd(temp.path)  
  
fn="dummy"
gl2genepop(x, outfile=paste0(fn,".gen"), outpath = temp.path, output_format = "3_digits", verbose = verbose)
gl2plink(x, outfile = fn, outpath = temp.path, verbose = verbose)
mm <- read.csv(file.path(temp.path,paste0(fn,".map")), sep=" ", header=F)
ff <- data.frame(locus=mm[,2], chromosome=mm[,1], position=((mm[,4]-(as.numeric(x@chromosome)-1)*L)/1005034))
write.table(ff, file.path(temp.path,paste0(fn,".map2")), row.names = FALSE, quote = F, sep="\t",eol="\r\n")
if (perl) progs <- "perl LinkNe.pl"
if (os!="Windows") progs <- paste0("./",progs)
system(paste0(progs," -i ",fn,".gen -map ", fn,".map2 -o ",outfile," ",other.args))
Ne_LinkNe <- read.csv(paste0(file.path(temp.path,outfile)), sep="\t", header=T)
setwd(old.path)
return(Ne_LinkNe)
}


