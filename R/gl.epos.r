#' Runs epos on a genlight/dartR obejct
#' 
#' This is an implementation/wrapper to run epos \url{https://github.com/EvolBioInf/epos}. It relies on compiled version of epos, epos2plot and if a bootstrap output is required bootSfs. For more details on how to install these programs via \url{https://github.com/EvolBioInf/epos} and look out for the manual epos.pdf (in the doc folder of the repository \url{https://github.com/EvolBioInf/epos/blob/master/doc/epos.pdf}. The binaries need to be provided in a single folder. There are compiled versions available and also the necessary dlls for windows need to be provided (Under Linux/MacOS gls, blas need to be installed on your system)
#' @param x dartR/genlight object
#' @param epos.path path to epos and other required programs (epos, epos2plot are always required and bootSfs in case a bootstrap and confidence estimate is required )
#' @param sfs if no sfs is provided function gl.sfs(x, minbinsize=1, singlepop=TRUE) is used to calculate the sfs that is provided to epos 
#' @param minbinsize remove bins from the left of the sfs. if you run epos from a genlight object the sfs is calculated by the function (using gl.sfs) and as default minbinsize is set to 1 (the monomorphic loci of the sfs are removed). This parameter is ignored if sfs is provide via the sfs parameter (see below). Be aware even if you genlight object has more than one population the sfs is calculated with singlepop set to true (one sfs for all individuals) as epos does not work with multidimensional sfs)
#' @param folded if set to TRUE (default) a folded sfs (minor allele frequency sfs) is returned. If set to FALSE then an unfolded (derived allele frequency sfs) is returned. It is assumed that 0 is homozygote for the reference and 2 is homozygote for the derived allele. So you need to make sure your coding is correct. option -U in epos.
#' @param l length of sequences (including monomorphic and polymorphic sites). If the sfs is provided with minbinsize=1 (default) then l needs to be specified. 
#' option -l in epos
#' @param u mutation rate. If not provided the default value of epos is used (5e-9).
#' option -u in epos
#' @param boot if set to a value >0 the programm bootSfs is used to provide multiple
#' bootstrapped sfs, which allows to calculate confidence intervals of the historic Ne
#' sizes. Be aware the runtime can be extended. default:0 no bootstrapped simulations are
#' run, otherwise boot number of bootstraps are run (option -i in bootSfs)
#' @param upper upper quantile of the bootstrap (only used if boot>0). default 0.975.
#' (option -u in epos2plot)
#' @param lower lower quantile of the bootstrap (only used if boot>0). default 0.025.
#' (option -l in epos2plot)
#' @param other.options additional options for epos (e.g -m, -x etc.)
#' @param cleanup if set to true intermediate tempfiles are deleted after the run
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report [default 2,
#' unless specified using gl.set.verbosity].
#' @return returns a table of the Ne over time and a ggplot using estimate of times
#' and Ne (see epos.doc)
#' @export
#' @examples 
#' require(dartR.data)
#' gl.epos(possums.gl[1:60,])
#' gl.msfs(possums.gl[c(1:5,31:33),], minbinsize=1)
#'@references Lynch M, Haubold B, Pfaffelhuber P, Maruki T. Inference of Historical
#' Population-Size Changes with Allele-Frequency Data. G3 (Bethesda). 2020 Jan 7;10(1):211-223. doi: 10.1534/g3.119.400854. PMID: 31699776; PMCID: PMC6945023.
#'@author Custodian: Bernd Gruber


gl.epos <- function(x,epos.path, sfs=NULL, minbinsize=1,folded=TRUE, l=NULL, 
                    u=5e-9, boot=0, upper=0.975, lower=0.025, other.options="",
                    cleanup=TRUE, verbose=NULL)
{
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  #utils.flag.start(func = funname,
  #                 build = "Jody",
  #                 verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # all sanity checks (e.g. if folded length(sfs)=nInd etc.)
  
  # check OS
  os <- tolower(Sys.info()['sysname'] )
  
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  #if (os=="Linux" | os=="Darwin" )  os=="Windows"
  # check if epos epos2plot and [bootSfs are there]
  progs <- c("epos", "epos2plot")
  if (boot>0) progs <- c(progs,"bootSfs")
  if (os=="windows") progs <- paste0(progs,".exe") 
  if (os=="windows") progs <- c(progs, "libblas.dll","libgsl.dll","libgslcblas.dll")
  fex <- file.exists(file.path(epos.path, progs))
  if (all(fex)) {
    file.copy(file.path(epos.path, progs),
              to = tempd,
              overwrite = TRUE)
  } else{
    cat(
      "  Cannot find",
        progs[!fex],
        "in the specified folder given by epos.path:",
        epos.path,
        "\n"
      
    )
    stop()
  }

  if (is.null(sfs)) sfs <- gl.sfs(x, minbinsize=minbinsize,folded=folded,singlepop = TRUE, plot.out = FALSE, verbose = 0 )
  df <- data.frame(r=1:length(sfs), fr=sfs)
  write.table(df,file=file.path(tempd,"dummy.sfs"), row.names = F, sep="\t", col.names = TRUE, quote = FALSE)
  bootcmd<-""
  
  if (boot>0) bootcmd <- paste("bootSfs -i",boot, "dummy.sfs ")
  if (os!="windows") bootcmd <- paste0("./", bootcmd)
  if (minbinsize==0) lcmd="" else lcmd=paste0(" -l ", l)
  ucmd <- paste0(" -u ",u)
  if (boot>0) sfsfile <-"bs.sfs" else sfsfile <- " dummy.sfs"
  eposcmd <- paste0("epos ",lcmd, ucmd,other.options, " ", sfsfile)
  if (os!="windows") eposcmd <- paste0("./", eposcmd)
  # DO THE JOB
  old.path <- getwd()
  setwd(tempd)

  if (boot>0) {
    if (os=="linux") system("chmod 777 bootSfs")
    bsdummy <- system(bootcmd, intern = TRUE)
    writeLines(bsdummy,file.path(tempd,"bs.sfs"))
  }
  if (os=="linux") system("chmod 777 epos")
  epdummy <- system(eposcmd, inter=T)
  writeLines(epdummy,file.path(tempd,"ep.dat"))
  eposplotcmd <-"epos2plot ep.dat"
  if (os=="linux") system("chmod 777 epos2plot")
  if (os!="windows") eposplotcmd <- paste0("./", eposplotcmd)
  eposout <- system(eposplotcmd, intern = TRUE)
  setwd(old.path)
  ep2 <- (do.call(rbind,(strsplit(eposout,split = "\t"))))
  epp <- data.frame(ep2[-1,])
  colnames(epp)<- ep2[1,]
  epp <- data.frame(apply(epp,2, as.numeric))
  
  
  if (cleanup) unlink(tempd, recursive = T)
  
  
  return(epp)
}
  
  
