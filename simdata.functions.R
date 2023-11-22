################################################
#==== ADD SIMULATED DATA TO ESTIMATIONS ========  
################################################





#=========== EPOS METHOD ======================
# epos time variable = X.Time, measured in ybp
#epos_simpop function performs on a single line, must be called in a for loop

#designed to be run over converted output file

epos_simpop <- function(df) {
  yr <- (1200 - df$X.Time)
  ts <- df$ts
  tl <- df$tl
  startyr <- 1200 - ts
  endyr <- 1200 - ts + tl
  if (df$model == "decline") {
    p1 <- df$pop_init
    pn <- df$crash_prop * df$pop_init
    r <- log(p1 - pn)/tl
    popout <- if_else(yr < startyr,
                      p1,
                      (p1-pn)*exp(-r * (yr - startyr)) + pn)
  }
  if (df$model == "expansion") {
    p1 <- df$crash_prop * df$pop_init
    pn <- df$pop_init
    r <- log(pn/p1)/tl
    popout <- if_else(yr < startyr, p1, 
                      if_else(p1*exp(r*(yr - startyr)) > pn, pn, p1*exp(r*(yr - startyr))))
  }
  if (df$model == "bottle") {
    p1 <- df$pop_init
    pn <- df$crash_prop * df$pop_init
    r1 <- log(p1 - pn)/(0.1*tl)
    r2 <- log(p1/pn)/(0.9*tl)
    popout <- if_else(yr < startyr, p1,
                      if_else(yr < (startyr + 0.1*tl), (p1-pn)*exp(-r1 * (yr - startyr)) + pn,
                              if_else(pn*exp(r2*(yr - startyr)) > p1), p1, pn*exp(r2*(yr - startyr))))
  }
  if (df$model == "stable") {
    popout <- df$pop_init
  }
  return(popout)
}

for (i in 1:nrow(trialres)) {
  trialres$sim[[i]] <- epos_simpop(trialres[i,])
}










#============ STAIRWAY METHOD ========================


sw_simpop <- function(df) {
  yr <- (1200 - df$X.Time)
  ts <- df$ts
  tl <- df$tl
  startyr <- 1200 - ts
  endyr <- 1200 - ts + tl
  if (df$model == "decline") {
    p1 <- df$pop_init
    pn <- df$crash_prop * df$pop_init
    r <- log(p1 - pn)/tl
    popout <- if_else(yr < startyr,
                      p1,
                      (p1-pn)*exp(-r * (yr - startyr)) + pn)
  }
  if (df$model == "expansion") {
    p1 <- df$crash_prop * df$pop_init
    pn <- df$pop_init
    r <- log(pn/p1)/tl
    popout <- if_else(yr < startyr, p1, 
                      if_else(p1*exp(r*(yr - startyr)) > pn, pn, p1*exp(r*(yr - startyr))))
  }
  if (df$model == "bottle") {
    p1 <- df$pop_init
    pn <- df$crash_prop * df$pop_init
    r1 <- log(p1 - pn)/(0.1*tl)
    r2 <- log(p1/pn)/(0.9*tl)
    popout <- if_else(yr < startyr, p1,
                      if_else(yr < (startyr + 0.1*tl), (p1-pn)*exp(-r1 * (yr - startyr)) + pn,
                              if_else(pn*exp(r2*(yr - startyr)) > p1), p1, pn*exp(r2*(yr - startyr))))
  }
  if (df$model == "stable") {
    popout <- df$pop_init
  }
  return(popout)
}







#=============== GONE METHOD ===========================
# GONE time variable = generations, must be converted to yrs

gone_simpop <- function(df) {
  yr <- df$generation)
  ts <- df$ts
  tl <- df$tl
  startyr <- 1200 - ts
  endyr <- 1200 - ts + tl
  if (df$model == "decline") {
    p1 <- df$pop_init
    pn <- df$crash_prop * df$pop_init
    r <- log(p1 - pn)/tl
    popout <- if_else(yr < startyr,
                      p1,
                      (p1-pn)*exp(-r * (yr - startyr)) + pn)
  }
  if (df$model == "expansion") {
    p1 <- df$crash_prop * df$pop_init
    pn <- df$pop_init
    r <- log(pn/p1)/tl
    popout <- if_else(yr < startyr, p1, 
                      if_else(p1*exp(r*(yr - startyr)) > pn, pn, p1*exp(r*(yr - startyr))))
  }
  if (df$model == "bottle") {
    p1 <- df$pop_init
    pn <- df$crash_prop * df$pop_init
    r1 <- log(p1 - pn)/(0.1*tl)
    r2 <- log(p1/pn)/(0.9*tl)
    popout <- if_else(yr < startyr, p1,
                      if_else(yr < (startyr + 0.1*tl), (p1-pn)*exp(-r1 * (yr - startyr)) + pn,
                              if_else(pn*exp(r2*(yr - startyr)) > p1), p1, pn*exp(r2*(yr - startyr))))
  }
  if (df$model == "stable") {
    popout <- df$pop_init
  }
  return(popout)
}