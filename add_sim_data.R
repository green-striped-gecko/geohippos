##############################################
#        SIMULATED MODEL FUNCTIONS           #
##############################################

library(tidyr)

###generate expansion model
expan_xy <- function(p0, cp, tl, ts, yr) {
  xn <- p0;
  x1 <- p0 * cp;
  if (yr > (ts - 1)) {
    x <- x1;
  }
  tend <- ts - tl;
  r <- log(xn/x1)/tl;
  if ((yr < ts) & (yr >= tend)) {
    x  <- x1*exp(r*(ts - yr));
  }
  if (yr < tend) {
    x <- xn;
  }
  return(x);
}

###generate decline model
decline_xy <- function(p0, cp, tl, ts, yr) {
  x1 <- p0;
  xn <- p0 * cp;
  if (yr > (ts - 1)) {
    x <- x1;
  }
  tend <- ts - tl;
  r <- log(x1 - xn)/tl;
  if ((yr < ts) & (yr >= tend)) {
    x <- (x1 - xn)*exp((-r)*(ts - yr)) + xn;
  }
  if (yr < tend) {
    x <- xn;
  }
  return(x);
}

###generate bottleneck model
bottle_xy <- function(p0, cp, tl, ts, yr) {
  x1 <- p0;
  xn <- p0 * cp;
  tl1 <- 0.1*tl;
  tl2 <- 0.9*tl;
  tend <- ts - tl;
  r1 <- log(x1-xn)/tl1;
  
  if ((yr >  ts) | (yr < tend)) {
    x <- p0;
  }
  
  if ((yr <= ts) & (yr >= (ts - tl1))){
    ##crash
    x  <- (x1 - xn)*exp(-r1 * (ts - yr)) + xn;
  }
  if ((yr < (ts - tl1)) & (yr > tend)) {
    x1.rec <- xn;
    xn.rec <- x1;
    r <- log(xn.rec/x1.rec)/tl2;
    x  <- x1.rec*exp(r*(ts - tl1 - yr));
  }
  return(x);
}


###run decline for list
model_plot <- function (mod, p0, c.prop, tl, ts, yr) {
  
  
  if (mod == 1) {
    x  <- decline_xy(p0, c.prop, tl, ts, yr)
  }
  else {
    if (mod == 2) {
      x  <- expan_xy(p0, c.prop, tl, ts, yr)
    }
    
    else {
      x  <- p0;
    }
  }
  return(x);
}

dd <- expand_grid(pop.init = c(1000, 500, 200, 100),
                  crash.prop = c(0.1, 0.5),
                  tl = c(100),
                  ts = c(200),
                  model = c(1:3))
                  #yr = 1:250)
dd$run <- 1:nrow(dd)

dd <- expand_grid(dd, yr = 1:250)

for (i in 1:nrow(dd)) {
  dd$Ne[i] <- model_plot(dd$model[i], dd$pop.init[i], dd$crash.prop[i], dd$tl[i], dd$ts[i], dd$yr[i])
}

ggplot(data = subset(dd, model == 3), aes(x = yr, y = Ne, colour = as.factor(crash.prop), group = run)) +
  geom_line() +
  facet_wrap(tl ~ ts)



