##############################################
#        SIMULATED MODEL FUNCTIONS           #
##############################################


x  <- (x1 - xn)*exp(-r * (yr)) + xn;





###generate function for decline model
###displayed in ybp: (ts - yr)

#==========sim script code
# x1 = pop_init;
# xn = cp * pop_init;
# yr = sim.cycle - startyr;
# r = log(x1 - xn)/tl;
# p1.setSubpopulationSize(asInteger(round(
#   (x1 - xn)*exp((-r) * yr) + xn)
# )); 

#==== forward simulation, x = yr
exp.x$decay <- if_else(exp.x$x < startyr,
                       y1,
                       (y1 - yn)*exp(-r * (exp.x$x - startyr)) + yn)

#==== reverse simulation, x = ybp
exp.x$decay <- if_else((1200 - exp.x$ybp) < startyr,
                       y1,
                       (y1 - yn)*exp(-r * ((1200 - exp.x$ybp) - (startyr))) + yn)

exp.x$growth

###for output plots, after data put into long form
###Epos uses X.Time = ybp

epos_simpop <- function(df) {
  for (i in 1:nrow(df)) {
    yr <- df$yr[i]
    ts <- df$ts[i]
    tl <- df$tl[i]
    startyr <- 1200 - ts
    endyr <- 1200 - ts + tl
    if (df$model == "decline") {
      p1 <- df$pop_init[i]
      pn <- df$crash_prop[i] * df$pop_init[i]
      r <- log(p1 - pn) / tl
      popout <- if_else(yr < startyr,
                        p1,
                        (p1 - pn) * exp(-r * (yr - startyr)) + pn)
    }
    if (df$model == "expansion") {
      p1 <- df$crash_prop[i] * df$pop_init[i]
      pn <- df$pop_init[i]
      r <- log(pn / p1) / tl
      popout <- if_else(yr < startyr, p1,
                        if_else(p1 * exp(r * (yr - startyr)) > pn, pn, p1 *
                                  exp(r * (yr - startyr))))
    }
    if (df$model == "bottle") {
      p1 <- df$pop_init[i]
      pn <- df$crash_prop[i] * df$pop_init[i]
      r1 <- log(p1 - pn) / (0.1 * tl)
      r2 <- log(p1 / pn) / (0.2 * tl)
      popout <- if_else(yr < startyr, p1,
                        if_else(
                          yr < (startyr + 0.1 * tl),
                          (p1 - pn) * exp(-r1 * (yr - startyr)) + pn,
                          if_else(
                            yr < (startyr + 0.8 * tl),
                            pn,
                            if_else(pn * exp(r2 * (yr - startyr)) > p1),
                            p1,
                            pn * exp(r2 * (yr - startyr))
                          )
                        ))
      
    }
    if (df$model == "stable") {
      popout <- df$pop_init[i]
    }
    df$sim_ne[i] <- popout
  }
  return(df)
}



for (i in 1:nrow(trialres)) {
  trialres$sim[[i]] <- epos_simpop(trialres[i,])
}


decline <- function(df) {
  outdf <- 
    
    if (yr <= ts) {
      x <- p0;
    }
  tend <- ts - tl;
  x1 <- p0;
  xn <- p0 * cp - 1;
  r <- log(x1 - xn)/tl;
  if (yr < ts) {
    x  <- (x1 - xn)*exp(-r * (ts-yr)) + xn;
  }
  return(x)
}


g <- data.frame(yr = 1:300)
g$d <- lapply(g$yr, function(x) decline_xy(200, 0.1, 100, 200, x))
plot(g$yr, g$d, type = "l")

##decay
exp.x <- data.frame(x = 0:300)
exp.x
tl <- 200
y1 <- 500
yn <- 50
r = log(y1 - yn)/tl

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
  return(x)
}

test$model[9876]

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
  return(x)
}
g <- data.frame(yr = 1:1000)


res_epos %>% select()

for (i in 1:nrow(g)) {
  g$b[i] <- bottle_xy(200, 0.1, 100, 200, g$yr[i])
}

?apply
botl <- as.list(x = g$yr)
botl <- apply(g$yr[], 1, bottle_xy, p0 = 200, cp = 0.1, tl = 100, ts = 200)
plot(g$yr, g$b, type = "l")


ggplot(data = g, aes(x = as.numeric(yr), y = as.numeric(b))) + 
  geom_line()

bottle <- lapply(g$yr, function(x) bottle_xy(200, 0.1, 100, 200, x))
plot(g$yr, g$e, type = "l")

g$e <- lapply(g$yr, function(x) expan_xy(200, 0.1, 100, 200, x))
plot(g$yr, g$e, type = "l")

itdf
plot(itdf$yr, itdf$x, type = "l")
testydf <- df

###run decline for list
model_plot <- function (mod, p0, c.prop, tl, ts, x) {
  
  
  if (mod == "decline") {
    y  <- decline_xy(p0, c.prop, tl, ts, x)
  }
  else {
    if (mod == "expan") {
      y  <- expan_xy(p0, c.prop, tl, ts, x)
    }
    else {
      if (mod == "bottle") {
        y  <- bottle_xy(p0, c.prop, tl, ts, x)
      } 
      else {
        y  <- p0;
      }
    }
  }
  return(y);
}

dd <- expand_grid(pop.init = c(200, 100),
                  crash.prop = c(0.1),
                  tl = c(100),
                  ts = c(200),
                  model = c("decline", "expan", "bottle", "stable"),
                  y = 1:250)



dd <- dd %>% 
  mutate(x = model_plot(model, pop.init, crash.prop, tl, ts, y))

####convert trajectory variables to NA for stable model


#============WARNING=============fucky code ahead
for (i in 1:nrow(test)) {
  a <- test[i,];
  b <- model_plot(a$model,
                  a$pop.init,
                  a$crash.prop,
                  a$tl,
                  a$ts,
                  a$y);
  test$x[i] <- as.numeric(b);
}


test$x <- lapply(test, model_plot(test$model, test$pop.init, test$crash.prop, test$tl, test$ts, test$y))
#==========okay its done=============

###working growth and decay models


##decay
exp.x <- data.frame(x = 0:1500, ybp = 1500 - (0:1500))
exp.x
tl <- 100
ts <- 200
y1 <- 500
yn <- 50
startyr = 1200 - ts
endyr = startyr + tl
r = log(y1 - yn)/tl

##code for forward simulation
exp.x$decay <- if_else((1200 - exp.x$ybp) < startyr,
                       y1,
                       (y1 - yn)*exp(-r * ((1200 - exp.x$ybp) - (startyr))) + yn)




exp.x$rx <- (nrow(exp.x) - exp.x$x)

plot(exp.x$decay ~ exp.x$ybp, type = "l", xlim = c(0, 1200), ylim = c(0, 1000))

bot <- data.frame(x = 0:500)


##growthlog==========displays in ybp
tl <- 20
y1 <- 100
yn <- 200
r = log(yn/y1)/tl

exp.x$growth <- if_else(exp.x$x < startyr, y1, 
                        if_else(y1*exp(r*(exp.x$x - startyr)) > yn, yn, y1*exp(r*(exp.x$x - startyr))))

exp.x$rx <- (nrow(exp.x) - exp.x$x)
plot(exp.x$growth ~ exp.x$x, type = "l", xlim = c(900, 1100), ylim = c(0, 300))



ggplot(data = test,aes(x = as.integer(y), y = x, colour = model)) +
  geom_smooth() +
  theme_minimal()

b <- data.frame(yr = seq(1, 100, 0.2))
b$p <- 100*exp(r*b$yr)
plot(b$yr, b$p, ylim = c(0,400), xlim = c(0,30), type = "l")



s <- res_epos %>% select(pop_init, crash_prop, tl, ts, model)
unique(s)

v <- crossing(yr = 1:1200, unique(s))

v
epos_simpop(v)
