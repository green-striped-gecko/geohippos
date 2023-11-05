####Population curve testing####

####Using logistic equation to achieve s-shaped growth or decay

#K = carrying capacity
k <- 100

#P0 = pop size at time 0
P0 

###exponential expression defines growth multiplier, ranges from 1:10
testdf <- data.frame(yr = 0:200)
testdf$tt <- 100
testdf$pop.in <- 200   
testdf$k <- 500
####rate = log(10)/(timechange timeframe /4)
testdf$r <- 4*log(10)/(testdf$tt)
testdf$pop.cur <- (testdf$pop.in*exp(testdf$r*(testdf$yr)))

testdf$exp <- exp(-testdf$r*(testdf$yr - 0.5*testdf$tt))


####For k < pop.in --- decline
testdf$log <- ((testdf$pop.in - testdf$k)/(1 + (((testdf$pop.in/testdf$k)-1)*exp(testdf$r*(testdf$yr - 0.5*testdf$tt))))) + testdf$k
plot(testdf$yr, testdf$log, type = "l")
testdf

testdf <- testdf %>%
  mutate("log" = (pop.in - k))

####For k < pop.in --- testexp
testdf$log <- (-(testdf$pop.in - testdf$k)/(1 + ((((testdf$pop.in/testdf$k)^-1)-1)*exp(-testdf$r*(testdf$yr - 0.5*testdf$tt))))) + testdf$k
plot(testdf$yr, testdf$log, type = "l")
testdf

####For k > pop.in --- expansion
testdf$log <- ((testdf$k - testdf$pop.in)/(1 + (((testdf$k/testdf$pop.in)-1)*exp(-testdf$r*((testdf$yr) - 0.5*testdf$tt))))) + testdf$pop.in
plot(testdf$yr, testdf$log, type = "l")
testdf



plotdf <- df

###generate function for decline model
decline_y_to_x <- function(p0, crash.prop, tlen, tstart, y) {
  if (y > (tstart - 1)) {
    x <- p0;
  }
  tend <- tstart - tlen;
  k <- p0 * crash.prop;
  r <- 4*log(10)/(tlen);
  if (y < tstart & y >= tend ) {
    x <- ((p0 - k)/(1 + (((p0/k)-1)*exp(r*((tstart - y) - 0.5*tlen))))) + k;
  }
  if (y < tend) {
    x <- k;
  } 
  
  
  return(x)
}

###generate expansion model
expan_y_to_x <- function(p0, crash.prop, tlen, tstart, y) {
  k <- p0;
  p0 <- k*crash.prop;
  if (y > (tstart - 1)) {
    x <- p0;
  }
  tend <- tstart - tlen;
  r <- 4*log(10)/(tlen);
  if (y < tstart & y >= tend ) {
    x <- ((k - p0)/(1 + ((k/p0) - 1)*exp(-r*((tstart - y) - 0.5*tlen)))) + p0;
  }
  return(x)
}


###generate bottleneck model
bottle_y_to_x <- function(p0 = 200, crash.prop = 0.5, tlen = 50, tstart = 100, y) {
  k <- p0 * crash.prop;
  tend <- tstart - tlen;
  r <- 4*log(10)/(tlen);
  if (y < (tstart - 1)) {
    x = ((p0 - k)/(1 + ((p0/k) - 1)*exp(-r*((tstart - y) - 0.5*tlen)))) + k;
  } else{
     if (y == tstart) {
      x = k
     } else { x <- p0}
  }
  return(x)
}

bottle_y_to_x(y = 150)


testydf <- df

###run decline for list
model_plot <- function (model, p0, c.prop, tl, ts, y) {
  if (model == "decline") {
    x <- decline_y_to_x(p0, c.prop, tl, ts, y)
  }
  else {
    if (model == "expan") {
      x <- expan_y_to_x(p0, c.prop, tl, ts, y)
    }
    else {
     if (model == "bottle") {
      x <- bottle_y_to_x(p0, c.prop, tl, ts, y)
     } 
      else {
          x <- p0;
      }
    }
  }
  return(x);
}


  
testydf
  
for (i in 1:nrow(testydf)) {
    testydf$x[i] <- bottle_y_to_x(200, 0.5, 50, 100, testydf$y[i]);
}

test <- expand.grid(pop.init = c(500, 200, 100, 50),
                    crash.prop = c(0.5, 0.1),
                    tl = c(100, 50, 25),
                    ts = c(200, 100, 50),
                    model = c("decline", "expan", "bottle", "stable"),
                    y = 1:250)
test

####convert trajectory variables to NA for stable model



for (i in 1:nrow(test)) {
  test$x[i] <- model_plot(test$model[i],
                          test$pop.init[i],
                          test$crash.prop[i],
                          test$tl[i],
                          test$ts[i],
                          test$y[i])
}

test <- test %>%
  group_by(model, pop.init, crash.prop, tl, ts) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()


stable <- which(test$model=="stable")

test[stable,]$crash.prop <- NA
test[stable,]$tl <- NA
test[stable,]$ts <- NA

test <- unique(test) 

length(test[1,]$y)
nrow(test)



ggplot(data = test
       , aes(x = y, y = x, colour = as.factor(pop.init), group = group_id)) +
  geom_line() +
  ylab("Effective population size (Ne)") +
  xlab("Time (years before present)") +
  facet_wrap(~model)


exp.x <- data.frame(x = 0:300)
exp.x
tl <- 50
y1 <- 500
yn <- 100
r = log(y1 - yn)/tl


exp.x$decay <- (y1 - yn)*exp(-r * exp.x$x) + yn
plot(exp.x$decay ~ exp.x$x, type = "l", xlim = c(0, 300), ylim = c(0, 1000))


tl <- 100
y1 <- 50
yn <- 200
r = log(yn/y1)/tl

exp.x$growth <- y1*exp(r*exp.x$x)
plot(exp.x$growth ~ exp.x$x, type = "l", xlim = c(0, 300), ylim = c(0, 1000))



