#' ---
#' title: "Exercise 6: State-space models"
#' author: "Alexey Shiklomanov"
#' ---

#' # Main activity

#' First, I load the flu data.

library(rjags)
set.seed(7111992)

gflu <- read.csv("http://www.google.org/flutrends/about/data/flu/us/data.txt",skip=11)
time <- as.Date(gflu$Date)
y <- gflu$Massachusetts

#' Then, I randomly remove 3 out of every 4 observations, and move the removed 
#' observations into a new vector for storage.

obs.to.remove <- sample(seq_along(y), round(length(y)*3/4))
y.subsamp <- y
y.removed <- y.subsamp
y.removed[-obs.to.remove] <- NA
y.subsamp[obs.to.remove] <- NA

#' Next, I setup the model code.

RandomWalk <- "
model{
  
  #### Data Model
  for(i in 1:n){
    y[i] ~ dnorm(x[i],tau_obs)
  }
  
  #### Process Model
  for(i in 2:n){
    x[i]~dnorm(x[i-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

#' Next, I write a function to call the same model with a particular input. 
#' This will allow me to quickly run both the original model (with the full 
#' data) and the model with a subset of the data removed.

run.model <- function(y){
    data <- list(y=log(y), n=length(y),x_ic=log(1000),
                 tau_ic=100,a_obs=1,r_obs=1,a_add=1,r_add=1)

    nchain <- 3
    init <- list()
    for(i in 1:nchain){
        y.samp <- sample(y,length(y),replace=TRUE)
        init[[i]] <- list(tau_add=1/var(diff(log(y.samp)), na.rm=TRUE),
                          tau_obs=5/var(log(y.samp), na.rm=TRUE))
    }

    j.model   <- jags.model (file = textConnection(RandomWalk),
                             data = data,
                             inits = init,
                             n.chains = 3)

    jags.out   <- coda.samples (model = j.model,
                                variable.names = c("x","tau_add","tau_obs"),
                                n.iter = 10000,
                                progress.bar = 'none')
    out <- as.matrix(jags.out)
    return(out)
}

#' Similarly, I write a function that will plot the time series with confidence 
#' intervals for each model.

plot.model <- function(out, y.in, y.out=NULL, zoom=c(1, length(time))){
    time.rng = zoom
    ciEnvelope <- function(x,ylo,yhi,...){
        polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                            ylo[1])), border = NA,...) 
    }
    ci <- apply(exp(out[,3:ncol(out)]),2,quantile,c(0.025,0.5,0.975))

    plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Flu Index",log='y',xlim=time[time.rng])
    ## adjust x-axis label to be monthly if zoomed
    if(diff(time.rng) < 100){ 
        axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
    }
    ciEnvelope(time,ci[1,],ci[3,],col="lightBlue")
    points(time,y.in,pch="+",cex=0.5)
    if(!is.null(y.out)) points(time, y.out, pch="o", cex=0.5, col="red")
}

#' With the functions in hand, I run them once for each model and plot the outputs in a single figure.

out.full <- run.model(y)
out.rm <- run.model(y.subsamp)

par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.model(out.full, y)
plot.model(out.rm, y.subsamp, y.removed)

#' The 95\% confidence intervals around the model with the full data are much 
#' narrower than those around the model with three quarters of the data 
#' removed. Furthermore, in the latter case, some of the data points actually 
#' fall outside of the 95\% confidence interval.
#'
#' Next, I look at how well our random walk model was able to predict the data 
#' points that I artificially removed.

y.predicted <- apply(exp(out.rm[,obs.to.remove+2]), 2, median)
plot(y.predicted, y.removed[obs.to.remove],
     xlab='Predicted', ylab='Observed', pch='+',
     log='xy')
abline(a=0, b=1)

#' It looks like our model was pretty good at predicting the flu at low values, 
#' but consistently underestimated high flu peaks.
#'
#' Based on the combined time series and 1:1 plots, model accuracy was 
#' primarily driven by where data were missing -- i.e. when data were missing 
#' at around the peak of the flu season, the model was more likely to 
#' underestimate the flu index -- while model precision was primarily driven by 
#' the length of the period of missing data -- i.e. longer stretches with 
#' missing data resulted in much wider confidence intervals than shorter 
#' stretches.
#'
#' # Extra credit
#'
#' Here, I remove the last 40 observations from `y`.

y.rm.tail <- y
rm.tail <- length(y.rm.tail) - (40:0)
y.rm.tail[rm.tail] <- NA
y.keep.tail <- y
y.keep.tail[-rm.tail] <- NA

#' Next, I use my existing functions to refit the model as a pseudo-forecast.

out.forecast <- run.model(y.rm.tail)
plot.model(out.forecast, y.rm.tail, y.keep.tail, zoom=c(length(time)-80, length(time)))

#' Although the random walk did capture almost all of the observations in its 
#' confidence interval, it did so by naively widening it's confidence interval 
#' as a funciton of time, giving the forecast neither accuracy (the median is 
#' the value of the last observation) nor precision (the uncertainty quickly 
#' grows to larger than the observed range of values). A first step to 
#' improving this forecast would be to recognize the inherent seasonality of 
#' flu outbreaks, which could be accomplished by introducing parameters for 
#' outbreak timing (fixed to a typical month of outbreak, or linked to climate 
#' and social variables), for the background flu rate reached outside of the 
#' outbreak seasons (seems like it drops to a few hundred cases), and for the 
#' maximum flu index (consistently between 1000 and 5000).

