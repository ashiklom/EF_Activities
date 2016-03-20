#' ---
#' title: "Exercise 7: Fusing Time-series data"
#' author: "Alexey Shiklomanov"
#' ---

#' Load packages
#+ load-pkgs, results='hide'
library(PEcAn.data.land)
library(rjags)

#' Prepare confidence interval function
#+ cifunc
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

#' Read and organize tree data
#+read-data, results='hide', message=FALSE, warning=FALSE
trees <- read.csv("../data/H2012AdultFieldData.csv")
rings <- Read_Tucson("../data/TUCSON/")
combined <- matchInventoryRings(trees,rings,nyears=15)
data <- buildJAGSdata_InventoryRings(combined)

#' Define model and plot function
#+ def-runmodel-function
run.model <- function(modelText, n.iter = 20000, n.chain = 3){

    ## state variable initial condition
    z0 <- t(apply(data$y,1,function(y){-rev(cumsum(rev(y)))})) + data$z[,ncol(data$z)] 
    
    ## JAGS initial conditions
    init <- list()
    for(i in 1:n.chain){
        y.samp = sample(data$y,length(data$y),replace=TRUE)
        init[[i]] <- list(x = z0,tau_add=runif(1,1,5)/var(diff(y.samp),na.rm=TRUE),
                        tau_dbh=1,tau_inc=500,tau_ind=50,tau_yr=100,ind=rep(0,data$ni),year=rep(0,data$nt))
    }

    ## compile JAGS model
    j.model   <- jags.model (file = textConnection(modelText),
                            data = data,
                            inits = init,
                            n.chains = n.chain)
    ## burn-in
    jags.out   <- coda.samples (model = j.model,
                                variable.names = c("tau_add","tau_dbh","tau_inc","mu","tau_ind","tau_yr"),
                                n.iter = min(n.iter,2000))
    ## run MCMC
    jags.out   <- coda.samples (model = j.model,
                                variable.names = c("x","tau_add","tau_dbh","tau_inc","mu",
                                                    "tau_ind","tau_yr","ind","year"),
                                n.iter = n.iter)
    out <- as.matrix(jags.out)
    x.cols = which(substr(colnames(out),1,1)=="x")   ## which columns are the state variable, x
    ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))
    ci.names = parse.MatrixNames(colnames(ci),numeric=TRUE)
    return(list(out = out, x.cols = x.cols, ci = ci, ci.names = ci.names))
}

#' Run both models

#+ run-norf, results='hide', warning=FALSE
n.iter <- 20000
TreeDataFusionMV = "
model{

  ### Loop over all individuals
  for(i in 1:ni){
  
    #### Data Model: DBH
    for(t in 1:nt){
        z[i,t] ~ dnorm(x[i,t],tau_dbh)
    }
    
    #### Data Model: growth
    for(t in 2:nt){
        inc[i,t] <- x[i,t]-x[i,t-1]
        y[i,t] ~ dnorm(inc[i,t],tau_inc)
    }
    
    #### Process Model
    for(t in 2:nt){
        Dnew[i,t] <- x[i,t-1] + mu
        x[i,t]~dnorm(Dnew[i,t],tau_add)
    }
    
    x[i,1] ~ dnorm(x_ic,tau_ic)
  }  ## end loop over individuals
  
  #### Priors
  tau_dbh ~ dgamma(a_dbh,r_dbh)
  tau_inc ~ dgamma(a_inc,r_inc)
  tau_add ~ dgamma(a_add,r_add)
  mu ~ dnorm(0.5,0.5)
}"

norf.result <- run.model(TreeDataFusionMV, n.iter=n.iter)

#+ run-rf, results='hide', warning=FALSE
TreeDataFusionMVRF = "
model{
  
  ### Loop over all individuals
  for(i in 1:ni){
  
  #### Data Model: DBH
  for(t in 1:nt){
    z[i,t] ~ dnorm(x[i,t],tau_dbh)
  }
  
  #### Data Model: growth
  for(t in 2:nt){
    inc[i,t] <- x[i,t]-x[i,t-1]
    y[i,t] ~ dnorm(inc[i,t],tau_inc)
  }
  
  #### Process Model
  for(t in 2:nt){
    Dnew[i,t] <- x[i,t-1] + mu + ind[i] + year[t]
    x[i,t]~dnorm(Dnew[i,t],tau_add)
  }
  
  ## individual effects
  ind[i] ~ dnorm(0,tau_ind)
  
  ## initial condition
  x[i,1] ~ dnorm(x_ic,tau_ic)
  }  ## end loop over individuals
  
  ## year effects
  for(t in 1:nt){
    year[t] ~ dnorm(0,tau_yr)
  }
  
  
  #### Priors
  tau_dbh ~ dgamma(a_dbh,r_dbh)
  tau_inc ~ dgamma(a_inc,r_inc)
  tau_add ~ dgamma(a_add,r_add)
  tau_ind ~ dgamma(1,0.1)
  tau_yr  ~ dgamma(1,0.1)
  mu ~ dnorm(0.5,0.5)
  
}"

rf.result <- run.model(TreeDataFusionMVRF, n.iter=n.iter)

#' ## Ouput plots and comparison

#+ plot-dbh-growth, fig.width=8, fig.height=8
smp <- c(sample.int(data$ni,3),49)
par(mfcol=c(4,2), mar=c(4,4,1,1))
norf.col <- rgb(0,0,1,0.3)
rf.col <- rgb(1,0,0,0.3)

### DBH
for(i in smp){
    sel = which(rf.result$ci.names$row == i)
    yrange = range(c(rf.result$ci[,sel], norf.result$ci[,sel], data$z[i,]),
                   na.rm=TRUE)
    plot(data$time, data$time, type='n',
         ylim=yrange,
         ylab="DBH (cm)", 
         main=i)
    ciEnvelope(data$time, norf.result$ci[1,sel], norf.result$ci[3,sel],
               col=norf.col)
    ciEnvelope(data$time, rf.result$ci[1,sel], rf.result$ci[3,sel],
               col=rf.col)
    points(data$time,data$z[i,],pch="+",cex=1.5)
    if(i == smp[1]) legend("bottomright", legend = c("no RE", "RE"), lty=1, col=c(norf.col, rf.col))
}

### Growth
for(i in smp){
    sel = which(rf.result$ci.names$row == i)
    inc.mcmc.norf = apply(norf.result$out[,norf.result$x.cols[sel]],1,diff)
    inc.ci.norf = apply(inc.mcmc.norf,1,quantile,c(0.025,0.5,0.975))*5
    inc.mcmc.rf = apply(rf.result$out[,rf.result$x.cols[sel]],1,diff)
    inc.ci.rf = apply(inc.mcmc.rf,1,quantile,c(0.025,0.5,0.975))*5
    yrange <- range(c(inc.ci.rf, inc.ci.norf, data$y[i,]*5))
    plot(data$time[-1], data$time[-1], type='n', ylim=yrange,
         ylab="Ring Increment (mm)")
    ciEnvelope(data$time[-1], inc.ci.norf[1,], inc.ci.norf[3,], col=norf.col)
    ciEnvelope(data$time[-1], inc.ci.rf[1,], inc.ci.rf[3,], col=rf.col)
    points(data$time,data$y[i,]*5,pch="+",cex=1.5,type='b',lty=2)
}
dev.off()

#' The random effects model performed better than the model without random effects. First of all, especially for the growth trees, the confidence intervals around individual tree growth trajectories are tighter. Second, without random effects, the model significantly underestimated the growth of tree 49 for most of the record, while the model with random effects was able to capture growth in every year except the most extreme peak in 2009.

#+ process-model-plot
hist(rf.result$out[,"mu"], main="Mu", col=rf.col, xlab="value")
hist(norf.result$out[,"mu"], col=norf.col, add=TRUE)
legend("topleft", c("No RE", "RE"), pch=19, pt.cex=2, col=c(norf.col, rf.col))

## Standard deviations
taus <- colnames(rf.result$out)[grep("tau", colnames(rf.result$out))]
taus.norf <- colnames(norf.result$out)[grep("tau", colnames(norf.result$out))]
prec.norf <- 1 / sqrt(norf.result$out[,taus.norf])
prec.rf <- 1 / sqrt(rf.result$out[,taus])
prec.to.sd <- function(x) 1/sqrt(x)
par(mfrow=c(2,3))
for(tau in taus){
    hastau <- tau %in% colnames(prec.norf)
    rf.hist <- hist(prec.rf[,tau], plot=FALSE)
    if(hastau){
        norf.hist <- hist(prec.norf[,tau], plot=FALSE)
    } else {
        norf.hist <- NULL
    }
    xrange <- range(c(rf.hist$breaks, norf.hist$breaks), na.rm=TRUE)
    yrange <- range(c(rf.hist$counts, norf.hist$counts), na.rm=TRUE)
    plot(rf.hist, col=rf.col, main=tau, xlim=xrange, ylim=yrange, xlab="SD")
    if(hastau) plot(norf.hist, col=norf.col, add=TRUE)
}

cor(prec.norf)
cor(prec.rf)

pairs(prec.norf, pch=".")
pairs(prec.rf, pch=".")

#' The mean for the model without random effects is much narrower than the mean for the model with random effects. This is because in the former, the mean has a lot more degrees of freedom (i.e. sample size) that are not being used to estimate individual parameters for every individual tree and year. In other words, if we assume (as that model does) that every tree shares a single mean and individuals vary only because of uncertainty, then collecting more trees should increase out confidence in the value of that overall mean. However, in that case, we see that the standard deviation about the mean (`tau_add`, though the same is true of `tau_inc`) is also larger, which makes sense since all variability is assumed to be from random variability about the overall mean.
#'
#' By contrast, in the random effects model, while we are less certain about the precise value of the mean, we are able to explain some of the variance around it via differences in individuals and years, thereby reducing the value of the standard deviations. By partitioning variability, we are also able to reduce the covariance in the estimates of the standard deviations. When we don't distinguish between individuals and years, our estimate of the standard deviation is understandably confused between uncertainty in the measurement (`tau_inc`) and variance in the process (`tau_add`), hence the negative covariance between them. However, because measurement error is more closely tied to individual trees while process error is more likely temporal in nature, accounting for individuals and years helps decouple these two uncertainties (in addition to reducing them).

out <- rf.result$out
par(mfrow=c(1,1))
### YEAR
year.cols = grep("year",colnames(out))
if(length(year.cols>0)){
    ci.yr <- apply(out[,year.cols],2,quantile,c(0.025,0.5,0.975))
    plot(data$time,ci.yr[2,],type='n',ylim=range(ci.yr,na.rm=TRUE),main="Year Effect",ylab="cm")
    ciEnvelope(data$time,ci.yr[1,],ci.yr[3,],col="lightBlue")
    lines(data$time,ci.yr[2,],lty=1,lwd=2)
    abline(h=0,lty=2)
}

### INDIV
ind.cols= which(substr(colnames(out),1,3)=="ind")
boxplot(out[,ind.cols],horizontal=TRUE,outline=FALSE,col=combined$PLOT,main="Individual Effects By Plot",xlab="cm")
abline(v=0,lty=2)

## calculate plot-level means for random effects
tapply(apply(out[,ind.cols],2,mean),combined$PLOT,mean)
table(combined$PLOT)

spp = combined$SPP
boxplot(out[order(spp),ind.cols],horizontal=TRUE,outline=FALSE,col=spp[order(spp)],main="Individual Effects By Species",xlab="cm")
abline(v=0,lty=2)
spp.code = levels(spp)[table(spp)>0]
legend("bottomright",legend=rev(spp.code),col=rev(which(table(spp)>0)),lwd=4)

## calculate species-level means for random effects
spp.means <- tapply(apply(out[,ind.cols],2,mean, na.rm=TRUE),combined$SPP,mean, na.rm=TRUE)
print(spp.means)
barplot(spp.means)

#' To me, the most apparent target for adding to the model is to add a fixed species effect. For one, this effect has larger mean random effect values than plot, but more importantly, this effect captures physiological differences in plant growth rates based on their differing responses to environmental conditions in a competitive setting.
#'
#' This effect would require minimal additional exploratory analysis -- basically, just to examine the species means, as we have done above. From this analysis, it's clear that several specific species emerge as having strong differences in growth rates from the overall mean -- specifically, *Acer rubrum* and *Betula lenta* both have large negative means, while *Pinus strobus* and *Betula spp.* have large positive means.
#'
#' The JAGS code for this analysis would look like this:
mod <- "
model{
  ### Loop over all individuals
  for(i in 1:ni){
  
  #### Data Model: DBH
  for(t in 1:nt){
    z[i,t] ~ dnorm(x[i,t],tau_dbh)
  }
  
  #### Data Model: growth
  for(t in 2:nt){
    inc[i,t] <- x[i,t]-x[i,t-1]
    y[i,t] ~ dnorm(inc[i,t],tau_inc)
  }
  
  #### Process Mode1:l
  for(t in 2:nt){
    Dnew[i,t] <- x[i,t-1] + mu + ind[i] + year[t] + sp[sp_list[i]]
    x[i,t]~dnorm(Dnew[i,t],tau_add)
  }
  
  ## individual effects
  ind[i] ~ dnorm(0,tau_ind)
  
  ## initial condition
  x[i,1] ~ dnorm(x_ic,tau_ic)
  }  ## end loop over individuals
  
  ## year effects
  for(t in 1:nt){
    year[t] ~ dnorm(0,tau_yr)
  }
  
  
  #### Priors
  tau_dbh ~ dgamma(a_dbh,r_dbh)
  tau_inc ~ dgamma(a_inc,r_inc)
  tau_add ~ dgamma(a_add,r_add)
  tau_ind ~ dgamma(1,0.1)
  tau_yr  ~ dgamma(1,0.1)
  mu ~ dnorm(0.5,0.5)
  for(s in 1:nspecies) sp[s] ~ dnorm(0, 0.5)
}"
