#' ---
#' title: "Exercise 9: Kalman Filter"
#' author: "Alexey Shiklomanov"
#' output: html_document
#' ---

#' # Prepare data

## load the Google flu data & select states
gflu = read.csv("http://www.google.org/flutrends/about/data/flu/us/data.txt",skip=11)
time = as.Date(gflu$Date)
states = c("Massachusetts","Connecticut","Rhode.Island","New.Hampshire","Vermont","Maine")
nstates = length(states)
y = t(gflu[,states])

## define adjacency between states slected
adj = matrix(c(0,1,1,1,1,0,    ### state-to-state spatial adjacency (self=0)
               1,0,1,0,0,0,
               1,1,0,0,0,0,
               1,0,0,0,1,1,
               1,0,0,1,0,0,
               0,0,0,1,0,0),nstates,nstates,byrow=TRUE)

## plot time-series from states
plot(time,1:length(time),type='n',ylab="Flu Index",lwd=2,log='y',ylim=range(y,na.rm=TRUE))
for(i in 1:nstates){
  lines(time,y[i,],col=i,lwd=2)
}
legend("topleft",legend=states,lwd=2,col=1:nstates)

#' # Define Kalman Filter function

##'  Kalman Filter
##' @param  M   = model matrix
##' @param  mu0 = initial condition mean vector
##' @param  P0  = initial condition covariance matrix
##' @param  Q   = process error covariance matrix
##' @param  R   = observation error covariance matrix
##' @param  Y   = observation matrix (with missing values as NAs), time as col's
##'
##' @return list
##'  mu.f, mu.a  = state mean vector for (a)nalysis and (f)orecast steps
##'  P.f, P.a    = state covariance matrix for a and f
KalmanFilter <- function(M,mu0,P0,Q,R,Y){
  
  ## storage
  nstates = nrow(Y)  
  nt = ncol(Y)
  mu.f  = matrix(NA,nstates,nt+1)  ## forecast mean for time t
  mu.a  = matrix(NA,nstates,nt)  ## analysis mean for time t
  P.f  = array(NA,c(nstates,nstates,nt+1))  ## forecast variance for time t
  P.a  = array(NA,c(nstates,nstates,nt))  ## analysis variance for time t

  ## initialization
  mu.f[,1] = mu0
  P.f[,,1] = P0
  I = diag(1,nstates)

  ## run updates sequentially for each observation.
  for(t in 1:nt){

    ## Analysis step: combine previous forecast with observed data
    obs = !is.na(Y[,t]) ## which Y's were observed?
    if(any(obs)){
      H <- I[obs,]                                                        ## observation matrix
      K <- P.f[,,t] %*% t(H) %*% solve(H%*%P.f[,,t]%*%t(H) + R[obs,obs])  ## Kalman gain
      mu.a[,t] <- mu.f[,t] + K%*%(Y[obs,t] - H %*% mu.f[,t])              ## update mean
      P.a[,,t] <- (1-K %*% H)*P.f[,,t]                                    ## update covariance
    } else {
      ##if there's no data, the posterior is the prior
      mu.a[,t] = mu.f[,t]
      P.a[,,t] = P.f[,,t]
    }

    ## Forecast step: predict to next step from current
    mu.f[,t+1] = M%*%mu.a[,t]
    P.f[,,t+1] = Q + M*P.a[,,t]*t(M)
  
  }
  
  return(list(mu.f=mu.f,mu.a=mu.a,P.f=P.f,P.a=P.a))
}

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

## log transform data
Y   = log10(y)

## load parameters (assume known)
load("data/KFalpha.params.Rdata")



#' # Set up analysis function

kf.analysis <- function(alpha, Q){

    ## options for process model 
    M = adj*alpha + diag(1-alpha*apply(adj,1,sum))  ## random walk with flux

    ## observation error covariance (assumed independent)  
    R = diag(tau_obs,nstates) 

    ## prior on first step, initialize with long-term mean and covariance
    mu0 = apply(Y,1,mean,na.rm=TRUE)
    P0 = cov(t(Y),use="pairwise.complete.obs")

    ## Run Kalman Filter
    KF = KalmanFilter(M,mu0,P0,Q,R,Y)

    return(KF)
}

#' # Run analysis

# Parameters for process model
alpha.simple = 0       ## assume no spatial flux
alpha.spatial = 0.05    ## assume a large spatial flux

# Parameters for process error covariance
Q.full = tau_proc            ## full covariance matrix
Q.diag = diag(diag(Q.full))       ## diagonal covariance matrix

KF.simple.full <- kf.analysis(alpha.simple, Q.full)
KF.spatial.full <- kf.analysis(alpha.spatial, Q.full)
KF.simple.diag <- kf.analysis(alpha.simple, Q.diag)
KF.spatial.diag <- kf.analysis(alpha.spatial, Q.diag)

#' # Define results function

kf.results <- function(KF){
    attach(KF)
    nt = length(time)

    ### plot ANALYSIS mean & CI time-series
    par(mfrow=c(3,1))
    for(i in 1:6){
        ci = rbind(mu.a[i,]-1.96*sqrt(P.a[i,i,]),mu.a[i,]+1.96*sqrt(P.a[i,i,]))
        plot(time,mu.a[i,],ylim=range(ci,na.rm=TRUE),type='n',main=states[i])
        ciEnvelope(time,ci[1,],ci[2,],col="lightBlue")
        lines(time,mu.a[i,],col=4)
        lines(time,Y[i,])
    }

    ## plot ANALYSIS and FORECAST variance time-series
    par(mfrow=c(3,1))
    for(i in 1:6){
        plot(time,sqrt(P.a[i,i,]),ylim=c(0,sqrt(max(c(P.a[i,i,],P.f[i,i,])))),main=states[i],xlab="Time",
             ylab="Std Error",type='l')
        lines(time,sqrt(P.f[i,i,1:nt]),col=2)
        points(time[is.na(Y[i,])],rep(0,nt)[is.na(Y[i,])],pch="*",col=3) ## flag's the zero's
        legend("topright",legend=c("Analysis","Forecast","NAs"),col=1:3,lty=c(1,1,NA),pch=c(NA,NA,1),cex=1.4)
    }
    detach(KF)
}

#' # Analyze results
#'
#' ## No spatial flux, diagonal covariance
#'
#+ fig.height=12
kf.results(KF.simple.diag)

#' This is the simplest model, and, predictably, is the least informative, failing to capture any of the predictable seasonality in flu trends when data are not available and having extremely high uncertainty about the magnitude of flu values.
#'
#' ## No spatial flux, full covariance
#'
#+ fig.height=12
kf.results(KF.simple.full)

#' The inclusion of the full covariance matrix means allows estimates for states with missing data to at least capture the shape of flu trends, which are predictably higher during winter than summer months. However, in the absence of spatial constraint, the certainty in the magnitudes of flu values is very low.
#'
#' ## Spatial flux, diagonal covariance
#'
#+ fig.height=12
kf.results(KF.spatial.diag)

#' The inclusion of spatial constraint dramatically reduces model uncertainty where data are missing, even moreso than the inclusion of the full covariance matrix. Even without the full covariance, this model predicts flu seasonality. However, noted in the assignment, the actual model estimate of $\alpha$ is much lower than that used here, meaning that this model overestimates the role of spatial proximity in flu predictability. Also, Maine--as a state with only one neighbor--had the least to gain from this approach in terms of uncertainty reduction, while other states with more neighbors fared better.
#'
#' ## Spatial flux, full covariance
#'
#+ fig.height=12
kf.results(KF.spatial.full)
#'
#' The inclusion of the full covariance to a model that already used spatial proximity did little to reduce the uncertainty around flu estimates where data are missing, especially in states that had multiple neighbors. Although the model did not much reduce the uncertainty for flu in Maine during the missing data period, it did enhance the "peakiness" of the flu trend during this time (compared to the diagonal covariance), reflecting the model's prediction that large, sudden flu outbreaks happen in all places simultaneously rather than diffusing gradually.
#'
#' ## Nonlinear process model
#'
#' The Extended Kalman Filter (EKF) uses a first-order Taylor series approximation to account for non-linearity in the model. Since the first-order Taylor approximation is based on the first derivative of the function at hand, performing the EKF on this model would require analytically solving for the full Jacobian matrix (derivative with respect to every parameter) of this model. This Jacobian derivative matrix (which changes numerical values at every time step) would replace the currently fixed slope matrix $M$.
