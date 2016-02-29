#' ---
#' title: "Exercise 5: Introduction to JAGS"
#' author: "Alexey Shiklomanov"
#' output: html_document
#' ---

#' ### Preface: Spin instead of Rmarkdown

#' To be contrarian, for this exercise, I'm going to demonstrate the use of R spin rather than Rmarkdown. As you may recall, Spin works the opposite of Markdown, in that the text surrounding the R code has specialized formatting (simply by starting lines with the `#'` characters, while the actual R code chunks have no formatting. 

#' Chunk settings are set via the `#+` flags -- all the commands that would normally go into the `{r, ...}` structure become simply `#+ ...` (e.g. `#+ echo=FALSE`).

#' To compile the script, the command is: `rmarkdown::render('path-to-file')`. 

#' Anyway, onto the activity.

#' # Activity Task 1

library(rjags)
library(coda)

#' Running the model:
NormalMeanN <- "
model {
  mu ~ dnorm(mu0,T) # prior on the mean 
  for(i in 1:N){
    X[i] ~ dnorm(mu,S) # data model
  }
}
"

data = list(N = 297, mu0=20, T=0.01, S = 1/27, X = c(20.9, 13.6, 15.7, 6.3, 2.7, 25.6, 4, 20.9, 7.8, 27.1, 25.2, 19, 17.8, 22.8, 12.5, 21.1, 22, 22.4, 5.1, 16, 20.7, 15.7, 5.5, 18.9, 22.9, 15.5, 18.6, 19.3, 14.2, 12.3, 11.8, 26.8, 17, 5.7, 12, 19.8, 19, 23.6, 19.9, 8.4, 22, 18.1, 21.6, 17, 12.4, 2.9, 22.6, 20.8, 18.2, 14.2, 17.3, 14.5, 8.6, 9.1, 2.6, 19.8, 20, 22.2, 10.2, 12.9, 20.9, 21.1, 7.3, 5.8, 23.1, 17, 21.5, 10.1, 18.4, 22.6, 21.2, 21.5, 22.4, 17.3, 16, 25, 22.4, 23.9, 23, 21.9, 19, 28.6, 16, 22.5, 23.2, 8.7, 23.4, 15.3, 25.6, 19.2, 17.4, 23.8, 20.4, 19, 3.6, 23.4, 19.6, 17.5, 16.5, 22, 19.7, 7.35, 18, 17.8, 9.6, 15, 12, 17.7, 21.4, 17, 22.1, 18.9, 15.05, 12.9, 19.3, 15.3, 13.6, 15.4, 10.6, 11.3, 11.8, 22.2, 22.2, 13.1, 7.4, 4.5, 11.7, 19.5, 19.9, 11.6, 13.9, 15.5, 11, 18.6, 17.6, 12.7, 20.9, 18.8, 22.4, 21.2, 18.2, 15.3, 13.6, 7.3, 17.4, 17.4, 10.5, 22.9, 23.2, 13.8, 14.8, 22.2, 20.9, 13, 18.9, 19, 15.2, 16.8, 18, 24.6, 15.4, 17.2, 23.2, 22.8, 25.5, 7.8, 6, 6.4, 19, 13.5, 23.7, 18, 22.2, 22.4, 9.3, 13.7, 18.9, 20.5, 23.3, 20.8, 18.4, 4.5, 12.2, 16.9, 13.5, 17.8, 16.9, 20.4, 19.5, 22.2, 24.5, 21.2, 16.5, 18, 16.4, 3.9, 17.9, 22, 12.9, 21, 18, 9.2, 15.9, 8.1, 8.3, 10.7, 12, 19.9, 13.6, 17.3, 11.5, 12.4, 15.1, 22, 19.3, 17.5, 14.5, 14.7, 17.5, 19.6, 12.9, 20.3, 17.9, 20.2, 18.3, 9.5, 19, 21, 13.1, 20.4, 16.3, 18.3, 11.8, 23.3, 15.2, 20, 17.9, 12, 19.6, 18.5, 16.2, 10.9, 17.8, 13.8, 10, 17.9, 15.6, 20.3, 14.9, 18.6, 12.5, 18.2, 16, 18.7, 18, 15.3, 19, 17.9, 15.8, 17.7, 14.4, 19.6, 18.3, 18.7, 17.8, 18, 10.1, 18.8, 16.4, 21.2, 16.6, 16.7, 17.8, 16.5, 19.3, 16.3, 14.2, 13, 9.4, 19.7, 13.4, 2.6, 17.6, 16.7, 17.6, 5.8, 17.6, 20.1, 18.2, 16.7, 14, 13.9, 5.1, 16.6, 3.9, 17.5, 18))

nsamples <- 3000
n.chains <- 3
normal.mean.n.model <- jags.model(file = textConnection(NormalMeanN),
                                  data = data,
                                  n.chains = n.chains, quiet=TRUE)
normal.mean.n.samples <- coda.samples(model = normal.mean.n.model,
                                      variable.names = 'mu',
                                      n.iter = nsamples,
                                      progress.bar = "none")

#' The MCMC was run with `r n.chains` chains for `r nsamples` samples. Here is a plot of the Gelman diagnostic...
gelman.plot(normal.mean.n.samples)
burnin <- 1000
normal.mean.n.burnedin <- window(normal.mean.n.samples,
                                 start = burnin)

#' ...which suggests that a burnin period of `r burnin` is more than sufficient to ensure chain convergence.

#' Similarly, a plot of autocorrelation...
autocorr.plot(normal.mean.n.burnedin)

#' ...shows that even at lag 1, autocorrelation is negligible, meaning that there is no need to apply a thinning interval.

#' Combining the 3 chains and accounting for the autocorrelation (which is negligible) yields the following number of samples:

effectiveSize(normal.mean.n.burnedin)

#' ...which exceeds the 5000 sample rule of thumb for representative sampling.

#' Here is a plot of the MCMC trace before burnin...
traceplot(normal.mean.n.samples)

#' ...and after burnin...
traceplot(normal.mean.n.burnedin)

#' ...and here is a plot of the posterior density after burnin.
densplot(normal.mean.n.burnedin)

#' Here is a table of summary statistics for the posterior estimates:
summary(normal.mean.n.burnedin)

#' # Activity task 2

#' First, we modify the model code to set a prior on the precision S. I am going to use a Gamma distribution because it is conjugate with the Normal (allowing for fast Gibbs sampling), has a lower bound of zero (as is required of a precision), and its parameters (shape $\alpha$ and rate $\beta$) have meaningful interpretations in the context of precision: Namely, that the precision was estimated from $2 \alpha$ samples with a variance of $\beta \over \alpha$.

VarianceModel <- "
model {
  mu ~ dnorm(mu0,T) # prior on the mean 
  S ~ dgamma(shape, rate)  # prior on the variance
  sd <- 1/sqrt(S)   # standard deviation
  for(i in 1:N){
    X[i] ~ dnorm(mu,S) # data model
  }
}
"

#' For an uninformative prior, I will use two assumptions... 

tree.max <- 500
tree.z <- 4.75

#' ... that [the largest trees in the world](https://en.wikipedia.org/wiki/List_of_superlative_trees) have DBH of around `r tree.max` cm, and that such trees are approximately one in a million (corresponding to a z-score of around 4.75). My list of superlative trees (from Wikipedia) contains 12 trees, but my confidence in the histogram of these trees is poor, so I'll reduce that by half, giving me an effective sample size of 6 and $\alpha$ of 3. For a standard deviation of 100, the variance becomes $100^2 = 10000$, and the resulting value of $\beta$ is $10000 \over 6$:

tree.sd <- tree.max/tree.z
alpha.value <- 3
beta.value <- tree.sd^2/alpha.value

#' To test that our fitted value is reasonable, here's a histogram of standard deviations sampled from this distribution.
hist(1/sqrt(rgamma(1e5, shape=alpha.value, rate=beta.value)))

#' This gives us a reasonably informative prior on the standard deviation.

#' Now that we have the values, we assign them to the `data` list and compile the model.

data.var <- data
data.var$S <- NULL
data.var$shape <- 3
data.var$rate <- beta.value

n.chains <- 5
n.samples <- 5000

variance.model <- jags.model(file=textConnection(VarianceModel),
                             data = data.var,
                             n.chains = n.chains, quiet=TRUE)
variance.samples <- coda.samples(model=variance.model,
                                 variable.names = c("mu", "sd"),
                                 n.iter = n.samples, 
                                 progress.bar="none")

#' To determine the burnin period, we plot the Gelman diagnostic.
gelman.plot(variance.samples)
burnin <- 1000
variance.burnedin <- window(variance.samples, start=burnin)

#' Looks like a burnin interval of `r burnin` is sufficient. Next, we check the effective sample size to ensure we have sufficient samples:
effectiveSize(variance.burnedin)

#' This is well above the rule-of-thumb of 5000 samples. Finally, we produce a plot of the MCMC history and densities...
plot(variance.burnedin)

#' ...as well as the summary statistics.
summary(variance.burnedin)

#' Now, let's compare the posterior densities of the $\mu$ parameter.
mu.samples.mean <- as.matrix(normal.mean.n.burnedin)
mu.samples.var <- as.matrix(variance.burnedin[,1])
mu.dens.mean <- density(mu.samples.mean)
mu.dens.var <- density(mu.samples.var)
plot(mu.dens.mean, type='l')
lines(mu.dens.var, col=2)
legend("topright", c("const. var.", "varied var."), col=1:2, lty=1)

#' Allowing the variance to vary results in a much fatter-tailed distribution, although the peaks of the two distributions line up pretty well.

