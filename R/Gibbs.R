##' Gibbs sampling for simple linear regression model with one parameter and one outcome
##'
##' Gibbs sampling or a Gibbs sampler is a Markov chain Monte Carlo (MCMC) algorithm for obtaining a 
##' sequence of observations which are approximated from a specified multivariate probability 
##' distribution. Gibbs sampling is commonly used as a means of statistical inference, especially 
##' Bayesian inference. From which we can get the approximate joint posterior distribution and marginal 
##' posterior distribution.
##' @title Gibbs sampling for simple linear regression model
##' @param X Input values used as independent variable in the linear model
##' @param Y Observed values used as dependent variable in the linear model
##' @param n Number of input values/number of rows in the data set used to fit the linear model
##' @return GibbsSampleResult: Sampling results for Intercept, Slope and sigma2 (sampling for 30000 times)
##' @return GibbsSampleSummary: Summary of the sampling results (Mean, SD, 2.5% quantile, 97.5% quantile)
##' @return SD_Comparison: Comparison between the standard deviation of the posterior samples and the standard deviation of OLS
##' @return Plot of X and Y with fitted line from OLS
##' @return Histogram of the posterior slope
##' @return Plot of X and Y with the fitted line of Intercept mean and Slope mean from sampling
##' @author Luer Zhong
##' @export

gibbs <- function(X, Y, n){
  OLS <- lm(Y~X)
  summary(OLS)
  plot(X,Y,xlab="Variate",ylab="Response", main = 'Plot with fitted line from OLS')
  print(paste0('The intercept for linear model is ', round(OLS$coef[1],3)))
  print(paste0('The slope for linear model is ', round(OLS$coef[2],3)))
  y_hat <- OLS$coef[1]+OLS$coef[2]*X
  lines(X,y_hat)
  
  ### Priors
  mu0 <- 0
  s20 <- 1000
  a   <- 0.01
  b   <- 0.01
  
  n.iters <- 30000
  keepers <- matrix(0,n.iters,3)
  colnames(keepers)<-c("alpha","beta","sigma2")
  
  # Initial values
  alpha       <- OLS$coef[1]
  beta        <- OLS$coef[2]
  s2          <- var(OLS$residuals)
  keepers[1,] <- c(alpha,beta,s2)
  
  for(iter in 2:n.iters){
    # sample alpha
    V     <- n/s2+mu0/s20
    M     <- sum(Y-X*beta)/s2+1/s20
    alpha <- rnorm(1,M/V,1/sqrt(V))
    # sample beta
    V     <- sum(X^2)/s2+mu0/s20
    M     <- sum(X*(Y-alpha))/s2+1/s20
    beta  <- rnorm(1,M/V,1/sqrt(V))
    # sample s2|mu,Y,Z
    A  <- n/2 + a
    B  <- sum((Y-alpha-X*beta)^2)/2 + b
    s2 <- 1/rgamma(1,A,B)
    # keep track of the results
    keepers[iter,] <- c(alpha,beta,s2)
  }
  colnames(keepers) <- c('Intercept', 'Slope', 'sigma2')
  GibbsSampleResult <<- data.frame(keepers)
  print(head(GibbsSampleResult))
  
  output <- matrix(0,3,4)
  rownames(output) <- c("Intercept","Slope","sigma2")
  colnames(output) <- c("Mean","SD","Q025","Q975")
  
  output[,1] <- apply(keepers,2,mean)
  output[,2] <- apply(keepers,2,sd)
  output[,3] <- apply(keepers,2,quantile,0.025)
  output[,4] <- apply(keepers,2,quantile,0.975)
  #Summarize the sampling results
  print(output)
  GibbsSampleSummary <<-output
  
  beta <- keepers[,2]
  hist(beta,main="Posterior of the slope, beta",breaks=100)
  
  fit_bayes <- output[1:2,1]
  plot(X,Y,xlab="Variate",ylab="Response", main = 'Plot with the fitted line from sampling')
  lines(X,fit_bayes[1]+fit_bayes[2]*X)
  
  posterior_sd <- output[1:2,2]
  OLS_sd <- sqrt(diag(vcov(OLS)))
  #Comparing the standard deviation of the posterior samples to the standard error of the regression
  compare <- cbind(posterior_sd, OLS_sd)
  print(compare)
  SD_Comparison <<- compare
}



