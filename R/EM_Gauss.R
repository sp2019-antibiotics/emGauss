##################
# Packages ######
#################
#library(invgamma)

########################
## Functions ###########
########################

#Probabilty of being in interval j under component k
pjk <- function(a, b, mu, sigma2){
  # a numeric value lower bound of a bin
  # b numeric value upper bound of a bin
  # mu numeric value mean of a normal
  # sigma2 numeric value variance of a normal
  #OUTPUT
  # number - Area unter the density

    pjk_est <- stats::pnorm(b, mu, sqrt(sigma2)) -
      stats::pnorm(a, mu, sqrt(sigma2))

    pjk_est <- ifelse(pjk_est == 0,
                      .Machine$double.eps,
                      pjk_est)


    return(pjk_est)
}

# Probabilty of being in interjal j under mixed distribution
pj <- function(p0, pi, a, b, mu , sigma2, get.p0= FALSE){
  # p0 - numeric value probability of being in interval B0
  # a numeric value lower bound of a bin
  # b numeric value upper bound of a bin
  # pi numeric vector
  # mu numeric vector mean of a normal
  # sigma2 numeric vector variance of a normal

  buffer <- NULL
  for (i in length(mu)){
    buffer <- c(buffer,
                pi[i]*pjk(a, b, mu[i], sigma2[i]))


  }

  if(get.p0){
    return(sum(buffer))
  }else{
    correction.term <- 1/(1-p0)
    return(correction.term * sum(buffer))
  }

}


#Loglikelihood
loglik <- function(n0, p0, J, K, pi, pjk, njk, mu, sigma2){
  # n0 - numeric value number of resistent observations
  # p0 - numeric value probability of being in interval B0
  # J -  numeric value Number of Bins
  # K -  numeric value Number of components
  # pi - vector of mixing proportions
  # pjk - matrix (jxk) occurence probabilty of a jth bin in the kth density
  # njk - matrix (jxk) - number of expected observation in jth bin under kth density
  # mu - vector of mean of normals
  # sigma2 - vector of variance of normals
  #OUTPUT
  #likelihood value (loglik)

  #
  term1 <- n0 * log(p0)

  if(K == 1){
    term2 <- sum(njk * log(pi)) # should be 0
  }else{
    logpi <- log(pi)
    logpi[logpi == -Inf] <- log(.Machine$double.eps)
    term2 <- sum(colSums(njk) * logpi)
  }

  log_pjk <- log(pjk)

  log_pjk <- apply(log_pjk, 2, function(x){
    ifelse(x == -Inf, log(.Machine$double.eps), x)
  })



  term3 <- sum(njk * log_pjk) #instead of sum sum




  #Term 2
  loglik_value <- term1 + term2 + term3

  return(loglik_value)
}

loglik2 <- function(y, pi, mu, sigma2, ab_bin){
  ll <- numeric(length(y))
  for(i in 2: length(y)){
    ll[i] <- y[i] * log(sum(
      sapply(1:length(mu), function(x){
      pi[x]*(pjk(a = ab_bin$a[i],
                 b = ab_bin$b[i],
                 mu = mu[x],
                 sigma2 = sigma2[x]))
    })))
  }

  return(sum(ll))
}

# Vector of mixing proportions
pi <- function(N, n0, njk){
  # N - numeric value number of all observations
  # n0 - numeric value number of resistance observtions
  # njk - dataframe ((j-1)xk) with expected values per bin and component
  #OUTPUT
  # numeric vector with pi for each component k

  if(is.null(dim(njk))){ # K = 1
    return(sum(njk)/(N-n0))
  }else{
    return(colSums(njk)/(N-n0))
  }
}

#Function for Optimize penalized Logliklihood with optim
optim.loglik.pen <- function(musigma, J, njk, ab_bin, alpha, beta){
  # musigma, vector with length 2*k, first k = mu, last k= sigma2
  # njk vector length J
  mu <- musigma[1]
  sigma2 <- musigma[2]

  pjk_exp <- numeric(J)
  for(j in 1:J){#do it for each bin (44)
    pjk_exp[j] <-
      pjk(a = ab_bin$a[j],
          b = ab_bin$b[j],
          mu = mu,
          sigma2 = sigma2)
  }

  log_pjk <- log(pjk_exp)
  log_pjk[log_pjk == -Inf]<- log(.Machine$double.xmin)

 # if(sigma2 < 0.01){
#    sigma2 <- 0.01
 # }

 # loglik.pen <-  sum(log_pjk* njk) - log(dinvgamma(sigma2, alpha, beta))
  loglik.pen <-  sum(log_pjk[-1]* njk[-1]) + invgamma::dinvgamma(sigma2, alpha, beta, log = TRUE)

  return((-1)*loglik.pen)
}

# Calculate expected number of observations in bin j nd k
njk <- function(nj, pi, mu, sigma2, a, b,  k){
  # nj - Observed number of observations in bin j
  # mu - vector of mean of normals
  # sigma2- vector of variance of normals
  # pi - vector of mixing proportions of normals
  # k - indizes of kth component
  #OUTPUT
  # numeric number - expected number of observations
  denum <- NULL
  for(i in 1:length(mu)){
   denum <- c(denum,
     pi[i] * pjk(a, b, mu[i], sigma2[i])
   )

    num <- nj* pi[k]* pjk(a, b, mu[k], sigma2[k])
  }

  return(num / sum(denum))

}


########################
## EM - Algorithm ######
########################
#' @title Estimating Ecoff by mixing Gaussians
#' @description
#' This function tries to find the Ecoff value of the input data, by aproximating the distribution
#' by an mixing distribution of Gaussian distributions
#'
#' @param y A data vector with observed values per bin.
#' @param mu A vector with startvalues for mu
#' @param sigma2 A vector with startvalues for the variance
#' @param pi A vector with startvalues for the mixing proportions
#' @param alpha Shape parameter of inverse gamma distribution
#' @param beta Scale parameter of inverse gamma distribution
#' @param epsilon Convergence criterium
#' @param ecoff.quantile Quantil which defines the ECOFF
#' @param max.iter Maximum number of iterations
#'
#' @return A list with components mu, var, pi, loglik and ecoff
#' @details
#' The data vector y is defined as one row of the EUCAST Data of the Zone Diameter Data.
#' The first value of y has to be the number of resistant observations.
#' Then the following elements are the number of observations in bin 6 to 6+length(y)-1.
#' A bin x includes all observations which had values form x-0.5 to x+0.5.
#'
#' The function uses the EM Algorithm to fit a mixing distribution of Normals on the data.
#' Based on the result of the converged mixing distribution the algorithm evaluates the ECOFF value.
#' The ECOFF value is defined by default as the 0.01 quantile of the rightest distribution.
#' The rightest density is defined as the distribution, where the sum of the pis firstly is greater than 0.3.
#' Therefore, the distribution get ordered by their mean in ascending order and then the pis are commulated of the right.
#'
#' Furthermore,  to avoid, that one of the sigma2 converges to 0, the likelihood get penalized in dependence of sigma2.
#' In detail, the penalization term follows a inverse gamma distribution with parameter alpha and beta.
#' @examples
#' y <- c(2, 4, 5,6,5,2,2, 1, 1, 2,  2, 1,6,7,8,7,6, 5, 2,1)
#' em.gauss(y = y,
#'          mu = c(8.36, 17.28),
#'          sigma2 = c(1.67,  8.4),
#'          pi = c(1/2, 1/2),
#'          alpha = 1,
#'          beta = 3,
#'          epsilon = 0.0001)

#' @export
em.gauss <- function(y, mu, sigma2, pi, alpha, beta,
                     epsilon=0.000001,ecoff.quantile=0.01, max.iter = 1000){

  # y - data numeric vector with observation per bin,
  # n0  - first value of y must be n0
  # mu - vector of mean of normals
  # sigma2- vector of variance of normals
  # pi - vector of mixing proportions of normals
  # alpha - numeric number alpha of inverse gauss
  # beta - numeric number beta of inverse gauss
  # epsilon - stopping criteria (change of likelihood)
  #OUTPUT
  # list(estimated mu, estimated sigma2, liklihood value)


#  if(class(y) != "numeric") stop("y is not a numeric vector")
  if(any(y<0)) stop("only positive y values are allowed")

  if(sum(y[-1])<10) stop("EM-Algorithm needs at least 10 observations")

  if(!is.numeric(mu)) stop("mu is not a numeric vector")
  if(any(mu<0)) stop("only positive mu values are allowed")

  if(!is.numeric(sigma2)) stop("sigma2 is not a numeric vector")
  if(any(sigma2<0)) stop("only positive sigma2 values are allowed")

  if(!is.numeric(pi)) stop("pi is not a numeric vector")
  if(any(pi<0)) stop("only positive pi values are allowed")

  if(length(mu) != length(sigma2) || length(sigma2) != length(pi)){
    stop("mu and sigma2 or sigma2 and pi have not the same length")
  }

  if(length(y) <= length(mu)){
    stop("y must be at least the same length as mu")
  }

  if(!is.numeric(alpha) || length(alpha) != 1 || any(alpha<0)){
    stop("alpha is not a positive numeric value")
  }

  if(!is.numeric(beta)|| length(beta) != 1 || any(beta<0)){
    stop("beta is not a positive numeric value")
  }

  if(!is.numeric(epsilon) || length(epsilon) != 1 || any(epsilon<0)){
    stop("epsilon is not a positive numeric value")
  }

  #Initialize
  K <- length(mu) # K number of components
  J <- length(y) # J number of Bins
  n0 <- y[1] # number of restistant observations
  N <- sum(y) # Total number of observations
  loglik_prev <- 0 #Initialize previous loglik for first iteration
  delta <- epsilon + 1 #Initialize to go at least one iteration in em_algorithm


  mu_est <- mu
  sigma2_est <- sigma2
  pi_est <- pi

  n0 <- y[1]

  is0 <- (y == 0)
  y_org <- y
  y <- y[!is0]

  # a = lower bound of bin , b= upper bound
  ab_bin<- data.frame(y=y_org,
                      a = (1:length(y_org))+ 4.5,
                      b= (1:length(y_org))+ 5.5)

  ab_bin$a[1] <- 0 #Set the first interval from 0 to 6.5
  ab_bin <- ab_bin[!is0, ]

  J <- length(y)
  iter <- 0
  while(delta > epsilon & iter < max.iter){
 # for(i in 1:10){
    iter <- iter +1
    #E- Step
    #Calculate Expected values in bin j under distribution k
    njk_exp <- matrix(0, nrow = J, ncol= K)
    for(j in 1:J){#do it for each bin (44)
      njk_exp[j,] <- sapply(1:K, FUN = function(k){
        njk (nj = y[j],
             pi = pi_est,
             mu = mu_est,
             sigma2 = sigma2_est,
             a = ab_bin$a[j],
             b = ab_bin$b[j],
             k = k)
      })
    }





    #M - Step

    # Estimate pi
    pi_est <- pi(N = N,
                 n0 = n0,
                 njk = njk_exp[-1, ]) #njk has to be without n0

 #should be 1
    # make a pjk matrix for likelihood

    #GRUEN: Calculate denom of pjk firstly
    # Do not recalculate each time
    pjk_exp <- matrix(0, nrow = J, ncol= K)
    for(j in 1:J){#do it for each bin (44)
      pjk_exp[j,] <- sapply(1:K, FUN = function(k){
        pjk(a = ab_bin$a[j],
            b = ab_bin$b[j],
            mu = mu_est[k],
            sigma2 = sigma2_est[k])
      })
    }


    # Optimize Likelihood (Estimate sigma, mu with optim)
    # Use as start values estimates from before

    p0 <- pj(p0 = -1,
             pi = pi_est,
             a = 0,
             b = 6.5,
             mu = mu_est ,
             sigma2 = sigma2_est,
             get.p0= TRUE)
    if(p0 == 0){
      p0 <- .Machine$double.eps
    }

    for(k in 1:K){
      musigma <- c(mu_est[k], sigma2_est[k])

      est1 <- stats::optim(par = musigma,
                    fn = optim.loglik.pen,
                    method="L-BFGS-B",
                    lower= c(-Inf, 0.1),
                    J = J,
                    njk = njk_exp[, k],
                    ab_bin = ab_bin,
                    alpha= alpha,
                    beta = beta
                    )

      mu_est[k] <- est1$par[1]
      sigma2_est[k] <- est1$par[2]


    }



    # Check loglikelihood

    loglik_curr <- loglik(n0 = n0,
                          p0 = p0,
                          J = J,
                          K = K,
                          pi =pi_est,
                          pjk = pjk_exp,
                          njk = njk_exp,
                          mu = mu_est,
                          sigma2 = sigma2_est)


    if(loglik_prev == -Inf & loglik_curr == -Inf){
      delta <- 0
    }else{
      delta <- abs(loglik_curr - loglik_prev)
    }
    loglik_prev <- loglik_curr


  }
  ecoff <- ecoff(mu_est=mu_est, pi_est = pi_est,
                 sigma2_est = sigma2_est, quantile = ecoff.quantile)

loglik.test <- loglik2(y = y,
                         mu = mu_est,
                         sigma2 = sigma2_est,
                         pi = pi_est,
                         ab_bin = ab_bin)


  return(list(mu= mu_est, sigma2 =sigma2_est, pi = pi_est, loglik = loglik.test, ecoff=ecoff))

}



#' @title Finding the optimal number of components
#' @description
#' This function provides the user information of the fitted distributions in order to choose the optimal number of components.
#'
#' @param y A data vector with observed values per bin.
#' @param k maximum numbers of components
#' @param alpha inverse gamma shape parameter
#' @param beta inverse gamma rate parameter
#' @param method method how startvalues should be evaluated. For more details see function CreateCluster
#' @param epsilon stopping criterion
#'
#' @return A list with  mu, var, pi, loglik, ecoff, AIC, BIC for each fitted distribution
#' @details
#' This function fits for each number of components (1:k) a mixing distribution of Gaussians by using
#' the function em.gauss. Furthermore, the fit of the distributions is measured by the two information criteria AIC and BIC.
#'
#' @examples
#' y <- c(2, 4, 5,6,5,2,2, 1, 1, 2,  2, 1,6,7,8,7,6, 5, 2,1)
#' em.gauss.opti.groups(y,
#'                      k= 2,
#'                      alpha = 1,
#'                      beta = 2)
#' @export
em.gauss.opti.groups <- function(y, k, alpha, beta, method = "quantile", epsilon=0.000001){

  # y - data numeric vector with observation per bin,
  # alpha - numeric number alpha of inverse gauss
  # beta - numeric number beta of inverse gauss
  # epsilon - stopping criteria (change of likelihood)
  # OUTPUT
  # list(estimated mu, estimated sigma2, likelihood value)
  # method quantiles and binbased

  if(any(y<0)) stop("only positive y values are allowed")

  if(!is.numeric(k)) stop("k is not a numeric vector")
  if(any(k<0)) stop("only positive k values are allowed")

  if(!is.numeric(alpha) || length(alpha) != 1 || any(alpha<0)){
    stop("alpha is not a positive numeric value")
  }

  if(!is.numeric(beta) || length(beta) != 1 || any(beta<0)){
    stop("beta is not a positive numeric value")
  }

  m <- matrix()
  goodness <- list()

  aic.vec <- vector("numeric",length = k)
  bic.vec <- vector("numeric",length = k)

  for(i in 1:k){


    y.df <- data.frame(bin = (1:length(y)+6), nrObs = y)

    m <- createCluster(y = y.df , k = i , method = method)

     goodness[[i]]<- em.gauss(y=y,
                  mu=m[,1],
                  sigma2=m[,2],
                  pi= rep(1/i, i),
                  alpha=alpha,
                  beta=beta,
                  epsilon=epsilon)

     aic <- AIC.gauss(goodness[[i]]$loglik, i)
     bic <- BIC.gauss(goodness[[i]]$loglik, i, sum(y))

     aic.vec[i] <- aic
     bic.vec[i] <- bic

  }
  goodness$AIC <- aic.vec
  goodness$BIC <- bic.vec

 # goodness[[k+1]] <- aic.vec
 # goodness[[k+2]] <- bic.vec

  return(goodness)

}

AIC.gauss <- function(lik, par){

  #lik - numeric value describing the likelihood
  #par - number of groups on the basis of normal distribution

  aic <- -2*lik+2*(3* par -1)

  return(aic)
}

BIC.gauss <- function(lik, par, n){

  #lik - numeric value describing the likelihood
  #par - number of groups on the basis of normal distribution
  #n - number of observations

  bic <- -2*lik + 2* (3* par -1) * log(n)

  return(bic)
}


ecoff <- function(mu_est, pi_est, sigma2_est,quantile=0.01) {
  val <- 0
  for(i in length(mu_est):1) {
    val <- val + pi_est[i]
    if(val > 0.3) break
  }
  return(stats::qnorm(quantile,mean=mu_est[i], sd = sqrt(sigma2_est[i])))
}

plot.dens <- function(x, mu, sigma2, pi){
  dens <- 0

  for(i in 1:length(mu)){
    dens <- dens + pi[i]* stats::dnorm(x, mean = mu[i], sd = sqrt(sigma2[i]))
  }
  return(dens)
}


#' @title Plotting the fitted distribution with ECOFF
#' @description This function plotts a histogram of the data and the fitted mixing distribubion of normals.
#' Furthermore, the Ecoff value is plotted by a line.
#'
#' @param y A data vector with observed values per bin
#' @param mu A vector with values for mu
#' @param sigma2 A vector with values for the variance
#' @param pi A vector with values for the mixing proportions
#' @param ecoff a numeric value of the ECOFF value
#'
#' @return a plot of the fitted distribution
#'
#' @export
plot_fct <- function(y, mu, sigma2, pi, ecoff) {
  if(any(y<0)) stop("only positive y values are allowed")

  if(!is.numeric(mu)) stop("mu is not a numeric vector")
  if(any(mu<0)) stop("only positive mu values are allowed")

  if(!is.numeric(sigma2)) stop("sigma2 is not a numeric vector")
  if(any(sigma2<0)) stop("only positive sigma2 values are allowed")

  if(!is.numeric(pi)) stop("pi is not a numeric vector")
  if(any(pi<0)) stop("only positive pi values are allowed")

  if(length(mu) != length(sigma2) || length(sigma2) != length(pi)){
    stop("mu and sigma2 or sigma2 and pi have not the same length")
  }

  if(!is.numeric(ecoff)) stop("ecoff is not a numeric value")


  y.data <- data.frame(name = 1:length(y)+5, y)
  lim=max(graphics::hist(y,breaks = 30, freq = F, plot = F)$density)
  graphics::hist(rep(y.data[,1],y.data[,2]),freq = F , col = "deepskyblue", xlab="mm",
       main = paste("Gaussian Mixtures with ", length(mu), " components"),breaks=30,
       xlim=c(5,max(y.data[1])))

  graphics::curve(plot.dens(x,
                  mu,
                  sigma2,
                  pi), from= 6, to = 50, add = T, ylab = 'density')
  graphics::abline(v=ecoff, col="red", lwd=3,lty=2)
#  legend("topleft",paste("CUTOFF = ", round(ecoff,2), " mm"), col="red", cex=1, lwd=2, lty=2)
}
