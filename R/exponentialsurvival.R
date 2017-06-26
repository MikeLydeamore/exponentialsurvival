#' Fit Parametric Exponential Mixture Model
#' 
#' Fits a parametric mixture model of the form \eqn{S(t) = \alpha + (1-\alpha)exp(-\lambdat)} to (right-censored) data.
#' 
#' This form of model allows for some proportion of individuals to experience zero risk (\eqn{\alpha}), and the remainder to experience constant hazard (\eqn{\lambda}).
#' @param surv A matrix with the first column being event times (or censor time), and the second being whether the event was censored. 
#' 
#' Note: For censoring, 0=alive (censored), 1=dead
#' @return Data frame with two elements: Zero risk proportion, \eqn{\alpha}, and estimated hazard, \eqn{\lambda}.
#' @export
#' @examples 
#' surv <- expsurv.simulate(alpha = 0.3, lambda = 1, N=1000, duration = 15)
#' fit <- expsurv.fit(surv)

expsurv.fit <- function(surv) {
  #Times are in column 1, censoring in column 2.
  #Censoring: 0=alive, 1=dead.
  
  #Need the negative likelihood as nlm will only minimize.
  negativelikelihood <- function(parameters, data) {
    return (-1 * expsurv.likelihood(parameters[1], parameters[2], data))
  }
  
  MLE <- stats::nlm(f = negativelikelihood, p = c(0.5, 1), data=surv, gradtol = 1e-10)
  
  return (data.frame(alpha = MLE$estimate[1], lambda = MLE$estimate[2]))
}

#' Likelihood for a Parametric Exponential Mixture Model
#' 
#' Calculates the likelihood function for a Parametric Exponential Mixture Model
#' 
#' @param alpha Proportion of individuals whom experience 0 risk
#' @param lambda Constant hazard
#' @param data A matrix with the first column being event times (or censor time), and the second being whether the event was censored. 
#' 
#' Note: For censoring, 0=alive (censored), 1=dead
#' @export
#' @examples 
#' surv <- expsurv.simulate(alpha = 0.3, lambda = 1, N = 1000, duration = 15)
#' likelihood <- expsurv.likelihood(alpha = 0.3, lambda = 1, data = surv)

expsurv.likelihood <- function(alpha, lambda, data) {
  if (lambda < 0) {
    return (lambda * 500000)
  }
  if (alpha < 0 | alpha >= 1)
    return (-1*abs((alpha) * 500000))
  ret <- 0
  for (i in 1:nrow(data)) {
    if (data[i, 2] == 1) {
      ret <- ret + log( (1-alpha) * lambda * exp(-lambda*data[i, 1]))
    }
    else {
      ret <- ret + log((alpha + ( (1-alpha) * exp(-lambda*data[i, 1]) )))
    }
  }
  return (ret)
}

#' Simulate realisations from a Parametric Exponential Mixture Model
#' 
#' @param alpha Proportion of individuals who experience 0 risk
#' @param lambda Constant hazard
#' @param N Number of realisations (individuals in population)
#' @param duration Duration of the study
#' @return Data frame with 2 columns \code{time} and \code{status}, with one row for each realisation.
#' @export
#' @examples 
#' surv <- expsurv.simulate(alpha = 0.3, lambda = 1, N = 1000, duration = 15)

expsurv.simulate <- function(alpha, lambda, N, duration) {
  #Prob of infection = 1-a
  status <- as.numeric(stats::runif(n=N, min=0, max=1) < (1-alpha))
  durations <- stats::runif(n=N, min=0, max=duration)
  
  tab <- lapply(1:length(status), function(x) {
    if (status[x] == 0) {
      time <- durations[x]
    } else {
      time <- stats::rexp(n=1, rate=lambda)
    }
    ret <- data.frame(time = time, status = status[x])
  })
  
  return (do.call(rbind, tab))
}

#' Bootstrap a parametric exponential mixture model
#' 
#' Performs bootstrapping on the parametric exponential mixture model, to provide approximate measures on the uncertainty.
#' 
#' @param alpha Proportion of individuals who experience 0 risk
#' @param lambda Constant hazard
#' @param N Number of realisations (individuals in population)
#' @param duration Duration of the study
#' @param samples Number of samples to bootstrap.
#' 
#' Recommended number of samples is at least \code{N}.
#' 
#' @return Data frame with \code{samples} rows, with each row containing an estimate of \eqn{\alpha} and \eqn{\lambda} from a simulation.
#' @export
#' @examples 
#' bootstrapping <- expsurv.bootstrap(alpha = 0.3, lambda = 1, N = 1000, duration = 15, samples = 100)
expsurv.bootstrap <- function(alpha, lambda, N, duration, samples) {
  r <- lapply(1:samples, function(x) {
      sample <- expsurv.simulate(alpha, lambda, N, duration)
      return (expsurv.fit(sample))
    })
  
  
  return (do.call(rbind, r))
}