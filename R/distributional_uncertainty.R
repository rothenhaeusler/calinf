#' Function for creating distributional seeds
#' @param n The sample size
#' @param delta The pertubation parameter
#' @examples
#' d_seed <- distributional_seed(1000,10)
#' qqnorm(drnorm(d_seed))
#' # Compare to standard normal
#' qqnorm(rnorm(1000))
#' @export
distributional_seed <- function(n,delta){
  ## Set distributional seed
  m <- round(n/(delta^2-1),0)
  #Exact choice of delta is often not possible.
  #print(paste("Rounding to delta = ",sqrt(1+n/m),"."),sep="")
  d_seed <- sample.int(m,size=n,replace=TRUE)
  return(d_seed)
}

#' Distributional random number generation (equivalent to rnorm)
#' @param d_seed A distributional seed as returned by distributional_seed()
#' @param mean vector of means of the unperturbed Gaussian distribution
#' @param sd vector of standard deviations of the unperturbed Gaussian distribution
#' @examples
#' d_seed <- distributional_seed(1000,10)
#' qqnorm(drnorm(d_seed))
#' # Compare to standard normal
#' qqnorm(rnorm(1000))
#' n <- 1000
#' delta <- 1.3
#' # Test whether it works for multi-dimensional data generation
#' draw_perturbed_data <- function(n,delta){
#'   #Set distributional seed
#'   d_seed <- distributional_seed(n,delta)
#'   #Draw from perturbed model
#'   x <- drnorm(d_seed)
#'   y <- drnorm(d_seed)
#'   c(mean(x),mean(x^2),mean(x^3),mean(y),mean(y^2),mean(x*y))
#' }
#' mat <- sapply(rep(n,10000),draw_perturbed_data,delta=delta)
#' vec <- apply(mat,1,var)
#' # target delta^2
#' delta^2
#' # achieved delta^2 on simulated data
#' n*vec[1]/var(rnorm(10000))
#' n*vec[2]/var(rnorm(10000)^2)
#' n*vec[3]/var(rnorm(10000)^3)
#' n*vec[4]/var(rnorm(10000))
#' n*vec[5]/var(rnorm(10000)^2)
#' n*vec[6]/var(rnorm(10000)*rnorm(10000))
#' @export
drnorm <- function(d_seed,mean=0,sd=1){
  unique <- rnorm(max(d_seed))
  rv <- unique[d_seed] + 1/sqrt(length(unique(d_seed)))*rnorm(length(d_seed))
  return(mean + rv*sd)
}


#' Calibrated inference for linear models
#' @param formulas A list of formulas
#' @param data A dataframe containing the variables in the model
#' @examples
#' n <- 1000
#' X <- rnorm(n)
#' Z1 <- rnorm(n)
#' Z2 <- rnorm(n)
#' Y <- 1*X + X^2 + Z1 + Z2 + rnorm(n)
#' df <- data.frame(cbind(X,Y,Z1,Z2))
#' data <- as.data.frame(cbind(Y,X))
#' formulas <- list(Y~X, Y ~ X + I(X^2), Y ~ X + Z1, Y ~ X + Z2 + X^2, Y ~ X + Z1 + Z2 + X^2)
#' calm(formulas,data=data, target="X")
#' summary(lm(Y~X + Z1 + Z2 + X^2,data=df))
#' @export
calm <- function( formulas, data, target, ...){

  n <- nrow(data)
  infls <- vector(mode="list",length = length(formulas))
  coefs <-  vector(mode="list",length = length(formulas))

  for (j  in 1:length(formulas)){
    model <- lm(formulas[[j]],data=data)
    # Remark: the influence() function also exists for glms! We should be able to create a calibratedglm function analogously
    infls[[j]] <- n*influence(model)$coefficients[,target]
    coefs[[j]] <-  coef(model)[target]
  }

  infls <- simplify2array(infls)
  coefs <- simplify2array(coefs)

  # Careful! I haven't tested the coverage of this extensively yet! This is just a prototype.
  estimate <- mean(coefs)
  sd_estimate <- sd(rowMeans(infls))
  delta <-   sqrt(n)*sqrt(mean( (coefs - mean(coefs))^2/ apply( infls - rowMeans(infls),2,var)))

  p_values <- 2*pnorm(abs(estimate)/(sd_estimate*delta/sqrt(n)),lower.tail=FALSE)

  ret_table <- cbind(estimate,sd_estimate*delta/sqrt(n),p_values)
  ret_table <- round(ret_table,4)
  colnames(ret_table) <- c("Estimate", "Std. Error", "Pr(>|z|)")
  rownames(ret_table) <- target
  cat("\n")
  cat("Quantification of both distributional and sampling uncertainty")
  cat("\n\n")
  print(ret_table)
  cat("\n")

}


