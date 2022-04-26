# Calibrated inference: inference that accounts for both sampling and distributional uncertainty


This package provides functions for calibrating linear models to account for both sampling and distributional uncertainty. In addition, helper functions are provided to simplify simulating from the distributional perturbation model. Details on the theoretical background can be found in [Jeong and Rothenhaeusler (2022)](https://arxiv.org/abs/2202.11886).

## How to install

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using  `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("rothenhaeusler/calinf")`.

## Usage

First, one has to set a distributional seed. Then, one can draw from the perturbed distribution.

```R
d_seed <- distributional_seed(n=1000,delta=10)
# Sample from perturbed distribution 
# drnorm is the equivalent of rnorm for distributional perturbations
x <- drnorm(d_seed)

# In the qqplot there are clear deviations from the standard Gaussian distribution
qqnorm(x)
# Compare to standard normal
qqnorm(rnorm(1000))
```

Drawing from a multivariate perturbed distribution.

```R
#Set distributional seed
d_seed <- distributional_seed(n=1000,delta=2)
#Draw from perturbed model
x <- drnorm(d_seed)
y <- drnorm(d_seed)
```

Evaluation whether it works for multidimensional data generation.

```R
n <- 1000
delta <- 1.3
draw_perturbed_data <- function(n,delta){
  #Set distributional seed
  d_seed <- distributional_seed(n,delta)
  #Draw from perturbed model
  x <- drnorm(d_seed)
  y <- drnorm(d_seed)
  c(mean(x),mean(x^2),mean(x^3),mean(y),mean(y^2),mean(x*y))
}
mat <- sapply(rep(n,10000),draw_perturbed_data,delta=delta)
vec <- apply(mat,1,var)
# target delta^2
delta^2
# achieved delta^2 on simulated data
n*vec[1]/var(rnorm(10000))
n*vec[2]/var(rnorm(10000)^2)
n*vec[3]/var(rnorm(10000)^3)
n*vec[4]/var(rnorm(10000))
n*vec[5]/var(rnorm(10000)^2)
n*vec[6]/var(rnorm(10000)*rnorm(10000))
```

Calibrated inference for linear models.

```R
n <- 1000
X <- rnorm(n)
Z1 <- rnorm(n)
Z2 <- rnorm(n)
Y <- 1*X + X^2 + Z1 + Z2 + rnorm(n)
df <- data.frame(cbind(X,Y,Z1,Z2))
data <- as.data.frame(cbind(Y,X))
formulas <- list(Y~X, Y ~ X + I(X^2), Y ~ X + Z1, Y ~ X + Z2 + I(X^2), Y ~ X + Z1 + Z2 + I(X^2))
calm(formulas,data=data, target="X")
#summary(lm(Y~X + Z1 + Z2 + X^2,data=df))
```

## Looking behind the curtain: a worked out example

In this example, we show in detail the algebra behind calibrated inference. We assume that the population mean of Z1,...,Z4 is known to the data scientist through background knowledge. This can be used to calibrate the uncertainty of estimating the mean of X.


```R
delta <- 10

# Sample from a multivariate perturbation model
sample_data <-  function(){
d_seed <- distributional_seed(n=100,delta=delta)
X <- drnorm(d_seed,mean=2,sd=2)
Z1 <- drnorm(d_seed,mean=2.2,sd=2)
Z2 <- drnorm(d_seed,mean=3.2,sd=3)
Z3 <- drnorm(d_seed,mean=5.3,sd=4)
Z4 <- drnorm(d_seed,mean=1,sd=2)
df <- as.data.frame(cbind(X,Z1,Z2,Z3,Z4))
return(df)
}

evaluate_distributional_ci <- function(i){
df <- sample_data() 
# constructing distributional confidence intervals with nominal coverage .95
# Here, we use background knowledge about the population mean of Z1...Z4.
sigmasq = var(df$X)*(mean(df$Z1-2.2)^2/var(df$Z1) + mean(df$Z2-3.2)^2/var(df$Z2) + mean(df$Z3-5.3)^2/var(df$Z3) + mean(df$Z4-1)^2/var(df$Z4))/4
tquantile <- qt(p=.975,df=4)
ci <- c(mean(df$X)-tquantile*sqrt(sigmasq),mean(df$X)+tquantile*sqrt(sigmasq))
}

confidence_intervals <- sapply(1:1000,evaluate_distributional_ci)
mean( confidence_intervals[1,] <= 2 & 2 <= confidence_intervals[2,])


evaluate_sampling_ci <- function(df){
  df <- sample_data() 
  # constructing sampling confidence intervals with nominal coverage .95
  sigmasq = var(df$X)/sqrt(n)
  tquantile <- qt(p=.975,df=999)
  ci <- c(mean(df$X)-tquantile*sqrt(sigmasq),mean(df$X)+tquantile*sqrt(sigmasq))
}
confidence_intervals <- sapply(1:1000,evaluate_sampling_ci)
mean( confidence_intervals[1,] <= 2 & 2 <= confidence_intervals[2,])
```
