# Calibrated inference: inference that accounts for both sampling and distributional uncertainty


This package provides functions for calibrating linear models to account for both sampling and distributional uncertainty. In addition, helper functions are provided to simplify simulating from the distributional perturbation model. Details on the theoretical background can be found in [Yujin and Rothenhaeusler (2022)](https://arxiv.org/abs/2202.11886).

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
formulas <- list(Y~X, Y ~ X + I(X^2), Y ~ X + Z1, Y ~ X + Z2 + X^2, Y ~ X + Z1 + Z2 + X^2)
calm(formulas,data=data, target="X")
#summary(lm(Y~X + Z1 + Z2 + X^2,data=df))
```
