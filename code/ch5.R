# Chapter 5 - The Many Variables & Spurious Waffles ----

## 5.1 - Spurious association

# load data and copy 
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

# standardize variables
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)

## 5.2 - check how big a std. dev. of age at marriage is:
sd(d$MedianAgeMarriage)


## 5.3 - compute approximate posterior:
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)

## 5.4 - simulate from the priors:
set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post = prior, data = list(A = c(-2, 2))) # plot lines over the range of 2 SD for both outcome and predictor
plot(NULL, xlim = c(-2, 2), ylim = c(-2, 2))
for(i in 1:50) lines(c(-2, 2), mu[i, ], col = col.alpha("black", 0.4))

## 5.5 - posterior predictions

# compute percentile interval of mean
A_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.1, data = list(A = A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# plot
plot(D ~ A, data = d, col=rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)

## 5.6 - fit regression:
d$M <- scale(d$Marriage)
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)
