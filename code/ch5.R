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

## 5.7 - Drawing a DAG
library(dagitty)
dag5.1 <- dagitty("dag{
                  A -> D
                  A -> M
                  M -> D}")
coordinates(dag5.1) <- list(x = c(A = 0, D = 1, M = 2), y = c(A = 0, D = 1, M = 0))
drawdag(dag5.1)

## 5.8 - Define the second DAG
DMA_dag2 <- dagitty('dag{D <- A -> M}')
impliedConditionalIndependencies(DMA_dag2) 

## 5.9 - First DAG has no conditional independencies (check):
DMA_dag1 <- dagitty('dag{D <- A -> M -> D}')
impliedConditionalIndependencies(DMA_dag1) # no output to display

## 5.10 - approximate the posterior (fit model to divorce rate):
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)
precis(m5.3)

## 5.11 - visualize posterior distributions for all three models focusing on slope params bA + bM:
plot(coeftab(m5.1, m5.2, m5.3), par = c("bA", "bM"))
