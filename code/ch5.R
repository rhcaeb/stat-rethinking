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

## 5.13 - predictor residual plots: approx. posterior
m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM * A,
    a ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)


## 5.14 - computer residuals by subtracting the obs M in each state from the pred. rate:
mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean

## 5.15 - posterior predictive check -> simulate predictions, averaging over the posterior
# call link without specifying new data
# so it uses original data
mu <- link(m5.3)

# summarize samples across cases
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

# simulate observations
# again, no new data, so uses original data
D_sim <- sim(m5.3, n = 1e4)
D_PI <- apply(D_sim, 2, PI)


## 5.16 - plot predictions against observed 
## and add a line to show prediction and line segments for the CI for each pred.
plot(mu_mean ~ d$D, col=rangi2, ylim = range(mu_PI),
     xlab = "Obs. divorce", ylab = "Pred. divorve")
abline(a = 0, b = 1, lty = 2)
for(i in 1:nrow(d)) lines(rep(d$D[i],2), mu_PI[, i], col = rangi2)

## 5.17 - label a few select points
identify(x = d$D, y = mu_mean, labels = d$Loc)

## 5.19 - add regression to the quap model running two regressions at the same time
data("WaffleDivorce")
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)

m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1),
    
    ## A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma_M ~ dexp(1)
  ), data = d)

precis(m5.3_A)

## 5.20 - define a range of values for A:
A_seq <- seq(from = -2, to = 2, length.out = 30)

## 5.21 - sim.observations
# prep data
sim_dat <- data.frame(A=A_seq)

# simulate M and then D, using A_seq
s <- sim(m5.3_A, data = sim_dat, vars = c("M", "D"))

## 5.22 - display counterfactual preds
plot( sim_dat$A , colMeans(s$D) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated A" , ylab="counterfactual D" )
shade( apply(s$D,2,PI) , sim_dat$A )
mtext( "Total counterfactual effect of A on D" )

## 5.23 - simulate a counterfactual for an avg. stage w/ A = 0..and see what changing M does
sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 )
s <- sim( m5.3_A , data=sim_dat , vars="D" )
plot( sim_dat$M , colMeans(s) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated M" , ylab="counterfactual D" )
shade( apply(s,2,PI) , sim_dat$M )
mtext( "Total counterfactual effect of M on D" )
