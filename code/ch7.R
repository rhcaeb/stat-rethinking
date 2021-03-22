# Chapter 7 - Ulysses' Compass

library(rethinking)

## E.g. overfitting

sppnames <- c("afarensis", "africanus", "habilis", "boisei",
              "rudolfensis", "ergaster", "sapiens")
brainvolcc <- c(438, 452, 612, 521, 752, 871, 1350)
masskg <- c(37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.0)
d <- data.frame(species = sppnames, brain = brainvolcc, mass = masskg)

## Fit series of complex models and see which functions fits the data 'best'

## scale params
d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain / max(d$brain)

## model fit
m7.1 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b*mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d) # R2 -> proportion of variance explained by the model


## E.g. ordinary least squares (OLS) and Bayesian anti-essentialism
## Use OLS to get posterior dist. from these brain size models:
## Note! You won't get a posterior for sigma....
m7.1_OLS <- lm(brain_std ~ mass_std, data = d)
post <- extract.samples(m7.1_OLS)


## Calculate R2
set.seed(12)
s <- sim(m7.1)
r <- apply(s, 2, mean) - d$brain_std
resid_var <- var2(r)
outcome_var <- var2(d$brain_std)
1 - resid_var/outcome_var

## Write a function to make this repeatable (for models to come...)
R2_is_bad <- function(quap_fit){
  s <- sim(quap_fit, refresh = 0)
  r <- apply(s, 2, mean) - d$brain_std
  1 - var2(r)/var2(d$brain_std)
}

## Build models (each will be a polynomial of higher degree) to compre to m7.1
m7.2 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 1),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b = rep(0,2)))

## The following (four) models are constructed
## Note, they are just third- fourth- fifth- and sixth degree polynomials

m7.3 <- quap(
  alist(
  brain_std ~ dnorm(mu, exp(log_sigma)),
  mu <- a + b[1]*mass_std + b[2]*mass_std^2 + b[3]*mass_std^3,
  a ~ dnorm(0.5, 1),
  b ~ dnorm(0, 10),
  log_sigma ~ dnorm(0, 1)
), data = d, start = list(b = rep(0,3)))

m7.4 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + b[3]*mass_std^3 +
    b[4]*mass_std^4,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b = rep(0,4)))

m7.5 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + b[3]*mass_std^3 +
      b[4]*mass_std^4 + b[5]*mass_std^5,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b = rep(0,5)))

m7.6 <- quap(
  alist(
    brain_std ~ dnorm(mu, 0.001),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + b[3]*mass_std^3 +
      b[4]*mass_std^4 + b[5]*mass_std^5 + b[6]*mass_std^6,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b = rep(0,6)))

## plot each model

## m7.1:
post <- extract.samples(m7.1)
mass_seq <- seq(from = min(d$mass_std), to = max(d$mass_std), length.out = 100)
l <- link(m7.1, data = list(mass_std = mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)
