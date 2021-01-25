# Chapter 4 - Geocentric Models ----

## 4.1 Normal by addition
pos <- replicate(1000, sum(runif(16, -1, 1)))
hist(pos) # plot

## 4.2 - Normal by multiplication
prod(1+runif(12, 0, 0.1)) # sample a random growth rate

## 4.3
growth <- replicate(10000, prod(1 + runif(12, 0, 0.1)))
dens(growth, norm.comp = TRUE)

## 4.4 - small effects that multiply togther are approx. additive (stabilize Gaussian)
big <- replicate(10000, prod(1 + runif(12, 0, 0.5)))
small <- replicate(10000, prod(1 + runif(12, 0 , 0.1)))

## 4.1.3 - Normal by log-multiplication
log.big <- replicate(10000, log(prod(1 + runif(12, 0, 0.5))))
dens(log.big)
