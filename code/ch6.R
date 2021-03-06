# Chapter 6 - The Haunted DAG and Causal Error

## 6.2  ---- 

# Multicollinear legs

N <- 100 # n indviduals
set.seed(909)
height <- rnorm(N, 10, 2) # sim total height of each
leg_prop <- runif(N, 0.4, 0.5) # leg as proportion of height
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # sim left leg as prop + error
leg_right <- leg_prop*height + rnorm(N, 0, 0.02) # sim right reg as prop + error
d <- data.frame(height, leg_left, leg_right)

# 6.3 ----

# predict outcome height with both predictors
m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2, 10),
    sigma ~ dexp(1)
  ),
  data = d)
precis(m6.1)

# 6.4 ----

# graphical view of precis output
plot(precis(m6.1))

## 6.5 ----

#bivariate posterior distribution
post <- extract.samples(m6.1)
plot(bl ~ br, post, col = col.alpha(rangi2, 0.1), pch = 16)

## 6.6 - compute posterior distribution of their sum (beta l + beta r)
sum_blbr <- post$bl + post$br
dens(sum_blbr, col = rangi2, lwd = 2, xlab = 'sum of bl and br')

## 6.7 - fit regression w/ one leg length variable
m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.2)

## 6.8 - multicollinear (primate milk e.g.)
data(milk)
d <- milk
d$K <- scale(d$kcal.per.g)
d$F <- scale(d$perc.fat)
d$L <- scale(d$perc.lactose)

## 6.9 - model energy content as a function of percent fat + percent lactose
# kcal.per.g regressed on percent fat
m6.3 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bF*F,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)

# kcal.per.g regressed on percent lactose
m6.4 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bL*L,
    a ~ dnorm(0, 0.2),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)

precis(m6.3)
precis(m6.4)

## 6.10 - place both predictor variables in the same regression model
m6.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bF*F + bL*L,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)

precis(m6.5) # sd are twice as large as in the bivariate models

## 6.11 - use pairs plot; variables form essentially a single axis of variation
pairs(~ kcal.per.g + perc.fat + perc.lactose, data = d, col = rangi2)

## 6.13 - post treatment e.g.

set.seed(71)
# number of plants
N <- 100

# simulate initial heights 
h0 <- rnorm(N, 10, 2)

# assign treatments and simulate fungus and growth
treatment <- rep(0:1, each = N/2)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment * 0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose a clean data frame
d <- data.frame(h0 = h0, h1 = h1, treatment = treatment, fungus = fungus)
precis(d)

## 6.14 - prior dist
sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))

## 6.15 - fit model
m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.6)

## 6.16 - approx. posterior
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.7)

## 6.17 - omit post-treatment variable 'fungus'
## 'blocked by consequence' 
m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.8)

## 6.2.3 - d-separation
## 6.18 - look at the above problem in terms of a DAG
library(dagitty)
plant_dag <- dagitty("dag{
                     H_0 -> H_1
                     F -> H_1
                     T -> F}")
coordinates( plant_dag ) <- list( x=c(H_0=0,T=2,F=1.5,H_1=1) ,
                                  y=c(H_0=0,T=0,F=0,H_1=0) )
drawdag( plant_dag )

## 6.19 - query the implied conditional independecies for this DAG
impliedConditionalIndependencies(plant_dag)

## 6.20 - include F in the model & modify sim -> no influence on growth
set.seed(71)
N <- 1000
h0 <- rnorm(N, mean = 10, sd = 2)
treatment <- rep(0:1, each = N/2)
M <- rbern(N)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment*0.4 + 0.4*M)
h1 <- h0 + rnorm(N, 5 + 3*M)
d2 <- data.frame(h0 = h0, h1 = h1, treatment = treatment, fungus = fungus)

##### Collidar bias ----

library(rethinking)
d <- sim_happiness(seed = 1977, N_years = 1000)
precis(d)

## 6.22

d2 <- d[d$age > 17, ] # only adults
d2$A <- (d2$age - 18) / (65-18)

## 6.23 

d2$mid <- d2$married + 1
m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
       
    a[mid] ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2)
precis(m6.9, depth = 2)

# compare the above model with one that omits marriage status

m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2)
precis(m6.10)

## 6.25 

## simple simulation -> invent some strength of association
N <- 200 # number of grandparent-parent-child triads
b_GP <- 1 # direct effect of G on P
b_GC <- 0 # direct effect of G on C
b_PC <- 1 # direct effect of P on C
b_U <- 2 # direct effect of U on P and C

## use the above 'slopes' to draw random observations
set.seed(1)
U <- 2*rbern(N, 0.5) - 1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U)
C <- rnorm(N, b_PC*P + b_GC*G + b_U*U)
d <- data.frame(C = C, P = P, G = G, U = U)

## simple regression model of C on P and G
m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.11)

## need to measure 'U' for collider bias; regress on U
m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0,1),
    c(b_PC, b_GC, b_U) ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.12)

