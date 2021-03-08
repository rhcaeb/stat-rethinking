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
