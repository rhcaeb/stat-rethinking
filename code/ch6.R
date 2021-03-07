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
