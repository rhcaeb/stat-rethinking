# Chapter 2 ----

## 2.1 - The Garden of Forking Data ----

ways <- c(0, 3, 8, 9, 0) # ways to produce data
ways/sum(ways) # plausibility

## 2.3 - Components Of the Model

dbinom(x = 6, size = 9, prob = 0.5) # binomial distribution;
# likelihood of data: '6' water in '9' tosses (N), under any value of p.

## 2.4 - Making the Model 'Go' ----

### 2.4.3 - Grid Approximation

# define grid
p_grid <- seq(from = 0, to = 1, length.out = 20) 

# define prior
prior <- rep(1, 20) 

# computer likelihood at each value in grid
likelihood <- dbinom(x = 6, size = 9, prob = p_grid) 

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior 

# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# display the posterior distribution
plot(p_grid, posterior, type = 'b',
     xlab = 'probability of water',
     ylab = 'posterior probability')
mtext('20 points')

# replicate different priors in fig 2.5
prior <- ifelse(p_grid < 0.5, 0, 1)
prior <- exp(-5*abs(p_grid - 0.5))

### 2.4.4 - Quadratic Approximation

library(rethinking)

globe.qa <- quap(
  alist(
    W ~ dbinom(W+L, p), # binomial likelihood
    p ~ dunif(0, 1) # uniform prior
  ),
  data = list(W = 6, L = 3)
)

# display the summary of quadratic approximation
precis(globe.qa)

# analystical calculation
W <- 6
L <- 3
curve(dbeta(x, W+1, L+1), from = 0, to = 1)

# quadratic approximation
curve(dnorm(x, 0.67, 0.16), lty = 2, add = TRUE)

### 2.4.5 - Markov chain Monte Carlo

n_samples <- 1000
p <- rep(NA, n_samples)
p[1] <- 0.5
W <- 6
L <- 3

for (i in 2:n_samples){
  p_new <- rnorm(1, p[i-1], 0.1)
  if(p_new < 0) p_new <- abs(p_new)
  if(p_new > 1) p_new <- 2-p_new
  q0 <- dbinom(W, W+L, p[i-1])
  q1 <- dbinom(W, W+L, p_new)
  p[i] <- ifelse(runif(1) < q1/q0, p_new, p[i-1])
}

dens(p, xlim = c(0,1))
curve(dbeta(x, W+1, L+1), lty=2, add = TRUE)
