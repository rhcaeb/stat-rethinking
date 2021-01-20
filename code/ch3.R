# Chapter 3 ----

# 3.1 ----

Pr_positive_vampire <- 0.95
Pr_positive_mortal <- 0.01
Pr_vampire <- 0.001

Pr_positive <- Pr_positive_vampire * Pr_vampire +
               Pr_positive_mortal * (1 - Pr_vampire)

( Pr_vampire_positive <- Pr_positive_vampire * Pr_vampire / Pr_positive )


# 3.2 ----
# Sampling from a grid-approximate posterior

p_grid <- seq(from = 0, to = 1, length.out = 1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(x = 6, size = 9, prob = p_grid)
posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)

# 3.3 ----
# draw 10,000 samples from the posterior
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)

# plot samples
plot(samples)

# plot density estimate computed from samples
library(rethinking)
dens(samples)

# 3.6 ----
# Interval of defined boundaries
# posterior probability that the prop'n of water is < 0.5
# add up posterior probability where p < 0.5:
sum(posterior[p_grid < 0.5]) 

# perform same calculation using samples from the posterior:
sum(samples < 0.5) / 1e4 # add samples below 0.5, divide by total sample size

# how much posterior probability lies between 0.5 and 0.75:
sum(samples > 0.5 & samples < 0.75) / 1e4

# 3.9 ----
# boundaries of the lower 80% posterior probability
quantile(x = samples, 0.8)

# middle 80% interval (lies between 10th percentile and 90th percentile)
quantile(x = samples, c(0.1, 0.9))

# 3.1.1 ----
# E.g. observing three waters in three tosses and a flat (uniform) prior:
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(x = 3, size = 3, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)

# 50% percentile compatibility interval:
PI(samples, prob = 0.5)

# computer HPDI (50%)
HPDI(samples, prob = 0.5) # captures params with highest posterior probability

# 3.1.4 ----
# maximum a posteriori (MAP) estimate (parameter value w/ highest posterior probability)
p_grid[which.max(posterior)]

# 3.1.5 ----
# samples from the posterior
chainmode(samples, adj=0.01)

# 3.1.7 ----
# calculate expected loss
sum(posterior*abs(0.5 - p_grid))

# 3.1.8 ----
# trick to repeat this calculation for every possible decision:
loss <- sapply(p_grid, function(d) sum(posterior*abs(d - p_grid)))

# 3.1.9 ----
# Find the value (from 'loss') which minimizes the loss:
p_grid[which.min(loss)]

# 3.20 ----
# binomial likelihood (N = 2, tosses of the globe)
dbinom(0:2, size = 2, prob = 0.7)

# 3.21 ----
# simulate observations:
rbinom(1, size = 2, prob = 0.7)

# 3.22 ----
# set of 10 simulations:
rbinom(10, size = 2, prob = 0.7)

# 3.23 ----
# generate 100,000 dummy observations to clarify each value (0,1,2) appears
# in proportion to its likelihood:
dummy_w <- rbinom(1e5, size = 2, prob = 0.7)
table(dummy_w)/1e5 # calculated likelihoods

# 3.24 ----
# simulate the sample size from earlier (N = 9):
dummy_w <- rbinom(1e5, size = 5, prob = 0.7)
simplehist(dummy_w, xlab = 'dummy water count')

# 3.25 ----
# simulate predicted observations for a single value of p (e.g., p = 0.6)
w <- rbinom(1e4, size = 9, prob = 0.6)

simplehist(w)

# 3.26 ----
# propogate parameter uncertainty -> replace 0.6 with samples from posterior
w <- rbinom(1e4, size = 9, prob = samples)
