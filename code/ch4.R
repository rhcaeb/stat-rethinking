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

## 4.6 - bayes theorem -> grid approximation calculation:
w <- 6; n <- 9;
p_grid <- seq(from = 0, to = 1, length.out = 100)
posterior <- dbinom(w, n, p_grid) * dunif(p_grid, 0, 1)
posterior <- posterior/sum(posterior)

## 4.7 - Howell (census data for the Dobe area !Kung San)
# library(rethinking)
data("Howell1")
d <- Howell1

## 4.8 - Inspect the structure of the data
str(d)

## 4.9 - rethinking summary function
precis(d, hist=FALSE)

## 4.10 - Access the height vector ($ = extract)
d$height

## 4.11 - filter non-adults (>= 18)
d2 <- d[d$age >= 18, ]

## 4.12 
curve(dnorm(x, 178, 20), from = 100, to = 250)

## 4.13 - std dev is a flat/uniform prior to constrain std dev to have posititive 
# plausibility between 0 - 50 cm:
curve(dunif(x, 0, 50), from = -10, to = 60)

## 4.14 - simulate heights from the prior
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

## 4.15 - E.g. prior with large sigma
sample_mu <- rnorm(1e4, mean = 178, sd = 100)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

## 4.16 - Grid approximation of the posterior distribution
mu.list <- seq(from = 150, to = 160, length.out = 100)
sigma.list <- seq(from = 7, to = 9, length.out = 100)
post <- expand.grid(mu = mu.list, sigma = sigma.list)
post$LL <- sapply(1:nrow(post), function(i) sum(
  dnorm(d2$height, post$mu[i], post$sigma[i], log = TRUE)))
post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) +
  dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))

## 4.17 - contour plot
contour_xyz(post$mu, post$sigma, post$prob)

## 4.18 - or plot with a simple heat map
image_xyz(post$mu, post$sigma, post$prob)

## 4.19 - sampling from the posterior (randomly sample rows in prop'n to post$prob)
sample.rows <- sample(1:nrow(post), size = 1e4, replace = TRUE, prob = post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]

## 4.20 - view the 1e4 samples (w replacement) from the height data
plot(sample.mu, sample.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.1))

## 4.21 - characterize the shapes of the marginal posterior densitives of mu and sigma
dens(sample.mu)
dens(sample.sigma)

## 4.22 - summarise the widths of these densities with compatibility intervals 
PI(sample.mu)
PI(sample.sigma)


# Sample size and the normality of sigmas posterior
## 4.23 - sample 20 random heights from the original list
d3 <- sample(d2$height, size = 20)

## 4.24 - repeat code from previous subsections modified to 20 samples
mu.list <- seq(from = 150, to = 170, length.out = 200)
sigma.list <- seq(from = 4, to = 20, length.out = 200)
post2 <- expand.grid(mu = mu.list, sigma = sigma.list)
post2$LL <- sapply(1:nrow(post2), function(i)
  sum(dnorm(d3, mean=post2$mu[i], sd = post2$sigma[i], log = TRUE)))
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )

## 4.25 - inspect the marginal posterior density for sigma, averaing over mu
dens(sample2.sigma, norm.comp = TRUE)

## 4.26 - Find the posterior distribution with quap()
data("Howell1")
d <- Howell1
d2 <- d[d$age >= 18, ] # filter to include adults

## 4.27 - formula list for quap()
flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)

## 4.28 - Fit the model to the data ('d2')
m4.1 <- quap(flist, data = d2)

## 4.29 - Posterior distribution
precis(m4.1)

## 4.31 - Changing sigma to a narrow prior
m4.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0, 50)
  ), data = d2)
precis(m4.2)

## 4.32 - View matrix of variances and covariances
vcov(m4.1)

## 4.33 - Vector of variances & Correlation matrix
diag(vcov(m4.1)) # list of variances (square root of this = sigma)
cov2cor(vcov(m4.1)) # correlation matrix

## 4.34 - Sample vectors of values from multi-dimensional Gaussian dist.
post <- extract.samples(m4.1, n = 1e4)
head(post)
precis(post)
