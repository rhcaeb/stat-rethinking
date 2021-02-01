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

## 4.37 - Linear prediction: plot height ~ weight -> how do they covary?
plot(d2$height ~ d2$weight)

## 4.38 - Simulate lines
set.seed(2971) # reproduce text example
N <- 100 # 100 lines
a <- rnorm(N, 178, 20)
b <- rnorm(N, 0, 10)

# now (above, 4.38) we have 100 pairs of alpha and beta values

## 4.39 - Plot lines
plot(NULL, xlim=range(d2$weight), ylim = c(-100, 400),
     xlab = "weight", ylab = "height")
abline(h = 0, lty = 2)
abline(h = 272, lty = 1, lwd = 0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for(i in 1:N) curve(a[i] + b[i]*(x - xbar),
                    from = min(d2$weight), to = max(d2$weight), add = TRUE,
                    col = col.alpha("black", 0.2))

## 4.40 - simulate relationship and see what this means for beta (using log-normal)
b <- rlnorm(1e4, 0, 1)
dens(b, xlim = c(0,5), adj = 0.1)

## 4.41 - prior predictive sim. using log-normal prior
set.seed(2971) # reproduce text e.g.
N <- 100 # 100 lines
a <- rnorm(N, 178, 20)
b <- rlnorm(N, 0, 1)

plot(NULL, xlim=range(d2$weight), ylim = c(-100, 400),
     xlab = "weight", ylab = "height")
abline(h = 0, lty = 2)
abline(h = 272, lty = 1, lwd = 0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for(i in 1:N) curve(a[i] + b[i]*(x - xbar),
                    from = min(d2$weight), to = max(d2$weight), add = TRUE,
                    col = col.alpha("black", 0.2))


## 4.42 - Build the posterior approximation
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]

# define the average weight, x-bar
xbar <- mean(d2$weight)

# fit model
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data = d2)

## 4.44 - tables of marginal posterior distributions of params.
precis(m4.3)

## 4.45 - variance-covariance matrix
round(vcov(m4.3), 3)
pairs(m4.3) # marginal posterios and covariance

## 4.46 - plotting posterior inference against the data
plot(height ~ weight, data = d2, col=rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map*(x-xbar), add = TRUE)

## 4.47 - visualizing uncertainty in the regression relationship
post <- extract.samples(m4.3)
post[1:5, ] # each row is a correlated random sample
# paired values of alpha and beta define a line

## 4.48 - add more data to see the change in scatter of lines (re-estimate model)
N <- 10
dN <- d2[1:N, ] # extract first 10 cases
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - mean(weight)),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = dN)

## 4.49 - plot 20 of these lines to view uncertainty
# extract 20 samples from the posterior
post <- extract.samples(mN, n = 20)

# display raw data and sample size
plot(dN$weight, dN$height,
     xlim = range(d2$weight), ylim = range(d2$height),
     col = rangi2, xlab = "weight", ylab = "height")
mtext(concat("N = ", N))

# plot the lines with transparency
for(i in 1:20)
  curve(post$a[i] + post$b[i]*(x - mean(dN$weight)),
        col = col.alpha("black", 0.3), add = TRUE)

## 4.50 - compute arbitrary interval(s) using cloud of regression lines 
# E.g. 50 kilograms ($weight)

post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50 - xbar)
dens(mu_at_50, col = rangi2, lwd = 2, xlab = "mu|weight = 50")


## 4.51 - 89% compatibility interval of mu at 50 kg
PI(mu_at_50, prob = 0.89)

## 4.53 - repeat above for every weight value (not just 50 kg)
mu <- link(m4.3)
str(mu) # matrix of values of mu (each row = sample from post. dist.)

## 4.54 - distribution of mu for each weight(i)
# define sequence of weights to computer predictions for
# these values will be on the horizontal axis
weight.seq <- seq(from = 25, to = 70, by = 1)

# use link to compute mu
# for each sample from posterior
# and for each  weight.seq
mu <- link(m4.3, data = data.frame(weight = weight.seq))
str(mu)

## 4.55 - plot distribution of mu values at each height 
# use type = 'n' to hide raw data
plot(height ~ weight, d2, type = 'n')

# loop over samples and plot each mu value
for(i in 1:100)
  points(weight.seq, mu[i, ], pch = 16, col = col.alpha(rangi2, 0.1))

## 4.56 - summarise the distribution for each weight value
mu.mean <- apply(mu, 2, mean) # avg. mu at each weight value
mu.PI <- apply(mu, 2, PI, prob=0.89) # lower/upper bounds for each weight value

## 4.57 - plot summaries on top of the data
# plot raw data
# fading out points to make line and interval more visible
plot(height ~ weight, data = d2, col = col.alpha(rangi2, 0.5))

# plot the MAP line, aka the mean mu for each weight
lines(weight.seq, mu.mean)

# plot a shaded region for 89% PI
shade(mu.PI, weight.seq)

## 4.59 - Simulate heights (predicting sigma)
sim.height <- sim(m4.3, data = list(weight = weight.seq))
str(sim.height) # simulated heights

## 4.60 - summarize simulated heights
height.PI <- apply(sim.height, 2, PI, prob = 0.89) # 89% posterior prediction interval of observable heights
# across values of weight in weight.seq

## 4.61 - Plot!
# plot raw data
plot(height ~ weight, d2, col = col.alpha(rangi2, 0.5))

# draw MAP line
lines(weight.seq, mu.mean)

# draw HPDI region for line
shade(mu.HPDI, weight.seq)

# draw PI region for simulated heights
shade(height.PI, weight.seq)

