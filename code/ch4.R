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

## Chapter 4.5 - Curves from lines ----

## 4.64 - polynomial regression
library(rethinking)
data(Howell1)
d <- Howell1
str(d)

plot(height ~ weight, data = d)

## 4.65 - pre-process variable transformations
d$weight_s <- (d$weight - mean(d$weight))/sd(d$weight)
d$weight_s2 <- d$weight_s^2

m4.5 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    sigma ~ dunif(0, 50)
  ), data = d
)

## 4.67 - calculate mean relationship & 89% intervals of the mean and preds:
weight.seq <- seq(from = -2.2, to = 2, length.out = 30)
pred_dat <- list(weight_s = weight.seq, weight_s2 = weight.seq^2)
mu <- link(m4.5, data = pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.89)
sim.height <- sim(m4.5, data = pred_dat)
height.PI <- apply(sim.height, 2, PI, prob = 0.89)

## 4.68 - plot above
plot(height ~ weight_s, data = d, col = col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

## 4.69 - fit model with mod (cubic regression on weight) of parabolic code
d$weight_s3 <- d$weight_s^3 # new var

m4.6 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,10),
    b3 ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ), data = d
)

## Chapter 4.5.2 - Splines ----

## 4.72 - Japanese cherry blossom data ('wiggly data')
library(rethinking)
data("cherry_blossoms")
d <- cherry_blossoms
precis(d)

plot(temp ~ year, data = d)

## 4.73 - define a list of knot positions
d2 <- d[complete.cases(d$temp), ] # complete cases on temp
num_knots <- 15
knot_list <- quantile(d2$year, probs = seq(0, 1, length.out = num_knots))

## 4.74 - construct basis functions for degree 3, cubic, spline:
library(splines)
B <- bs(d2$year,
        knots = knot_list[(-c(1, num_knots))],
        degree = 3, intercept = TRUE)

## 4.75 - display basis functions -> plot columns against years
plot( NULL , xlim=range(d2$year) , ylim=c(0,1) , xlab="year" , ylab="basis value")
for ( i in 1:ncol(B) ) lines( d2$year , B[,i] )

## 4.76 - build model in quap()
m4.7 <- quap(
  alist(
    T ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(6, 10),
    w ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data = list(T = d2$temp, B=B), start = list(w = rep(0, ncol(B)))
)

precis(m4.7, depth = 2)

## 4.77 - need to plot the posterior predictions
# first -> weighted functions
post <- extract.samples(m4.7)
w <- apply(post$w, 2, mean)
plot( NULL , xlim=range(d2$year) , ylim=c(-2,2) ,
      xlab="year" , ylab="basis * weight" )
for ( i in 1:ncol(B) ) lines( d2$year , w[i]*B[,i] )

## 4.78 - knots added for reference
mu <- link(m4.7)
mu_PI <- apply(mu, 2, PI, 0.97)
plot( d2$year , d2$temp , col=col.alpha(rangi2,0.3) , pch=16 )
shade( mu_PI , d2$year , col=col.alpha("black",0.5) )
