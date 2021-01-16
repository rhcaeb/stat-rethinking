# Chapter 2 ----

## 2.1 - The Garden of Forking Data 

ways <- c(0, 3, 8, 9, 0) # ways to produce data
ways/sum(ways) # plausibility

## 2.3 - Components Of the Model

dbinom(x = 6, size = 9, prob = 0.5) # binomial distribution;
# likelihood of data: '6' water in '9' tosses (N), under any value of p.

## 2.4 - Making the Model 'Go'