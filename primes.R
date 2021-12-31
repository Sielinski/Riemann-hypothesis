library(pracma)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)

#########################
## prime number theorem
#########################

# prime number theorem
# Legendre, Gauss, Enke
PNT <- function(x) {
  x / log(x)
}

# Dirichlet approximation
# Using Eulerian logarithmic integral
# https://en.wikipedia.org/wiki/Logarithmic_integral_function
Li <- function(x) {
  #integral(function(x) (1 / log(x)), 2, x)
  li(x) - li(2)
}

# Prime counting function
# From *Prime Obsession*, John Derbyshire,  p 297
PCF <- function(x) {
  p_x <- primes(x)
  length(p_x) - (0.5 * (tail(p_x, 1) == x))
}

# Reimann's prime counting function J(x)
# From *Prime Obsession*, p 299 (Expression 19-1)
J <- function(x) {
  i <- 1        # ith term
  j <- 0
  repeat {
    y <- x ^ (1 / i)  # ith root of x
    if (y < 2) break
    j <- j + PCF(y) / i
    i <- i + 1
  }
  j
}

# Mobius inversion
mobius_mu <- function(x) {
  if (x == 1)
    return(1)
  
  f <- factors(x)
  
  if (max(table(f)) >= 2)
    return(0)
  if (length(f) %% 2 == 0)
    return(1)
  else
    return(-1)
}

# Reimann's counting function J(x)
# in terms of mobius function
# https://medium.com/cantors-paradise/the-riemann-hypothesis-explained-fa01c1f75d3f
J_mu <- function(x) {
  i <- 1        # ith term
  j <- 0
  repeat {
    y <- x ^ (1 / i)  # ith root of x
    if (y < 2) break
    # Li version subtracts log(2)
    j <- j + Li(y) * mobius_mu(i) / i
    #j <- j + li(y) * mobius_mu(i) / i
    i <- i + 1
  }
  j
}

i <- 1000

PNT(i)
Li(i)
J(i)
J_mu(i)

PCF(i)

PNT(i) - PCF(i) # Legendre/Gauss
li(i) - PCF(i) # Dirichlet
J(i) - PCF(i)
J_mu(i) - PCF(i) # Reimann


# Visualize

prime_cnt <- data.frame(x = seq(2, 1000, 1))
prime_cnt$actual_cnt <-
  prime_cnt$x %>% map(primes) %>% map_dbl(length)
prime_cnt$pcf <- map_dbl(prime_cnt$x, PCF)
prime_cnt$Li <- map_dbl(prime_cnt$x, Li)
prime_cnt$li <- li(prime_cnt$x)
prime_cnt$pnt <- PNT(prime_cnt$x)
prime_cnt$j <- map_dbl(prime_cnt$x, J)
prime_cnt$j_mu <- map_dbl(prime_cnt$x, J_mu)


ggplot(prime_cnt[1:50, ], aes(x = x)) +
  geom_step(aes(y = actual_cnt)) +
  #geom_line(aes(y = li), color = 'purple') +
  geom_line(aes(y = Li), color = 'black') +
  geom_line(aes(y = pnt), color = 'red') +
  #geom_line(aes(y = j), color = 'blue') +
  geom_step(aes(y = pcf), color = 'orange') +
  geom_line(aes(y = j_mu), color = 'green')


##########################
## Riemann's hypothesis
##########################

my_zeta <- function(s, n = 10 ^ 6) {
  # this function works when the real part greater than 1
  if (s <= 1) return(NA)
  z  <- 0
  for (i in 1:n) {
    z <- z + i ^ -s
  }
  s=(complex(r=0.5, i=1*t));
  return(z)
}

# *Prime Obsession*, p 146
# this function works when the real part != 1
# Using Dirichlet eta function
eta_zeta <- function(s) {
  eta(s) / (1 - (1 / (2 ^ (s - 1))))
}


###############
# zeta zeroes
###############

# from https://en.wikipedia.org/wiki/Gamma_function
complex_gamma <- function(z) {
  integral(function(x) (x ^ (z - 1) * exp(-x)), 0, Inf)
}

# from https://www.ams.org/journals/tran/1988-309-02/S0002-9947-1988-0961614-2/S0002-9947-1988-0961614-2.pdf
# p 799
theta <- function(t) {
  Arg(pi ^ (1i * -t/2) * complex_gamma(0.25 + (1i * t/2)))
}

theta_v <- Vectorize(theta)

# t is the complex component of zeta
Z <- function(t) {
  exp(1i * theta_v(t)) * zeta(0.5 + (1i * t))
}

# from http://www.dtc.umn.edu/~odlyzko/doc/arch/zeta.zero.spacing.pdf
S <- function(t) {
  (pi ^ -1) * Arg(zeta(0.5 + (1i * t)))
}


# This is reasonable for the first ~30 zeroes (t ~ 102)
# Explore why: Perhaps the precision degrades as values get larger
critical_line <- data.frame(t = seq(10, 102, 0.001)) %>%
  mutate(s = S(t)) #%>%
#mutate(z = Z(t))

plt <- ggplot(critical_line, aes(x = t, y = s)) +
  geom_line()

my_zeroes <- critical_line %>%
  filter(sign(s) == 1 & sign(s) != sign(lag(s))) %>%
  select(t)

plt +
  geom_point(data = my_zeroes, aes(x = t, y = 0), color = 'red')


#
# For more accurate values, read the zeroes
#

# from http://www.dtc.umn.edu/~odlyzko/zeta_tables/index.html
zeroes <- read.delim('zeros_100k', header = F) %>%
  as.data.frame()
colnames(zeroes)[1] <- 't'

cbind(my_zeroes, zeroes$t[1:dim(my_zeroes)[1]])


# Look at some examples
dat <- data.frame(input = complex(real = 0.5, imaginary = seq(0, 35, .1))) %>%
  mutate(zeta = eta_zeta(input))

# plot the result of the zeta function for Re(0.5) on the complex plane
ggplot(dat, aes(x = Re(zeta), y = Im(zeta))) +
  geom_path() +
  geom_point(x = 0, y = 0, color = 'red')

# plot the real and imaginary parts of the zeta function for Re(0.5)
# with the Im() part of the critical line on x-axis
ggplot(dat, aes(x = Im(input))) +
  geom_path(aes(y = Re(zeta))) +
  geom_path(aes(y = Im(zeta)), color = 'blue') +
  geom_point(data = filter(zeroes, t <= max(Im(dat$input))), aes(x = t, y = 0), color = 'red') +
  scale_y_continuous(
    sec.axis = dup_axis(name = 'Im(zeta)')
  ) +
  theme(axis.title.y.right = element_text(color = 'blue'))


##############
# Error term
##############

# See Prime Obsession, p 328
# y is the number of non-trivial zeroes to include
#secondary_terms <- function(x, y = 10) {
#  x ^ (0.5 + 1i * zeroes$t[1:y]) %>%
#    li() %>%
#    Re() * 2
#}

# See "An Improved Analytic Method for Calculating π(x)", Jan Buthe
# https://arxiv.org/pdf/1410.7008.pdf p 2
secondary_terms <- function(x, y = 10) {
  expint_Ei((0.5 + 1i * zeroes$t[1:y]) * log(x)) %>%
    Re() * 2
}

# see *Prime Obsession*, p 340 - 341
x <- 20
secondary_terms(x, 50) %>% sum()
secondary_terms(x, 50) %>% plot(type = 'l')

integral_term <- function(x) {
  quadinf(function(t)
    1 / (t * (t ^ 2 - 1) * log(t)), xa = x, xb = Inf)$Q
}

# See *Riemann's Zeta Function*, H.M. Edwards, p 34, formula (3)
# see *Prime Obsession*, p 344
# See https://mathworld.wolfram.com/RiemannPrimeCountingFunction.html
J_error <- function(x, depth = 10) {
  i <- 1        # ith term
  j <- 0
  repeat {
    y <- x ^ (1 / i)  # ith root of x
    if (y < 2)
      break
    j <-
      j + mobius_mu(i) / i * (li(y) - sum(secondary_terms(y, depth)) - log(2) - integral_term(y))
    i <- i + 1
  }
  j
}

# Take a look at the results
prime_cnt$j_mu <- map_dbl(prime_cnt$x, J_mu)
prime_cnt$j_error <- map_dbl(prime_cnt$x, J_error, 500)

prime_cnt[1:100, ] %>% 
  select(x, actual_cnt, Li, pnt, j_mu, j_error) %>% 
  pivot_longer(cols = -'x', names_to = 'counting_function', values_to = 'primes') %>% 
  ggplot(aes(x = x, y = primes, color = counting_function)) + 
  geom_line() +
  #geom_step(aes(y = actual_cnt), color = 'blue', alpha = 0.5) +
  theme(legend.position="right") +
  labs(title = 'Prime Counting Functions', y = 'Number of Primes')


# Code that I explored trying to calculate the error term
# Not really needed

# See *Prime Obsession*, p 334
x ^ (0.5 + 1i * zeroes$t[1:20]) %>%
  plot(xlim = c(-5, 5), ylim = c(-5, 5))


#See http://www.chemistrylearning.com/logarithm-of-complex-number/
log_complex <- function(z) {
  a <- Re(z)
  b <- Im(z)
  
  r <- sqrt(a ^ 2 + b ^ 2)
  theta <- atan(b / a)
  log(r) + 1i * theta
}


# compare various methods
i <- complex(real = 0.5, imaginary = zeroes$t[1])
li(i)
expint_Ei(log(i))
expint_Ei(log_complex(i))


# https://mathworld.wolfram.com/Euler-MascheroniConstant.html
# https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant
euler_gamma <- 0.577215665


#####################
## Goldbach’s comet
#####################

# from https://stackoverflow.com/questions/19767408/prime-number-function-in-r
is_prime <- function(n) {
  if (n %in% 0:1) return(FALSE)
  n == 2L || all(n %% 2L:max(2, floor(sqrt(n))) != 0)
}

is_primes <- Vectorize(is_prime)

which(is_primes(1:10))


# from https://www.r-bloggers.com/playing-with-primes-in-r-part-ii/
# conjecture that every even number greater than 2 can be expressed
# as the sum of two primes.

gbp <- function(x) {
  i <- 1
  parts <- NULL
  while (i <= x / 2) {
    if (is_prime(i) && is_prime(x - i)) {
      if (i > 1) {
        parts <- append(parts, sprintf("%d+%d", i, x - i))
      }
    }
    i <- ifelse(x == 4, i + 1, i + 2)
  }
  return(parts)
}

gbp(26)


# How many distinct combinations of primes exist for each even number?
prime_cnt <- prime_cnt %>%
  #select(x) %>%
  mutate(gpb = map(x, gbp)) %>%
  mutate(gpb_l = map_dbl(gpb, length))

prime_cnt %>%
  filter(gpb_l != 0) %>%
  ggplot(aes(x, gpb_l)) +
  geom_point()


#
# Any even number >2 can be expressed as the sum of two primes
#

x <- 64
primes_all <- primes(x)
primes_cnt <- length(primes_all)
primes_tbl <- matrix(nrow = primes_cnt, ncol = primes_cnt)
colnames(primes_tbl) <- primes_all
rownames(primes_tbl) <- primes_all

# calculate the sum of each possible pair of primes
for (i in 1:primes_cnt) {
  for (j in i:primes_cnt) { # to avoid duplication, start at i (not 1)
    primes_tbl[i, j] <- primes_all[i] + primes_all[j]
  }
}

# find the pairs that add up to x
matches <- which(primes_tbl == x)
# pick a match from the middle (because the primes are more interesting)
match_idx <- matches[length(matches) / 2]

for (i in seq_along(matches)) {
  #print(matches[i])
  match_idx <- matches[i]
  # get the column/row index of the match from the table
  match_col <-
    ifelse(mod(match_idx, primes_cnt) == 0, primes_cnt, mod(match_idx, primes_cnt))
  match_row <- ceiling(match_idx / primes_cnt)
  
  # get the prime values
  two_primes <- c(primes_all[match_col], primes_all[match_row])
  
  print(two_primes)
  # add them up to ensure they = x
  # sum(two_primes) == x
}
