library(data.table)
library(dplyr)
library(emdbook)

##### Check projection distribution #####

# Generate sample from two-component normal location model
n_obs <- 1000000
d <- 2
p <- 0.5
sigma <- sqrt(3)
true_sample <- matrix(NA, nrow = n_obs, ncol = d)

for(i in 1:n_obs) {
  mixture_comp <- rbinom(n = 1, size = 1, prob = p)
  if(mixture_comp == 0) {
    true_sample[i, ] <- rnorm(n = d, mean = 0, sd = 1)
  } else if(mixture_comp == 1) {
    true_sample[i, ] <- rnorm(n = d, mean = 0, sd = sigma)
  }
}

# Get random vector for projection
random_vector <- rnorm(n = d, mean = 0, sd = 1)
random_vector <- random_vector / sqrt(sum(random_vector^2))

# Projection of data
true_sample_proj <- true_sample %*% random_vector

# Empirical CDF
q <- -1
mean(true_sample_proj <= q)

# Normal mixture CDF
0.25 * pnorm(q, mean = 0, sd = 1) +
  0.25 * pnorm(q, mean = 0, sd = sqrt(random_vector[1]^2 + 3*random_vector[2]^2)) +
  0.25 * pnorm(q, mean = 0, sd = sqrt(3*random_vector[1]^2 + random_vector[2]^2)) +
  0.25 * pnorm(q, mean = 0, sd = sqrt(3))

##### Projection distribution #####

# Create projection density
f <- function(y, random_vector) {

  sigma_1 <- 1
  sigma_2 <- sqrt(random_vector[1]^2 + 3*random_vector[2]^2)
  sigma_3 <- sqrt(3*random_vector[1]^2 + random_vector[2]^2)
  sigma_4 <- sqrt(3)

  (1/4) * 
    ((1 / (sqrt(2*pi) * sigma_1)) * exp(-0.5 * y^2 / sigma_1^2) + 
     (1 / (sqrt(2*pi) * sigma_2)) * exp(-0.5 * y^2 / sigma_2^2) +
     (1 / (sqrt(2*pi) * sigma_3)) * exp(-0.5 * y^2 / sigma_3^2) +
     (1 / (sqrt(2*pi) * sigma_4)) * exp(-0.5 * y^2 / sigma_4^2))

}

# Check projection density
f(y = 1, random_vector = random_vector)

0.25 * dnorm(x = 1, mean = 0, sd = 1) +
  0.25 * dnorm(x = 1, mean = 0, sd = sqrt(random_vector[1]^2 + 3*random_vector[2]^2)) +
  0.25 * dnorm(x = 1, mean = 0, sd = sqrt(3*random_vector[1]^2 + random_vector[2]^2)) +
  0.25 * dnorm(x = 1, mean = 0, sd = sqrt(3))

##### First derivative of f wrt y #####

sigma_1 <- 1
sigma_2 <- sqrt(random_vector[1]^2 + 3*random_vector[2]^2)
sigma_3 <- sqrt(3*random_vector[1]^2 + random_vector[2]^2)
sigma_4 <- sqrt(3)

df_dy <- function(y, random_vector_1, random_vector_2) eval(D(expression(
  (1/4) * 
    ((1 / (sqrt(2*pi) * 1)) * exp(-0.5 * y^2 / 1) + 
     (1 / (sqrt(2*pi) * sqrt(random_vector_1^2 + 3*random_vector_2^2))) * 
       exp(-0.5 * y^2 / (random_vector_1^2 + 3*random_vector_2^2)) +
     (1 / (sqrt(2*pi) * sqrt(3*random_vector_1^2 + random_vector_2^2))) * 
       exp(-0.5 * y^2 / (3*random_vector_1^2 + random_vector_2^2)) +
     (1 / (sqrt(2*pi) * sqrt(3))) * exp(-0.5 * y^2 / 3))
), "y"))

df_dy_written <- function(y, random_vector){

  sigma_1 <- 1
  sigma_2 <- sqrt(random_vector[1]^2 + 3*random_vector[2]^2)
  sigma_3 <- sqrt(3*random_vector[1]^2 + random_vector[2]^2)
  sigma_4 <- sqrt(3)

  - (y / (4* sqrt(2*pi) * sigma_1^3)) * exp(-0.5 * y^2 / sigma_1^2) -
    (y / (4* sqrt(2*pi) * sigma_2^3)) * exp(-0.5 * y^2 / sigma_2^2) -
    (y / (4* sqrt(2*pi) * sigma_3^3)) * exp(-0.5 * y^2 / sigma_3^2) -
    (y / (4* sqrt(2*pi) * sigma_4^3)) * exp(-0.5 * y^2 / sigma_4^2)
}

df_dy(y = -1, random_vector_1 = random_vector[1], random_vector_2 = random_vector[2])

df_dy_written(y = -1, random_vector = random_vector)

##### Second derivative of log(f) wrt y #####

d2logf_dyy <- function(y, random_vector_1, random_vector_2) eval(D(D(expression(
  log((1/4) * 
    ((1 / (sqrt(2*pi) * 1)) * exp(-0.5 * y^2 / 1) + 
     (1 / (sqrt(2*pi) * sqrt(random_vector_1^2 + 3*random_vector_2^2))) * 
       exp(-0.5 * y^2 / (random_vector_1^2 + 3*random_vector_2^2)) +
     (1 / (sqrt(2*pi) * sqrt(3*random_vector_1^2 + random_vector_2^2))) * 
       exp(-0.5 * y^2 / (3*random_vector_1^2 + random_vector_2^2)) +
     (1 / (sqrt(2*pi) * sqrt(3))) * exp(-0.5 * y^2 / 3)))
), "y"), "y"))

d2logf_dyy_written <- function(y, random_vector) {

  sigma_1 <- 1
  sigma_2 <- sqrt(random_vector[1]^2 + 3*random_vector[2]^2)
  sigma_3 <- sqrt(3*random_vector[1]^2 + random_vector[2]^2)
  sigma_4 <- sqrt(3)

  term_11 <- (-1/(32*pi*sigma_1^4)) * exp(-y^2/sigma_1^2)

  term_12 <- (y^2/(32*pi*sigma_1*sigma_2^5) - y^2/(32*pi*sigma_1^3*sigma_2^3) - 
    1/(32*pi*sigma_1*sigma_2^3)) * exp(-y^2/(2*sigma_1^2) - y^2/(2*sigma_2^2))

  term_13 <- (y^2/(32*pi*sigma_1*sigma_3^5) - y^2/(32*pi*sigma_1^3*sigma_3^3) - 
    1/(32*pi*sigma_1*sigma_3^3)) * exp(-y^2/(2*sigma_1^2) - y^2/(2*sigma_3^2))

  term_14 <- (y^2/(32*pi*sigma_1*sigma_4^5) - y^2/(32*pi*sigma_1^3*sigma_4^3) - 
    1/(32*pi*sigma_1*sigma_4^3)) * exp(-y^2/(2*sigma_1^2) - y^2/(2*sigma_4^2))

  term_21 <- (y^2/(32*pi*sigma_2*sigma_1^5) - y^2/(32*pi*sigma_2^3*sigma_1^3) - 
    1/(32*pi*sigma_2*sigma_1^3)) * exp(-y^2/(2*sigma_2^2) - y^2/(2*sigma_1^2))

  term_22 <- (-1/(32*pi*sigma_2^4)) * exp(-y^2/sigma_2^2)

  term_23 <- (y^2/(32*pi*sigma_2*sigma_3^5) - y^2/(32*pi*sigma_2^3*sigma_3^3) - 
    1/(32*pi*sigma_2*sigma_3^3)) * exp(-y^2/(2*sigma_2^2) - y^2/(2*sigma_3^2))

  term_24 <- (y^2/(32*pi*sigma_2*sigma_4^5) - y^2/(32*pi*sigma_2^3*sigma_4^3) - 
    1/(32*pi*sigma_2*sigma_4^3)) * exp(-y^2/(2*sigma_2^2) - y^2/(2*sigma_4^2))

  term_31 <- (y^2/(32*pi*sigma_3*sigma_1^5) - y^2/(32*pi*sigma_3^3*sigma_1^3) - 
    1/(32*pi*sigma_3*sigma_1^3)) * exp(-y^2/(2*sigma_3^2) - y^2/(2*sigma_1^2))

  term_32 <- (y^2/(32*pi*sigma_3*sigma_2^5) - y^2/(32*pi*sigma_3^3*sigma_2^3) - 
    1/(32*pi*sigma_3*sigma_2^3)) * exp(-y^2/(2*sigma_3^2) - y^2/(2*sigma_2^2))

  term_33 <- (-1/(32*pi*sigma_3^4)) * exp(-y^2/sigma_3^2)

  term_34 <- (y^2/(32*pi*sigma_3*sigma_4^5) - y^2/(32*pi*sigma_3^3*sigma_4^3) - 
    1/(32*pi*sigma_3*sigma_4^3)) * exp(-y^2/(2*sigma_3^2) - y^2/(2*sigma_4^2))

  term_41 <- (y^2/(32*pi*sigma_4*sigma_1^5) - y^2/(32*pi*sigma_4^3*sigma_1^3) - 
    1/(32*pi*sigma_4*sigma_1^3)) * exp(-y^2/(2*sigma_4^2) - y^2/(2*sigma_1^2))

  term_42 <- (y^2/(32*pi*sigma_4*sigma_2^5) - y^2/(32*pi*sigma_4^3*sigma_2^3) - 
    1/(32*pi*sigma_4*sigma_2^3)) * exp(-y^2/(2*sigma_4^2) - y^2/(2*sigma_2^2))

  term_43 <- (y^2/(32*pi*sigma_4*sigma_3^5) - y^2/(32*pi*sigma_4^3*sigma_3^3) - 
    1/(32*pi*sigma_4*sigma_3^3)) * exp(-y^2/(2*sigma_4^2) - y^2/(2*sigma_3^2))

  term_44 <- (-1/(32*pi*sigma_4^4)) * exp(-y^2/sigma_4^2)

  numerator <- term_11 + term_12 + term_13 + term_14 + 
    term_21 + term_22 + term_23 + term_24 +
    term_31 + term_32 + term_33 + term_34 + 
    term_41 + term_42 + term_43 + term_44

  denominator <- (f(y = y, random_vector = random_vector))^2

  numerator / denominator
}

# These do not agree!
d2logf_dyy(y = 1, random_vector_1 = random_vector[1], random_vector_2 = random_vector[2])

d2logf_dyy_written(y = 1, random_vector = random_vector)

dlogf_dx1_written <- function(x1, x2){
  (1/f(x1, x2)) * (-x1/(4*pi)) *
    (exp(-(x1^2 + x2^2)/2) + (1/9)*exp(-(x1^2 + x2^2)/6))
}

dlogf_dx1(x1 = 2, x2 = 1.5)

dlogf_dx1_written(x1 = 2, x2 = 1.5)

##### Second derivative wrt (x1, x1) #####

d2logf_dx1x1 <- function(x1, x2) eval(D(D(expression(log( (1 / (4 * pi)) *
  (exp(-(x1^2 + x2^2)/2) + (1/3) * exp(-(x1^2 + x2^2) / 6)) )), "x1"), "x1"))

d2logf_dx1x1_written <- function(x1, x2){
  (-1/(4*pi*f(x1, x2))) * (exp(-(x1^2 + x2^2)/2) + (1/9)*exp(-(x1^2 + x2^2)/6)) -
    (x1/(4*pi*f(x1, x2)))^2 * (exp(-(x1^2 + x2^2)/2) + (1/9)*exp(-(x1^2 + x2^2)/6))^2 +
    (x1^2 / (4*pi*f(x1,x2))) * exp(-(x1^2 + x2^2)/2) +
    (x1^2 / (108*pi*f(x1,x2))) * exp(-(x1^2 + x2^2)/6)
}

d2logf_dx1x1(x1 = 5, x2 = 3)

d2logf_dx1x1_written(x1 = 5, x2 = 3)

##### Second derivative wrt (x2, x2) #####

d2logf_dx2x2 <- function(x1, x2) eval(D(D(expression(log( (1 / (4 * pi)) *
  (exp(-(x1^2 + x2^2)/2) + (1/3) * exp(-(x1^2 + x2^2) / 6)) )), "x2"), "x2"))

d2logf_dx2x2_written <- function(x1, x2){
  (-1/(4*pi*f(x1, x2))) * (exp(-(x1^2 + x2^2)/2) + (1/9)*exp(-(x1^2 + x2^2)/6)) -
    (x2/(4*pi*f(x1, x2)))^2 * (exp(-(x1^2 + x2^2)/2) + (1/9)*exp(-(x1^2 + x2^2)/6))^2 +
    (x2^2 / (4*pi*f(x1,x2))) * exp(-(x1^2 + x2^2)/2) +
    (x2^2 / (108*pi*f(x1,x2))) * exp(-(x1^2 + x2^2)/6)
}

d2logf_dx2x2(x1 = 5, x2 = 3)

d2logf_dx2x2_written(x1 = 5, x2 = 3)

##### Second derivative wrt (x1, x2) #####

d2logf_dx1x2 <- function(x1, x2) eval(D(D(expression(log( (1 / (4 * pi)) *
  (exp(-(x1^2 + x2^2)/2) + (1/3) * exp(-(x1^2 + x2^2) / 6)) )), "x1"), "x2"))

d2logf_dx1x2_written <- function(x1, x2){
  -(x1*x2/(f(x1, x2)*4*pi)^2) * (exp(-(x1^2 + x2^2)/2) + (1/9)*exp(-(x1^2 + x2^2)/6))^2 +
    (x1*x2/(4*pi*f(x1, x2))) * exp(-(x1^2 + x2^2)/2) +
    (x1*x2/(108*pi*f(x1, x2))) * exp(-(x1^2 + x2^2)/6)
}

d2logf_dx1x2(x1 = 1.5, x2 = 5)

d2logf_dx1x2_written(x1 = 1.5, x2 = 5)

##### Hessian matrix #####

zt_H_z <- function(x1, x2, z1, z2) {

  H <- matrix(c(d2logf_dx1x1(x1 = x1, x2 = x2),
                d2logf_dx1x2(x1 = x1, x2 = x2),
                d2logf_dx1x2(x1 = x1, x2 = x2),
                d2logf_dx2x2(x1 = x1, x2 = x2)),
              nrow = 2)

  z <- matrix(c(z1, z2), ncol = 1)

  return(t(z) %*% H %*% z)
}

##### Determinant of Hessian matrix #####

det_H <- function(x1, x2) {

  H <- matrix(c(d2logf_dx1x1(x1 = x1, x2 = x2),
                d2logf_dx1x2(x1 = x1, x2 = x2),
                d2logf_dx1x2(x1 = x1, x2 = x2),
                d2logf_dx2x2(x1 = x1, x2 = x2)),
              nrow = 2)

  return(det(H))
}

##### Search z^T H z for choice of x, z such that this is > 0 #####

search_matrix <-
  expand.grid(x1_val = seq(0, 5, by = 0.2), x2_val = seq(0, 5, by = 0.2),
              z1_val = seq(-5, 5, by = 0.2), z2_val = seq(-5, 5, by = 0.2))

search_matrix <- as.data.table(cbind(search_matrix, value = NA_real_))

search_matrix[, value :=
                zt_H_z(x1 = x1_val, x2 = x2_val,
                       z1 = z1_val, z2 = z2_val), by = seq_len(nrow(search_matrix))]

# Log of density is not concave at point (0, 3)
zt_H_z(x1 = 0, x2 = 3, z1 = 0, z2 = 1)

##### Search z^T H z for choices of x, z such that this is > 0, #####
##### zoomed in around x = (0, 3) #####

search_matrix_zoom <-
  expand.grid(x1_val = seq(-0.5, 0.5, by = 0.1), x2_val = seq(2.5, 3.5, by = 0.1),
              z1_val = seq(-5, 5, by = 1), z2_val = seq(-5, 5, by = 1))

search_matrix_zoom <- as.data.table(cbind(search_matrix_zoom, value = NA_real_))

search_matrix_zoom[, value :=
                     zt_H_z(x1 = x1_val, x2 = x2_val,
                            z1 = z1_val, z2 = z2_val), by = seq_len(nrow(search_matrix_zoom))]

# x1 \in [-0.5, 0.5] and x2 \in [2.5, 3.0] seems to be a non-LC region
search_matrix_zoom %>%
  dplyr::filter(value > 0) %>%
  dplyr::group_by(x1_val, x2_val) %>%
  slice(1) %>%
  dplyr::arrange(x1_val, x2_val) %>%
  View()

##### Search determinant of Hessian for convex points (det < 0) #####

search_determinant_zoom <- expand.grid(x1_val = seq(-0.5, 0.5, by = 0.1), 
                                       x2_val = seq(2.5, 3.5, by = 0.1))

search_determinant_zoom <- as.data.table(cbind(search_determinant_zoom, 
                                               value = NA_real_))

search_determinant_zoom[, value := det_H(x1 = x1_val, x2 = x2_val), 
                        by = seq_len(nrow(search_determinant_zoom))]

search_determinant_zoom %>%
  dplyr::filter(value < 0) %>%
  dplyr::group_by(x1_val, x2_val) %>%
  slice(1) %>%
  dplyr::arrange(x1_val, x2_val) %>% View

##### Plot function for x1 \in [-0.5, 0.5] and x2 \in [2.5, 3.0] #####

curve3d(expr = log(f(x1 = x, x2 = y)), from = c(-0.5, 2.5), to = c(0.5, 3))

##### Plot function at (2, 2) #####

curve3d(expr = log(f(x1 = x, x2 = y)), from = c(1.9, 1.9), to = c(2.1, 2.1))

##### Plot function (zoomed out) #####
curve3d(expr = log(f(x1 = x, x2 = y)), from = c(-10, -10), to = c(10, 10))

