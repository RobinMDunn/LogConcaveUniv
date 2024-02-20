library(data.table)
library(dplyr)
library(emdbook)

f <- function(x1, x2) {
  (1 / (4 * pi)) *
    (exp(-(x1^2 + x2^2)/2) + (1/3) * exp(-(x1^2 + x2^2) / 6))
}

##### First derivative wrt x1 #####

dlogf_dx1 <- function(x1, x2) eval(D(expression(log( (1 / (4 * pi)) *
  (exp(-(x1^2 + x2^2)/2) + (1/3) * exp(-(x1^2 + x2^2) / 6)) )), "x1"))

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

