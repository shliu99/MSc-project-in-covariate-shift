# supporting functions for KLIEP
# run these functions when a kliep object is ready

compute_kernel_Gaussian <- function(x, centers, sigma) {
  apply(centers, 1, function(center) {
    apply(x, 1, function(row) {
      kernel_Gaussian(row, center, sigma)
    })
  })
}

kernel_Gaussian <- function(x, y, sigma) {
  exp(- squared_euclid_distance(x, y) / (2 * sigma * sigma))
}

squared_euclid_distance <- function(x, y) {
  sum((x - y) ^ 2)
}
