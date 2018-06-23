library(microbenchmark)
library(Rcpp)
sourceCpp("src/helpers.cpp")


Y <- matrix(rnorm(5 * 3), nrow = 5)
r <- 2

mb <- microbenchmark(rank_r_approx_cpp(Y, 2),
                     init_approx(Y, 2), times = 1e2, unit = "ms")
boxplot(mb, outline = F)
