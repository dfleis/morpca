#include "utils.h"

double sign(double x) {
  if(x > 0.00000000001) return 1.0;
  else if(x < -0.00000000001) return -1.0;
  else return 0.0;
}

// [[Rcpp::export]]
SEXP hello_world2(SEXP x) {
  Rprintf("Hello, World!2");
  return List::create(x);
}
