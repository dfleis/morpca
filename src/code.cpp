#include "utils.h"

RcppExport SEXP hello_world(SEXP x) {
  Rprintf("Hello, World!");
  return List::create(x);
}
