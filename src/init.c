#include <R.h>
#include <Rinternals.h> // for SEXP
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>  // optional

extern SEXP hello_world(SEXP x);

extern SEXP _morpca_hello_world2(SEXP x);

static R_CallMethodDef callMethods[] = {
  // {"cdfit_mgaussian_ssr", (DL_FUNC) &cdfit_mgaussian_ssr, 15},
  {"hello_world", (DL_FUNC) &hello_world, 1},
  {"_morpca_hello_world2", (DL_FUNC) &_morpca_hello_world2, 1},
  {NULL, NULL, 0}
};

void R_init_morpca(DllInfo *dll) {
  R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(dll, FALSE);
}
