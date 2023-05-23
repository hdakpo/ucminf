#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP mfopt(SEXP rho);

static const R_CallMethodDef CallEntries[] = {
    {"mfopt", (DL_FUNC) &mfopt, 1},
    {NULL, NULL, 0}
};

void R_init_ucminf(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
