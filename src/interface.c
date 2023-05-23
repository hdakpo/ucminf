/*
  Stig Bousgaard Mortensen
  DTU Informatics
  sbm@imm.dtu.dk

  Based on R-package 'FortranCallsR' by
  Diethelm Wuertz, ETH Zurich, www.rmetrics.org

  Modifications by Douglas Bates <bates@stat.wisc.edu>, Nov. 2010
  Modifications by Tomas Kalibera <tomas.kalibera@gmail.com> Aug. 2016.
*/


#include <R.h>
#include <Rinternals.h>  //R internal structures
#include <R_ext/RS.h>    //F77_CALL etc.

// Declare FORTRAN routine for use in C
extern void F77_NAME(ucminf)(int*, double[], double*, double[],
			     int*,double[],int*,int*,int*,double[],SEXP);

/*-------------------------------------------------------------------------------
  Define C functions that calls user defined function in R
*/

void installPar(int nn, double x[], SEXP rho) {
    int i;
    SEXP PAR = findVarInFrame(rho, install(".x"));
    double *xpt = REAL(PAR);
    if (LENGTH(PAR) != nn)
	error("Dimension mismatch, length(.x) = %d != n = $d", LENGTH(PAR), nn);
    for (i = 0; i < nn; i++) xpt[i] = x[i] ;
}

void F77_SUB(func)(int *n, double x[], double *value, SEXP rho) {
    installPar(*n, x, rho);
    SEXP dotf = findVarInFrame(rho, install(".f"));
    PROTECT(dotf);
    *value = asReal(PROTECT(eval(dotf, rho))) ;
    UNPROTECT(2);
}

void F77_SUB(usrgr)(int *n, double x[], double grval[], SEXP rho) {
    SEXP OUT;
    int i, nn = *n;
    double *grv;

    installPar(nn, x, rho);
    SEXP dotgr = PROTECT(findVarInFrame(rho, install(".gr")));
    PROTECT(OUT = eval(dotgr, rho));
    if (LENGTH(OUT) != nn || !isReal(OUT))
	error("gradient evaluation must return a numeric vector of length %d", nn);
    grv = REAL(OUT);
    for (i = 0; i < nn; i++) grval[i] = grv[i];
    UNPROTECT(2) ;
}

/*--------------------------------------------------------------------------------
  Define C function to be called from R
*/

SEXP mfopt(SEXP rho) {
    int n      = asInteger(PROTECT(findVarInFrame(rho, install(     ".n")))),
	iw     = asInteger(PROTECT(findVarInFrame(rho, install(    ".iw")))),
	grad   = asInteger(PROTECT(findVarInFrame(rho, install(  ".grad"))));
    UNPROTECT(3);
    SEXP
	EPS    = PROTECT(findVarInFrame(rho, install(   ".eps"))),
	GRSTEP = PROTECT(findVarInFrame(rho, install(".grstep"))),
	PAR    = PROTECT(findVarInFrame(rho, install(   ".par"))),
	icontr = PROTECT(findVarInFrame(rho, install(".icontr"))),
	maxfun = PROTECT(findVarInFrame(rho, install(".maxfun"))),
	dx     = PROTECT(findVarInFrame(rho, install(".stepmax"))),
	W      = PROTECT(findVarInFrame(rho, install(     ".w")));

    if (LENGTH(EPS) < 2 || !isReal(EPS))
	error(".eps must be a numeric vector of length >= 2");
    if (LENGTH(GRSTEP) < 2 || !isReal(GRSTEP))
	error(".eps must be a numeric vector of length >= 2");
    if (LENGTH(PAR) != n || !isReal(PAR))
	error("Dimension mismatch, length(.par) = %d != n = $d", LENGTH(PAR), n);
    if (LENGTH(W) != iw || !isReal(W))
	error("Dimension mismatch, length(.w) = %d != .iw = $d", LENGTH(W), iw);

    // duplicate dx, maxfun, .w because they are input/output arguments
    maxfun = PROTECT(duplicate(maxfun));
    defineVar(install(".maxfun"), maxfun, rho);
    dx = PROTECT(duplicate(dx));
    defineVar(install(".stepmax"), dx, rho);
    W = PROTECT(duplicate(W));
    defineVar(install(".w"), W, rho);
    UNPROTECT(3); // now protected via rho

				//  Call the FORTRAN routine 'ucminf'
    F77_CALL(ucminf)(&n, REAL(PAR), REAL(dx), REAL(EPS), INTEGER(maxfun), REAL(W), &iw, INTEGER(icontr),
		     &grad, REAL(GRSTEP), rho) ;
    UNPROTECT(7);
    return R_NilValue;
}

void F77_SUB(prtrac)(int *neval, double *fx, double *nmg, int *n, double x[]) {
    int i, nn = *n;
    Rprintf(" neval = %3d, F(x) =%11.4e, max|g(x)| =%11.4e\n", *neval, *fx, *nmg);
    Rprintf(" x =%11.4e", x[0]);
    for (i = 1; i < nn; i++) Rprintf(",%11.4e", x[i]);
    Rprintf("\n");
}

void F77_SUB(prline)(double *a, double sl[]) {
    Rprintf(" Line search: alpha =%11.4e, dphi(0) =%11.4e, dphi(1) =%11.4e\n",
	    *a, sl[0], sl[1]);
}

void F77_SUB(prconv)(void) {
    Rprintf(" Optimization has converged\n");
}

void F77_SUB(prfail)(int *neval) {
    Rprintf(" Optimization stopped after %d function evaluations\n", *neval);
}
