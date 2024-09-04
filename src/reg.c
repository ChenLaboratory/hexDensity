#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(hbin  )(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
static const R_FortranMethodDef FortranEntries[] = {
	{"hbin", (DL_FUNC) &F77_NAME(hbin), 14},
	{NULL, NULL, 0}
};

SEXP meanderingTrianglesC(SEXP xleft,SEXP xright,SEXP y,SEXP z,SEXP levels);
static const R_CallMethodDef CallEntries[] = {
  {"meanderingTrianglesC", (DL_FUNC) &meanderingTrianglesC,5},
  {NULL, NULL, 0}
};


void R_init_hexDensity(DllInfo *info) {
	R_registerRoutines(info, NULL, CallEntries, FortranEntries, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}
