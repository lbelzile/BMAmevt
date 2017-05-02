#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void d_dirimix(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void d_pairbeta(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void d_pairbeta_grid(void *, void *, void *, void *, void *, void *, void *);
extern void d_trinestlog(void *, void *, void *, void *, void *, void *, void *);
extern void d_trinestlog_grid(void *, void *, void *, void *, void *, void *, void *);
extern void ddirimix_grid(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ddirimix_grid1D(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"d_dirimix",         (DL_FUNC) &d_dirimix,         11},
    {"d_pairbeta",        (DL_FUNC) &d_pairbeta,         9},
    {"d_pairbeta_grid",   (DL_FUNC) &d_pairbeta_grid,    7},
    {"d_trinestlog",      (DL_FUNC) &d_trinestlog,       7},
    {"d_trinestlog_grid", (DL_FUNC) &d_trinestlog_grid,  7},
    {"ddirimix_grid",     (DL_FUNC) &ddirimix_grid,      9},
    {"ddirimix_grid1D",   (DL_FUNC) &ddirimix_grid1D,    7},
    {NULL, NULL, 0}
};

void R_init_BMAmevt(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
