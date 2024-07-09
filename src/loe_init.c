#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void fmajo(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getorder(void *, void *, void *);
extern void LOEgrad(void *, void *, void *, void *, void *, void *, void *);
extern void SDMLOE(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void SOEgrad(void *, void *, void *, void *, void *, void *, void *, void *);
extern void SOEobjt(void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP LOEobjt(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"fmajo",    (DL_FUNC) &fmajo,     9},
    {"getorder", (DL_FUNC) &getorder,  3},
    {"LOEgrad",  (DL_FUNC) &LOEgrad,   7},
    {"SDMLOE",   (DL_FUNC) &SDMLOE,   11},
    {"SOEgrad",  (DL_FUNC) &SOEgrad,   8},
    {"SOEobjt",  (DL_FUNC) &SOEobjt,   6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"LOEobjt", (DL_FUNC) &LOEobjt, 3},
    {NULL, NULL, 0}
};

void R_init_loe(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
