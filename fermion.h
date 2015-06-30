#ifndef _DIRAC_H
#define _DIRAC_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator, its square and derivatives over the gauge fields ***/
/**********************************************************************************************/

#include "linalg.h"

extern void fermion(complex double *out, complex double *in);

extern void fermion_fp(complex double *out, complex double *temp, complex double *in);

extern void fermion_herm(complex double *out, complex double *in);

extern void fermion_sqr(complex double *out, complex double *temp, complex double *in);

#endif
