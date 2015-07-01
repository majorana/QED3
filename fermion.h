#ifndef _FERMION_H
#define _FERMION_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator, 
 * its square and derivatives over the gauge fields ***/
/**********************************************************************************************/

#include "linalg.h"

extern void fermion(complex double *out, complex double *in);

extern void fermion_fp(complex double *out, complex double *temp, complex double *in);

extern void fermion_herm(complex double *out, complex double *in);

extern void fermion_sqr(complex double *out, complex double *temp, complex double *in);

extern void fermion_DGt(complex double *out, complex double *in, int j);

extern void fermion_DGx(complex double *out, complex double *in, int j);

extern void fermion_DGy(complex double *out, complex double *in, int j);

#endif
