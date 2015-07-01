#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif
#include "fields.h"
#include "complex/complex.h"
#include "lattice.h"
#include "linalg.h"
#include "rand/ranlxd.h"
#include "fermion.h"

void fermion(complex double *out, complex double *in) 
{
	int i;
  	for(i=0; i<GRIDPOINTS; i++) {
		out[i] = exp(-g_mu)*Ut[i]*in[tp[i]] - in[i]
			- g_t*(Ux[i]*in[xp[i]] + cconj(Ux[xm[i]])*in[xm[i]] + Uy[i]*in[yp[i]] + cconj(Uy[ym[i]])*in[ym[i]]);
	}
	return;
}

void fermion_fp(complex double *out, complex double *temp, complex double *in)
{
	fermion(out, in);
	return;
}

void fermion_herm(complex double *out, complex double *in) 
{
	int i;
  	for(i=0; i<GRIDPOINTS; i++) {
		out[i] = exp(-g_mu)*cconj(Ut[tm[i]])*in[tm[i]] - in[i]
			- g_t*(cconj(Ux[xm[i]])*in[xm[i]] + cconj(Ux[ym[i]])*in[ym[i]] + Ux[i]*in[xp[i]] + Uy[i]*in[yp[i]]);
	}
	return;
}

void fermion_sqr(complex double *out, complex double *temp, complex double *in)
{
	fermion_herm(temp, in);
	fermion(out, temp);
	return;
}

void fermion_DGx(complex double *out, complex double *in, int j)
{
	set_zero(out);
	out[j] = -I*g_t*Ux[j]*in[xp[j]];
	out[xp[j]] = I*g_t*cconj(Ux[j])*in[j];
	return;
}

void fermion_DGy(complex double *out, complex double *in, int j)
{
	set_zero(out);
	out[j] = -I*g_t*Uy[j]*in[yp[j]];
	out[yp[j]] = I*g_t*cconj(Uy[j])*in[j];
	return;
}

void fermion_DGt(complex double *out, complex double *in, int j)
{
	set_zero(out);
	out[j] = I*Ut[j]*in[tp[j]];
	return;
}


