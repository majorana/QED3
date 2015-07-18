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
	int i,j;
  	for(i=0; i<GRIDPOINTS; i++) {
		out[i] = in[i] -exp(-g_mu*dt)*Ut[i]*in[tp[i]] 
			+ dt*Ut[i]*(Ux[tp[i]]*in[tp[xp[i]]] + cconj(Ux[tp[xm[i]]])*in[tp[xm[i]]] + Uy[tp[i]]*in[tp[yp[i]]] + cconj(Uy[tp[ym[i]]])*in[tp[ym[i]]]);
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
		out[i] =  in[i] - exp(-g_mu*dt)*cconj(Ut[tm[i]])*in[tm[i]]
			+ dt*(cconj(Ut[tm[xm[i]]]*Ux[xm[i]])*in[tm[xm[i]]] + cconj(Ut[tm[ym[i]]]*Uy[ym[i]])*in[tm[ym[i]]] + cconj(Ut[tm[xp[i]]])*Ux[i]*in[tm[xp[i]]] + cconj(Ut[tm[yp[i]]])*Uy[i]*in[tm[yp[i]]]);
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
	out[tm[j]] = I*Ux[j]*in[xp[j]];
	out[tm[xp[j]]] = -I*cconj(Ux[j])*in[j];
	return;
}


void fermion_DGy(complex double *out, complex double *in, int j)
{
	set_zero(out);
	out[tm[j]] = I*Uy[j]*in[yp[j]];
	out[tm[yp[j]]] = -I*cconj(Uy[j])*in[j];
	return;
}

void fermion_DGt(complex double *out, complex double *in, int j)
{
	set_zero(out);
	out[j] = -I*exp(-g_mu*dt)*Ut[j]*(in[tp[j]] + dt*(Ux[tp[j]]*in[tp[xp[j]]] + Uy[tp[j]]*in[tp[yp[j]]] + cconj(Ux[tp[xm[j]]])*in[tp[xm[j]]] + cconj(Uy[tp[ym[j]]])*in[tp[ym[j]]] ) );
	return;
}


// Calculate fermion force. Assume g_eta is the inverse of MM^\dag
// eta^\dag*dM/dA*(M^\dag)*eta
double fermion_forcet(const int i)
{
	double f;
	fermion_herm(g_temp1, g_eta);
	fermion_DGt(g_temp2, g_temp1, i);
	f = 2*scalar_prod_r(g_eta, g_temp2);
	return(f);
}

double fermion_forcex(const int i)
{
	double f;
	fermion_herm(g_temp1, g_eta);
	fermion_DGx(g_temp2, g_temp1, i);
	f = 2*scalar_prod_r(g_eta, g_temp2);
	return(f);
}

double fermion_forcey(const int i)
{
	double f;
	fermion_herm(g_temp1, g_eta);
	fermion_DGy(g_temp2, g_temp1, i);
	f = 2*scalar_prod_r(g_eta, g_temp2);
	return(f);
}


