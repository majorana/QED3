#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif
#include <complex.h>
#include "lattice.h"
#include "linalg.h"
#include "rand/ranlxd.h"
#include "fields.h"

complex double g_fermion[GRIDPOINTS]; //Pseudofermion field
complex double g_eta[GRIDPOINTS];
complex double g_R[GRIDPOINTS];       //g_X = (gamma_5 D)^{-2} g_fermion
complex double g_temp1[GRIDPOINTS], g_temp2[GRIDPOINTS], g_temp3[GRIDPOINTS];

double At[GRIDPOINTS], Ax[GRIDPOINTS], Ay[GRIDPOINTS];         //Non-compact real-valued gauge fields
double At_old[GRIDPOINTS], Ax_old[GRIDPOINTS], Ay_old[GRIDPOINTS]; //Old values of the gauge field
complex double Ut[GRIDPOINTS], Ux[GRIDPOINTS], Uy[GRIDPOINTS];   //Compact lattice gauge fields: link = exp(I*gauge)
double gpt[GRIDPOINTS], gpx[GRIDPOINTS], gpy[GRIDPOINTS];               //Momenta corresponding to the gauge fields

double s_g, s_g_old;

double S_G(int i)
{
	return (-beta0/dt*(cos(At[i] + Ax[tp[i]] - At[xp[i]] - Ax[i]) 
			+ cos(At[i] + Ay[tp[i]] - At[yp[i]] - Ay[i])) 
			- beta*dt*cos(Ax[i] + Ay[xp[i]] - Ax[yp[i]] - Ay[i]) );
}

double DS_Gt(int i)
{
	return beta0/dt*(
			sin(At[i]+Ax[tp[i]] - At[xp[i]] - Ax[i]) -
			sin(At[xm[i]]+Ax[tp[xm[i]]] - At[i] - Ax[xm[i]]) +
			sin(At[i] + Ay[tp[i]] - At[yp[i]] - Ay[i]) - 
			sin(At[ym[i]]+Ay[tp[ym[i]]] - At[i] - Ay[ym[i]]) 
			);
}

double DS_Gx(int i)
{
	return beta*dt*(
			sin(Ax[i]+Ay[xp[i]] - Ax[yp[i]] - Ay[i]) -
			sin(Ax[ym[i]]+Ay[xp[ym[i]]] - Ax[i] - Ay[ym[i]]) 
			) + beta0/dt*(-sin(At[i] + Ax[tp[i]] - At[xp[i]] - Ax[i]) + 
			sin(At[tm[i]] + Ax[i] - At[xp[tm[i]]]  - Ax[tm[i]]) 
			);
}

double DS_Gy(int i)
{
	return beta*dt*(
			-sin(Ax[i]+Ay[xp[i]] - Ax[yp[i]] - Ay[i]) +
			sin(Ax[xm[i]]+ Ay[i] - Ax[yp[xm[i]]] - Ay[xm[i]]) 
			) + beta0/dt*(-sin(At[i] + Ay[tp[i]] - At[yp[i]] - Ay[i]) + 
			sin(At[tm[i]] + Ay[i] - At[yp[tm[i]]]  - Ay[tm[i]]) 
			);

}


void coldstart()
{
 	int i;
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		At[i]     = 0.0;
  		At_old[i] = 0.0;
  		Ax[i]     = 0.0;
  		Ax_old[i] = 0.0;
		Ay[i]     = 0.0;
  		Ay_old[i] = 0.0;

 	};
 	calculatelinkvars();
 	s_g=0;
 	for(i=0; i<GRIDPOINTS; i++)
  		s_g += S_G(i);
 	s_g_old = s_g;
}

void hotstart()
{
 	int i;
 	double r[GRIDPOINTS*3];
 	ranlxd(r, GRIDPOINTS*3);
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		At[i]=2*M_PI*r[i]-M_PI;
  		Ax[i]=2*M_PI*r[i+GRIDPOINTS]-M_PI;
		Ay[i]=2*M_PI*r[i+2*GRIDPOINTS]-M_PI;
		//Ax[i] = 0.0;
		//Ay[i] = 0.0;
 	}
 	calculatelinkvars();
 	s_g=0;
 	for(i=0; i<GRIDPOINTS; i++)
  		s_g += S_G(i);
 	s_g_old=s_g;
}

void calculatelinkvars()
{
 	int i;
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		Ut[i] = cos(At[i]) + I*sin(At[i]);
  		Ux[i] = cos(Ax[i]) + I*sin(Ax[i]);
		Uy[i] = cos(Ay[i]) + I*sin(Ay[i]);
 	};
}
