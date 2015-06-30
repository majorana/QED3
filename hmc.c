#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "complex/complex.h"
#include "rand/ranlxd.h"
#include "rand/gauss.h"
#include "linalg.h"
#include "fields.h"
#include "lattice.h"
#include "fermion.h"
#include "integrator.h"
#include "hmc.h"

int R;
int g_cgiterations;
int g_cgiterations1;
int g_cgiterations2;
double ham, ham_old;

double hamiltonian()
{
	ham = 0;
 	s_g = 0;
 	for (i=0; i<GRIDPOINTS; i++)
 	{
  		s_g += S_G(i);
  		ham += 0.5*(gpx[i]*gpx[i] + gpy[i]*gpy[i] + gpt[i]*gpt[i]);
 	};
 	ham += s_g;
	return(ham)
}

int update() //Basic HMC update step
{
 	int i, acc;
 	double exphdiff;
 
 	ham_old = 0.0;
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		gpx[i] = gauss();  //Momenta conjugate to the gauge fields 
  		gpy[i] = gauss();
		gpt[i] = gauss();
  		ham_old += 0.5*(gpx[i]*gpx[i] + gpy[i]*gpy[i] + gpt[i]*gpt[i]);
 	};
 	ham_old += s_g_old; //s_g_old contains the action of the gauge fields
 
	for(i=0; i<GRIDPOINTS; i++)
 	{
  		g_R[i] = (gauss() + I*gauss())/sqrt(2); //Pseudofermion fields times M^{-1} 
 	};
	squnrm = square_norm(g_R);
  	fermion(g_fermion, g_R); //g_fermion the pseudofermion field, i.e. phi = M R
  	ham_old += squnrm;

 	// Do the molecular dynamic chain
 	integrator(g_steps, g_stepsize);
 
 	// Calculate the new action and hamiltonian
 	ham = 0;
 	s_g = 0;
 	for (i=0; i<GRIDPOINTS; i++)
 	{
  		s_g += S_G(i);
  		ham += 0.5*(gpx[i]*gpx[i] + gpy[i]*gpy[i] + gpt[i]*gpt[i]);
 	};
 	ham += s_g;

	// Calculate phi (M M^\dag)^{-1} phi = R^\dag R, R= M^{-1}phi
	g_cgiterations1 += cg(g_R, g_fermion, ITER_MAX, DELTACG, &fermion_fp);
	ham += scalar_prod_r(g_R, g_R);

 	exphdiff = exp(ham_old-ham);
 	acc = accept(exphdiff);
 
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		At_old[i] = At[i];
  		Ax_old[i] = Ax[i];
		Ay_old[i] = Ay[i];
 	};
 	s_g_old = s_g;
 
 	//Increase the counter of the total number of MC updates...
 	hmc_iter++;
 	//Return 1 if the configuration was accepted, 0 otherwise
 	return(acc);
}

int accept(const double exphdiff)
{
  	int acc=0, i;
  	double r[1];

  	// the acceptance step
  	if(exphdiff>=1) {
    	acc = 1; 
    	R += 1;
  	}
  	else {
    	ranlxd(r,1);
    	if(r[0]<exphdiff) {
      		acc = 1;
      		R += 1;
    }
    else {
      // get the old values for phi, cause the configuration was not accepted
    	for (i=0; i<GRIDPOINTS; i++)
		{
	  		gauge1[i]=gauge1_old[i];
	  		gauge2[i]=gauge2_old[i];
		};
      	calculatelinkvars();
      	s_g = s_g_old;
	}
}
  return acc;
}



