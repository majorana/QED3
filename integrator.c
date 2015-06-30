#include <stdlib.h>
#include "hmc.h"
#include "integrator.h"
#include "fermion.h"
#include "fields.h"

/*  leap frog integrator */
void integrator(const double dtau) {
  int l;

  	/* first phase: \Delta\Tau / 2 step for p */
  	update_momenta(0.5*dtau); 

  	update_gauge(dtau);

  	/*  last phase: \Delta\Tau / 2 step for p */
  	update_momenta(dtau*0.5);
}

// Calculate fermion force. Assume g_eta is the inverse of MM^\dag
double fermion_forcet(const int i)
{
	double f;
	fermion_herm(g_temp1, g_eta);
	fermion_DGt(g_temp2, g_temp1, i);
	f = 2*scalar_prod_r(g_eta, g_temp2);
	return(f);
}

// eta^\dag*dM/dA*(M^\dag)*eta
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

void update_momenta(const double dtau) 
{
  	int i;
  	g_cgiterations1 += cg(g_eta, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
  	for(i = 0; i < GRIDPOINTS; i++) {
		gpt[i] = gpt[i] - dtau*(DS_Gt(i) - fermion_forcet(i));
    	gpx[i] = gpx[i] - dtau*(DS_Gx(i) - fermion_forcex(i));
    	gpy[i] = gpy[i] - dtau*(DS_Gy(i) - fermion_forcey(i));
  	}
  	return;
}

void update_gauge(const double dtau) {
  int i;
  for(i = 0; i < GRIDPOINTS; i++) {
    Ax[i] = Ax[i] + dtau*gpx[i];
    Ay[i] = Ay[i] + dtau*gpy[i];
	At[i] = At[i] + dtau*gpt[i];
  }
  calculatelinkvars();
  return;
}
