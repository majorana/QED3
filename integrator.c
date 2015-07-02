#include <stdlib.h>
#include <stdio.h>
#include "hmc.h"
#include "integrator.h"
#include "fermion.h"
#include "fields.h"
#include "test.h"

/*  leap frog integrator */
void leapfrog(const double dtau) {
  	/* first phase: \Delta\Tau / 2 step for p */
  	update_gauge(0.5*dtau); 


  	update_momenta(dtau);

  	/*  last phase: \Delta\Tau / 2 step for p */
  	update_gauge(dtau*0.5);

}

void integrator(const int steps, const double stepsize) {
	int i, j;
	double ham1, sg1;

	for(i = 0; i < steps; i++)
	{
		sg1 = 0.0;
		ham1 = 0.0;
 		for (j=0; j<GRIDPOINTS; j++)
 		{
  			sg1 += S_G(j);
  			ham1 += 0.5*(gpx[j]*gpx[j] + gpy[j]*gpy[j] + gpt[j]*gpt[j]);
 		};
 		ham1 += sg1;
	
		g_cgiterations1 += cg(g_temp3, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
		ham1 += scalar_prod_r(g_fermion, g_temp3);
		printf("Hamiltonian: %f\n", ham1);

		leapfrog(stepsize);
	}
}

void update_momenta_gauge(const double dtau) 
{
  	int i;
  	for(i = 0; i < GRIDPOINTS; i++) {
		gpt[i] = gpt[i] - dtau*DS_Gt(i);
    	gpx[i] = gpx[i] - dtau*DS_Gx(i);
    	gpy[i] = gpy[i] - dtau*DS_Gy(i);
  	}
  	return;
}


void update_momenta(const double dtau) 
{
  	int i;
	double fft[GRIDPOINTS];
	double ffx[GRIDPOINTS];
	double ffy[GRIDPOINTS];
	double gfx[GRIDPOINTS];

  	g_cgiterations1 += cg(g_eta, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
	//print_vector(g_eta);
  	for(i = 0; i < GRIDPOINTS; i++) {
		//printf("%d \n", i);
		//printf("Brute-force %f\n", stupid_fermion_force_x(i));
		//printf("Smart %f\n", fermion_forcex(i));
		gpt[i] = gpt[i] - dtau*(DS_Gt(i) - fermion_forcet(i));
    	gpx[i] = gpx[i] - dtau*(DS_Gx(i) - fermion_forcex(i));
    	gpy[i] = gpy[i] - dtau*(DS_Gy(i) - fermion_forcey(i));
#ifdef MONITOR_MD
		ffx[i] = stupid_fermion_force_x(i);
		gfx[i] = DS_Gx(i);
		fft[i] = fermion_forcet(i);
		ffy[i] = fermion_forcey(i);
#endif
  	}
#ifdef MONITOR_MD
	//print_vector_r(ffx);
	//print_vector_r(fft);
	printf("Max x fermion force: %.3f\n", max_r(ffx));
	printf("Max y fermion force: %.3f\n", max_r(ffy));
	printf("Max t fermion force: %.3f\n", max_r(fft));
	if (max_r(ffx) > 1000.0 || max_r(fft) > 1000.0) {
		//print_fermion_mat();
	}
#endif
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
