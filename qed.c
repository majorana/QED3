#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "lattice.h"
#include "hmc.h"
#include "fields.h"

/* global variables */
double g_mu = 1.0;
double g_t = 1.0;

int g_thermalize   = 100;   //Number of MC updates for thermalization
int g_measurements = 200;    //Number of measurements (statistically independent configurations)
int g_intermediate =  2;    //Number of MC updates between the measurements

/* extern in fields.h   */

double beta0  = 1.0;
double beta   = 1.0;        //Coupling constant for the gauge field

/* extern in hmc.h      */
int    g_steps    = 100;      //Number of steps in the molecular dynamics chain
double g_stepsize = 0.01;    //Size of each step

void echo_sim_params();

int main(int argc, char **argv) 
{
	int i, l;
  	int accepted = 0;        //Total number of accepted configurations
  	int total_updates = 0;   //Total number of updates
  	/* Initialize the random number generator */
  	rlxd_init(2, time(NULL)); 
  	/* Initialize the lattice geometry */
  	init_lattice(Lx, Ly, Lt);
  	/* Initialize the fields */
  	hotstart();
  	/* Print out the run parameters */
  	echo_sim_params();
 	
	printf("\n Test run:\n");
	update();
	return 0;
   
  	/* thermalization */
  	hmc_iter = 0; //Counts the total number of calls to the update() routine
  	printf("\n Thermalization: \n\n");
  	for(i=0; i<g_thermalize; i++)
  	{
   		update();
		//printf("\t Step %04i,\t mp = %2.4lf,\t pl = %2.4lf,\t cc = %2.4lf\n", i, mean_plaquette(), polyakov_loop());
  	};
  	
	/* measure the iterations only during real simulation, not thermalization */
  	R              = 0; //Counts the total number of accepted configurations
  	g_cgiterations = 0; //Counts the total number of CG iterations
  	hmc_iter       = 0; //Counts the total number of calls to the update() routine

  	printf("\n Generation: \n\n");
  	for(i=0; i<g_measurements; i++) 
  	{
   	/* do g_intermediate updates before measurement */
   		for (l=0; l<g_intermediate; l++)
   		{
    		accepted += update();
   		};
   		accepted += update();
   	/* Measurements should go here... */
  
  	};
 
  	/* Some output for diagnostics */
  	total_updates = g_measurements*(g_intermediate + 1);
  	printf("\n\n Algorithm performance:\n");
  	printf("\t Acceptance rate:             %2.2lf\n", (double)accepted/(double)total_updates);
  	printf("\t CG iterations per update:    %2.2lf\n", (double)g_cgiterations/(double)total_updates);

  	system("PAUSE");
  	return 0;
}

void echo_sim_params()
{
 	printf("Hybrid Monte-Carlo for U(1) gauge theory with spinor Fermi surface\n\n");
 	printf("Run parameters (adjust in qed.c !!!):\n");
 	printf("\t Beta:                            %2.4lf\n",  beta);
 	printf("\t Lattice size:                    %i x %i x %i\n", Lx, Ly, Lt);
 	printf("\t HMC step size:                   %2.4lf\n",  g_stepsize);
 	printf("\t HMC no. of steps:                %i\n",      g_steps);
 	printf("\t Thermalization steps:            %i\n",      g_thermalize);
 	printf("\t Number of measurements:          %i\n",      g_measurements);
 	printf("\t MC updates between measurements: %i\n",      g_intermediate);
 	printf("\n\n");
}
