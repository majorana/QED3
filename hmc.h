#ifndef _HMC_H
#define _HMC_H

/***********************************************************************************/
/*** This unit implement the basic HMC update step and the necessary procedures ****/
/***********************************************************************************/

/* Maximal number of iterations for CG method */
#define ITER_MAX 1000
/* Tolerance for CG method */
#define DELTACG 1.e-22
/* Step size and the number of steps for the leapfrog */
/* declared in qed.c */
extern int    g_steps;
extern double g_stepsize;

extern double ham, ham_old;

extern int R;  // Counter of all accepted configurations
extern int g_cgiterations, g_cgiterations1, g_cgiterations2;
extern int hmc_iter;

int  update(); //Basic HMC update step
int  accept(const double exphdiff); //Accept or reject the trajectory depending on exphdiff

#endif
