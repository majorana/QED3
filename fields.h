#include "linalg.h"  
#include "rand/ranlxd.h"


/* declared in qed.c */
extern double g_mu;
extern double g_t;
extern double beta0, beta;

extern complex double g_fermion[GRIDPOINTS]; //Pseudofermion field
extern complex double g_eta[GRIDPOINTS];
extern complex double g_R[GRIDPOINTS];       //g_X = (gamma_5 D)^{-2} g_fermion
extern complex double g_temp1[GRIDPOINTS], g_temp2[GRIDPOINTS];

extern double At[GRIDPOINTS], Ax[GRIDPOINTS], Ay[GRIDPOINTS];         //Non-compact real-valued gauge fields
extern double At_old[GRIDPOINTS], Ax_old[GRIDPOINTS], Ay_old[GRIDPOINTS]; //Old values of the gauge field
extern complex double Ut[GRIDPOINTS], U[GRIDPOINTS], U[GRIDPOINTS];   //Compact lattice gauge fields: link = exp(I*gauge)
extern double gpt[GRIDPOINTS], gpx[GRIDPOINTS], gpy[GRIDPOINTS];               //Momenta corresponding to the gauge fields

extern double s_g, s_g_old; //Action for the gauge fields in the beginning and in the end of the MD trajectory

void coldstart(); //Cold start for the gauge fields
void hotstart();  //Hot (random) start for the gauge fields

void calculatelinkvars(); //Sets the values of the compact gauge fields from noncompact ones

double S_G(int i);   //Action of the gauge fields
double DS_Gt(int i); //Derivative of the action over the time component of the noncompact gauge field
double DS_Gx(int i); //Derivative of the action over the x component of the noncompact gauge field
double DS_Gy(int i); //Derivative of the action over the y component of the noncompact gauge field
#endif
