#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

//This is the implementation of the Leapfrog integrator ... 
void integrator(const int steps, const double stepsize);

void omelyan(const double dtau);

void update_momenta(const double dtau);

void update_momenta_gauge(const double dtau);

void update_gauge(const double dtau);

#endif
