#ifndef MAIN_H
#define MAIN_H


#include "../constants.h"

#ifdef FLOW

    double *rho, *u, *v;

#endif

#if defined(FLOW) && !defined(DUALCOMPONENT)

    double *Fx, *Fy;

    double *f1, *f2;

#endif

#ifdef DUALCOMPONENT

    double *rho_lava, *rho_air;

    double *Fx_lava, *Fy_lava, *Fx_air, *Fy_air;

    double *f1_lava, *f2_lava, *f1_air, *f2_air;

#endif

#ifdef TEMPERATURE

    double *T, *g1, *g2;

#endif

#ifdef PHASECHANGE

    double *phi, *phi_old;

#endif

#endif