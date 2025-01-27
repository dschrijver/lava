#ifndef GLOBALS_H
#define GLOBALS_H


#include "../constants.h"

#ifdef FLOW

    extern double *rho, *u, *v;

#endif

#if defined(FLOW) && !defined(DUALCOMPONENT)

    extern double *Fx, *Fy;

    extern double *f1, *f2;

#endif

#ifdef DUALCOMPONENT

    extern double *rho_lava, *rho_air;

    extern double *Fx_lava, *Fy_lava, *Fx_air, *Fy_air;

    extern double *f1_lava, *f2_lava, *f1_air, *f2_air;

#endif

#ifdef TEMPERATURE

    extern double *T, *g1, *g2;

#endif

#ifdef PHASECHANGE

    extern double *phi, *phi_old;

#endif

#endif