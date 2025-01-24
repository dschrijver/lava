#ifndef GLOBALS_H
#define GLOBALS_H

#include "../constants.h"

// --- DECLARE ARRAYS ---
#ifdef FLOW
    extern double *rho, *u, *v;
    #ifdef DUALCOMPONENT
        extern double *rho_lava, *rho_air;
        extern double *Fx_lava, *Fy_lava, *Fx_air, *Fy_air;
        extern double *f1_lava, *f2_lava, *f1_air, *f2_air;
    #else
        extern double *Fx, *Fy;
        extern double *f1, *f2;
    #endif
#endif

#ifdef TEMPERATURE
    extern double *T, *g1, *g2;
    #ifdef PHASECHANGE
        extern double *phi, *phi_old;
    #endif
#endif

#endif