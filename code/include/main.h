#ifndef MAIN_H
#define MAIN_H


// --- DECLARE ARRAYS ---
#ifdef FLOW
    double *rho, *u, *v;
    #ifdef DUALCOMPONENT
        double *rho_lava, *rho_air;
        double *Fx_lava, *Fy_lava, *Fx_air, *Fy_air;
        double *f1_lava, *f2_lava, *f1_air, *f2_air;
    #else
        double *Fx, *Fy;
        double *f1, *f2;
    #endif
#endif

#ifdef TEMPERATURE
    double *T, *g1, *g2;
    #ifdef PHASECHANGE
        double *phi, *phi_old;
    #endif
#endif

#endif