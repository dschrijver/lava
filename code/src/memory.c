#include <stdlib.h> // malloc(), free()

#include "../constants.h"
#include "../include/globals.h"
#include "../include/memory.h"


void allocate_memory() {

#ifdef FLOW

    rho = (double*) malloc(NX*NY*sizeof(double));
    u = (double*) malloc(NX*NY*sizeof(double));
    v = (double*) malloc(NX*NY*sizeof(double));

#endif

#ifdef DUALCOMPONENT

    rho_lava = (double*) malloc(NX*NY*sizeof(double));
    rho_air = (double*) malloc(NX*NY*sizeof(double));

    Fx_lava = (double*) malloc(NX*NY*sizeof(double));
    Fy_lava = (double*) malloc(NX*NY*sizeof(double));
    Fx_air = (double*) malloc(NX*NY*sizeof(double));
    Fy_air = (double*) malloc(NX*NY*sizeof(double));

    f1_lava = (double*) malloc(NX*NY*NP*sizeof(double));
    f2_lava = (double*) malloc(NX*NY*NP*sizeof(double));
    f1_air = (double*) malloc(NX*NY*NP*sizeof(double));
    f2_air = (double*) malloc(NX*NY*NP*sizeof(double));

#else

    Fx = (double*) malloc(NX*NY*sizeof(double));
    Fy = (double*) malloc(NX*NY*sizeof(double));

    f1 = (double*) malloc(NX*NY*NP*sizeof(double));
    f2 = (double*) malloc(NX*NY*NP*sizeof(double));

#endif

#ifdef TEMPERATURE

    T = (double*) malloc(NX*NY*sizeof(double));

    g1 = (double*) malloc(NX*NY*NP*sizeof(double));
    g2 = (double*) malloc(NX*NY*NP*sizeof(double));

#endif

#ifdef PHASECHANGE

    phi = (double*) malloc(NX*NY*sizeof(double));
    phi_old = (double*) malloc(NX*NY*sizeof(double));

#endif

}


void free_memory() {

#ifdef FLOW

    free(rho);
    free(u);
    free(v);

#endif

#ifdef DUALCOMPONENT

    free(rho_lava);
    free(rho_air);

    free(Fx_lava);
    free(Fy_lava);
    free(Fx_air);
    free(Fy_air);

    free(f1_lava);
    free(f2_lava);
    free(f1_air);
    free(f2_air);

#else

    free(Fx);
    free(Fy);

    free(f1);
    free(f2);


#endif

#ifdef TEMPERATURE

    free(T);

    free(g1);
    free(g2);

#endif

#ifdef PHASECHANGE

    free(phi);
    free(phi_old);

#endif

}