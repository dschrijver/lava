#include "../constants.h"
#include "../include/globals.h"
#include "../include/forcing.h"



void calculate_body_forces() {

    double Fx_i, Fy_i;

#ifdef PHASECHANGE

    double phi_i, phi2, rho_i;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            Fx_i = 0.0;
            Fy_i = 0.0;

#ifdef FLOW

            Fx_i += Fx_body;
            Fy_i += Fy_body;

    #ifdef TEMPERATURE
            Fy_i += alpha*(T[INDEX_2D(i,j)] - T_top)*g;
    #endif

    #ifdef PHASECHANGE
            phi_i = phi[INDEX_2D(i,j)];
            phi2 = phi_i*phi_i;
            rho_i = rho[INDEX_2D(i,j)];
            Fx_i += -(1.0-phi2)*rho_i*u[INDEX_2D(i,j)];
            Fy_i += -(1.0-phi2)*rho_i*v[INDEX_2D(i,j)];
    #endif

#endif

            Fx[INDEX_2D(i,j)] = Fx_i;
            Fy[INDEX_2D(i,j)] = Fy_i;

        }
    }

}