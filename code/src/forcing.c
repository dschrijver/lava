#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
#include "../include/forcing.h"



void calculate_body_forces() {


#ifndef DUALCOMPONENT

    double Fx_i, Fy_i;

#endif

#ifdef DUALCOMPONENT

    int x_i, y_i;

    double rho_i, rho_lava_i, rho_air_i;

    double Fx_lava_i = 0.0, Fy_lava_i = 0.0, Fx_air_i = 0.0, Fy_air_i = 0.0;

#endif

#if defined(PHASECHANGE) && !defined(DUALCOMPONENT)

    double phi_i, phi2, rho_i;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

#ifndef DUALCOMPONENT

            Fx_i = Fx_body;
            Fy_i = Fy_body;

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

            Fx[INDEX_2D(i,j)] = Fx_i;
            Fy[INDEX_2D(i,j)] = Fy_i;

#endif

#ifdef DUALCOMPONENT

            rho_i = rho[INDEX_2D(i,j)];
            rho_lava_i = rho_lava[INDEX_2D(i,j)];
            rho_air_i = rho_air[INDEX_2D(i,j)];

            for (int p = 0; p < NP; p++) {

                x_i = i-cx_i[p];
                y_i = j-cy_i[p];

    #ifdef XPERIODIC
                if (x_i < 0) x_i = NX-1;
                else if (x_i == NX) x_i = 0;
    #endif

    #ifdef SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i < 0) {
                    x_i = i;
                    y_i = j;
                }
    #endif

    #ifdef NORTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i == NY) {
                    x_i = i;
                    y_i = j;
                }
    #endif

                Fx_lava_i += w[p]*rho_air[INDEX_2D(x_i, y_i)]*cx[p];
                Fy_lava_i += w[p]*rho_air[INDEX_2D(x_i, y_i)]*cy[p];

                Fx_air_i += w[p]*rho_lava[INDEX_2D(x_i, y_i)]*cx[p];
                Fy_air_i += w[p]*rho_lava[INDEX_2D(x_i, y_i)]*cy[p];

            }

            Fx_lava_i = -G*rho_lava_i*Fx_lava_i + rho_lava_i/rho_i*Fx_body;
            Fy_lava_i = -G*rho_lava_i*Fy_lava_i + rho_lava_i/rho_i*Fy_body;

            Fx_air_i = -G*rho_air_i*Fx_air_i + rho_air_i/rho_i*Fx_body;
            Fy_air_i = -G*rho_air_i*Fy_air_i + rho_air_i/rho_i*Fy_body; 

            Fx_lava[INDEX_2D(i,j)] = Fx_lava_i;
            Fy_lava[INDEX_2D(i,j)] = Fy_lava_i;

            Fx_air[INDEX_2D(i,j)] = Fx_air_i;
            Fy_air[INDEX_2D(i,j)] = Fy_air_i;
            
#endif

        }
    }

}