#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
#include "../include/forcing.h"



void calculate_body_forces() {

#ifdef DUALCOMPONENT

    int x_i, y_i;

    double rho_i, rho_lava_i, rho_air_i;

    double Fx_lava_i, Fy_lava_i, Fx_air_i, Fy_air_i;

    #ifdef TEMPERATURE
        double F_buoy;
    #endif

#else

    double Fx_i, Fy_i;

    #ifdef PHASECHANGE
        double rho_i;
    #endif

#endif

#ifdef PHASECHANGE

    double phi_i, phi2;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

#ifdef DUALCOMPONENT

            rho_i = rho[INDEX_2D(i,j)];
            rho_lava_i = rho_lava[INDEX_2D(i,j)];
            rho_air_i = rho_air[INDEX_2D(i,j)];

            Fx_lava_i = 0.0; 
            Fy_lava_i = 0.0; 
            Fx_air_i = 0.0;
            Fy_air_i = 0.0;

            for (int p = 0; p < NP; p++) {

                x_i = i+cx_i[p];
                y_i = j+cy_i[p];

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

            Fx_lava_i = -G*rho_lava_i*Fx_lava_i;
            Fy_lava_i = -G*rho_lava_i*Fy_lava_i;

            Fx_air_i = -G*rho_air_i*Fx_air_i;
            Fy_air_i = -G*rho_air_i*Fy_air_i; 

            Fx_lava_i += rho_lava_i/rho_i*Fx_body;
            Fy_lava_i += rho_lava_i/rho_i*Fy_body;

            Fx_air_i += rho_air_i/rho_i*Fx_body;
            Fy_air_i += rho_air_i/rho_i*Fy_body; 


    #ifdef TEMPERATURE
            F_buoy = (rho_0_lava*alpha_lava*(T[INDEX_2D(i,j)] - T_0_lava) + 
                      rho_0_air*alpha_air*(T[INDEX_2D(i,j)] - T_0_air))*g;
            Fy_lava_i += rho_lava_i/rho_i*F_buoy;
            Fy_air_i += rho_air_i/rho_i*F_buoy;
    #endif

    #ifdef PHASECHANGE
            phi_i = phi[INDEX_2D(i,j)];
            phi2 = phi_i*phi_i;
            Fx_lava_i += -(1.0-phi2)*rho_lava_i*u[INDEX_2D(i,j)];
            Fy_lava_i += -(1.0-phi2)*rho_lava_i*v[INDEX_2D(i,j)];
    #endif

            Fx_lava[INDEX_2D(i,j)] = Fx_lava_i;
            Fy_lava[INDEX_2D(i,j)] = Fy_lava_i;

            Fx_air[INDEX_2D(i,j)] = Fx_air_i;
            Fy_air[INDEX_2D(i,j)] = Fy_air_i;
            
#else

            Fx_i = Fx_body;
            Fy_i = Fy_body;

    #ifdef TEMPERATURE
            Fy_i += rho_0*alpha*(T[INDEX_2D(i,j)] - T_0)*g;
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

        }
    }

}