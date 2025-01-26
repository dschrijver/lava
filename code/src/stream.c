#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
#include "../include/stream.h"


#ifdef FLOW

void stream_hydrodynamic_populations() {

    int x_i, y_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int p = 0; p < NP; p++) {

                x_i = i-cx_i[p];
                y_i = j-cy_i[p];

#ifdef XPERIODIC
                if (x_i < 0) x_i = NX-1;
                else if (x_i == NX) x_i = 0;
#endif

#ifdef WEST_PRESSURE_NEBB
                if (x_i < 0) continue;
#endif

#ifdef EAST_PRESSURE_NEBB
                if (x_i == NX) continue;
#endif

#ifdef SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i < 0) {
                    f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(i, j, p_bounceback[p])];
                    continue;
                }
#endif

#ifdef NORTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i == NY) {
                    f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(i, j, p_bounceback[p])];
                    continue;
                }
#endif

                f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(x_i, y_i, p)];
                
            }
        }
    }

#if defined WEST_PRESSURE_NEBB || defined EAST_PRESSURE_NEBB

    double u_i;

    for (int j = 0; j < NY; j++) {

    #ifdef WEST_PRESSURE_NEBB
        u_i = 1.0 - 1.0/rho_left*(f1[INDEX_3D(0,j,0)] + f1[INDEX_3D(0,j,3)] + f1[INDEX_3D(0,j,4)] + 2.0*(f1[INDEX_3D(0,j,2)] + f1[INDEX_3D(0,j,6)] + f1[INDEX_3D(0,j,8)]));

        f1[INDEX_3D(0,j,1)] = f1[INDEX_3D(0,j,2)] + 2.0/3.0*rho_left*u_i;

        f1[INDEX_3D(0,j,5)] = f1[INDEX_3D(0,j,6)] - 0.5*(f1[INDEX_3D(0,j,3)] - f1[INDEX_3D(0,j,4)]) + 1.0/6.0*rho_left*u_i;

        f1[INDEX_3D(0,j,7)] = f1[INDEX_3D(0,j,8)] + 0.5*(f1[INDEX_3D(0,j,3)] - f1[INDEX_3D(0,j,4)]) + 1.0/6.0*rho_left*u_i;
    #endif
        
    #ifdef EAST_PRESSURE_NEBB

        u_i = -1.0 + 1.0/rho_right*(f1[INDEX_3D(NX-1,j,0)] + f1[INDEX_3D(NX-1,j,3)] + f1[INDEX_3D(NX-1,j,4)] + 2.0*(f1[INDEX_3D(NX-1,j,1)] + f1[INDEX_3D(NX-1,j,5)] + f1[INDEX_3D(NX-1,j,7)]));

        f1[INDEX_3D(NX-1,j,2)] = f1[INDEX_3D(NX-1,j,1)] - 2.0/3.0*rho_right*u_i;

        f1[INDEX_3D(NX-1,j,6)] = f1[INDEX_3D(NX-1,j,5)] + 0.5*(f1[INDEX_3D(NX-1,j,3)] - f1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*rho_right*u_i;

        f1[INDEX_3D(NX-1,j,8)] = f1[INDEX_3D(NX-1,j,7)] - 0.5*(f1[INDEX_3D(NX-1,j,3)] - f1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*rho_right*u_i;
    #endif

    }

#endif

}

#endif

#ifdef TEMPERATURE

void stream_thermal_populations() {
    
    int x_i, y_i, p_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int p = 0; p < NP; p++) {

                x_i = i-cx_i[p];
                y_i = j-cy_i[p];

#ifdef XPERIODIC
                if (x_i < 0) x_i = NX-1;
                else if (x_i == NX) x_i = 0;
#endif

#ifdef WEST_TEMPERATURE_NEBB_CONSTANT_VALUE
                if (x_i < 0) continue;
#endif

#ifdef EAST_TEMPERATURE_NEBB_CONSTANT_GRADIENT
                if (x_i == NX) continue;
#endif

#ifdef SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i < 0) {
                    p_i = p_bounceback[p];
                    g1[INDEX_3D(i,j,p)] = -g2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*T_bottom;
                    continue;
                }
#endif

#ifdef NORTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i == NY) {
                    p_i = p_bounceback[p];
                    g1[INDEX_3D(i,j,p)] = -g2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*T_top;
                    continue;
                }
#endif
                
                g1[INDEX_3D(i,j,p)] = g2[INDEX_3D(x_i, y_i, p)];

            }

        }
    }

#if defined WEST_TEMPERATURE_NEBB_CONSTANT_VALUE || defined EAST_TEMPERATURE_NEBB_CONSTANT_GRADIENT

    double T_i, u_i;

#ifdef PHASECHANGE

    double phi_i;

#endif

    for (int j = 0; j < NY; j++) {

    #ifdef WEST_TEMPERATURE_NEBB_CONSTANT_VALUE
        #ifdef PHASECHANGE
            phi_i = phi[INDEX_2D(0,j)];
        
            if (phi_i == 0.0) {
                T_i = T[INDEX_2D(1,j)];
            }
            else if (phi_i == 1.0) {
                T_i = T_left;
            }
            else {
                T_i = phi_i*T_left + (1.0-phi_i)*T[INDEX_2D(1,j)];
            }
        #else
            T_i = T_left;
        #endif

        u_i = 1.0 - 1.0/T_i*(g1[INDEX_3D(0,j,0)] + g1[INDEX_3D(0,j,3)] + g1[INDEX_3D(0,j,4)] + 2.0*(g1[INDEX_3D(0,j,2)] + g1[INDEX_3D(0,j,6)] + g1[INDEX_3D(0,j,8)]));

        g1[INDEX_3D(0,j,1)] = g1[INDEX_3D(0,j,2)] + 2.0/3.0*T_i*u_i;

        g1[INDEX_3D(0,j,5)] = g1[INDEX_3D(0,j,6)] - 0.5*(g1[INDEX_3D(0,j,3)] - g1[INDEX_3D(0,j,4)]) + 1.0/6.0*T_i*u_i;

        g1[INDEX_3D(0,j,7)] = g1[INDEX_3D(0,j,8)] + 0.5*(g1[INDEX_3D(0,j,3)] - g1[INDEX_3D(0,j,4)]) + 1.0/6.0*T_i*u_i;
    #endif

    #ifdef EAST_TEMPERATURE_NEBB_CONSTANT_GRADIENT
        T_i = T[INDEX_2D(NX-2, j)];

        u_i = -1.0 + 1.0/T_i*(g1[INDEX_3D(NX-1,j,0)] + g1[INDEX_3D(NX-1,j,3)] + g1[INDEX_3D(NX-1,j,4)] + 2.0*(g1[INDEX_3D(NX-1,j,1)] + g1[INDEX_3D(NX-1,j,5)] + g1[INDEX_3D(NX-1,j,7)]));

        g1[INDEX_3D(NX-1,j,2)] = g1[INDEX_3D(NX-1,j,1)] - 2.0/3.0*T_i*u_i;

        g1[INDEX_3D(NX-1,j,6)] = g1[INDEX_3D(NX-1,j,5)] + 0.5*(g1[INDEX_3D(NX-1,j,3)] - g1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*T_i*u_i;

        g1[INDEX_3D(NX-1,j,8)] = g1[INDEX_3D(NX-1,j,7)] - 0.5*(g1[INDEX_3D(NX-1,j,3)] - g1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*T_i*u_i;
    #endif

    }

#endif

}

#endif