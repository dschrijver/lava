#include <math.h>   // sin()

#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
#include "../include/initialize.h"


void initialize() {

#ifdef FLOW

    double u2, uc;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            
#ifdef FLOW
    
    #ifdef DENSITY_UNITY
            rho[INDEX_2D(i,j)] = 1.0;
    #endif

    #ifdef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
            rho[INDEX_2D(i,j)] = rho_ini_left + (rho_ini_right - rho_ini_left)/(double)(NX-1) * (double)i;
    #endif

    #ifdef VELOCITY_ZERO
            u[INDEX_2D(i,j)] = 0.0;
            v[INDEX_2D(i,j)] = 0.0;
    #endif

    #ifdef VELOCITY_SINE
            u[INDEX_2D(i,j)] = u_ini_amplitude*sin(2.0*M_PI*u_ini_frequency*((double)i + 0.5)/(double)NX);
            v[INDEX_2D(i,j)] = v_ini_amplitude*sin(2.0*M_PI*v_ini_frequency*((double)j + 0.5)/(double)NY);
    #endif

#endif

#ifdef TEMPERATURE

    #ifdef TEMPERATURE_CONSTANT_VALUE
            T[INDEX_2D(i,j)] = T_ini;
    #endif

    #ifdef TEMPERATURE_CONSTANT_VERTICAL_GRADIENT
            T[INDEX_2D(i,j)] = (T_ini_top-T_ini_bottom)/(double)NY*((double)j + 0.5) + T_bottom;
    #endif

    #ifdef TEMPERATURE_CHANNEL
            if (((double)j + 0.5 < 0.5*((double)NY-channel_width)) || 
                ((double)j + 0.5 > 0.5*((double)NY-channel_width) + channel_width)) {
                T[INDEX_2D(i,j)] = T_ini_solid;
            }
            else {
                T[INDEX_2D(i,j)] = T_ini_liquid; 
            }
    #endif

#endif

#ifdef PHASECHANGE

    #ifdef PHI_CONSTANT_VALUE
            phi[INDEX_2D(i,j)] = phi_ini;
            phi_old[INDEX_2D(i,j)] = phi_ini;
    #endif
    
    #ifdef PHI_CHANNEL
            if (((double)j + 0.5 < 0.5*((double)NY-channel_width)) || 
                ((double)j + 0.5 > 0.5*((double)NY-channel_width) + channel_width)) {
                phi[INDEX_2D(i,j)] = 0.0;
                phi_old[INDEX_2D(i,j)] = 0.0;
            }
            else {
                phi[INDEX_2D(i,j)] = 1.0;
                phi_old[INDEX_2D(i,j)] = 1.0;
            }
    #endif

#endif

#ifdef FLOW
            u2 = u[INDEX_2D(i,j)]*u[INDEX_2D(i,j)] + v[INDEX_2D(i,j)]*v[INDEX_2D(i,j)];
#endif

            for (int p = 0; p < NP; p++) {

#ifdef FLOW
                uc = u[INDEX_2D(i,j)]*cx[p] + v[INDEX_2D(i,j)]*cy[p];
#endif

#ifdef FLOW

                f1[INDEX_3D(i,j,p)] = w[p]*rho[INDEX_2D(i,j)]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));;

#endif

#ifdef TEMPERATURE

    #ifdef FLOW
                g1[INDEX_3D(i,j,p)] = w[p]*T[INDEX_2D(i,j)]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
    #else
                g1[INDEX_3D(i,j,p)] = w[p]*T[INDEX_2D(i,j)];
    #endif

#endif
                
            }
            
        }
    }
    
}

#ifdef FLOW

void shift_velocity() {

    double rho_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            rho_i = rho[INDEX_2D(i,j)];
            u[INDEX_2D(i,j)] -= Fx[INDEX_2D(i,j)]/(2.0*rho_i);
            v[INDEX_2D(i,j)] -= Fy[INDEX_2D(i,j)]/(2.0*rho_i);

        }
    }

}

#endif