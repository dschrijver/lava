#include <math.h>   // sin()

#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
#include "../include/initialize.h"


void initialize_macroscopic_quantities() {

#ifdef DENSITY_TWO_COMPONENT_POISEUILLE

    double y;

#endif

#ifdef VELOCITY_CHANNEL

    double y, y_shift;

    double d = channel_width / 2.0;
    double nu = cs2*(tau - 0.5);
    double offset = ((double)NY-channel_width) / 2.0;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
    
#ifdef DENSITY_UNITY
    
            rho[INDEX_2D(i,j)] = 1.0;

#endif

#ifdef DENSITY_CONSTANT_HORIZONTAL_GRADIENT

            rho[INDEX_2D(i,j)] = rho_ini_left + (rho_ini_right - rho_ini_left)/(double)(NX-1) * (double)i;

#endif

#ifdef DENSITY_TWO_COMPONENT_POISEUILLE

            rho[INDEX_2D(i,j)] = 1.0;
            y = (double)j + 0.5;
            
            if((y > ((double)NY-center_width) / 2.0) && 
               (y < ((double)NY-center_width) / 2.0 + center_width)) {
                rho_lava[INDEX_2D(i,j)] = 0.0;
                rho_air[INDEX_2D(i,j)] = 1.0;
            }
            else {
                rho_lava[INDEX_2D(i,j)] = 1.0;
                rho_air[INDEX_2D(i,j)] = 0.0;
            }

#endif

#ifdef VELOCITY_ZERO

            u[INDEX_2D(i,j)] = 0.0;
            v[INDEX_2D(i,j)] = 0.0;

#endif

#ifdef VELOCITY_SINE

            u[INDEX_2D(i,j)] = u_ini_amplitude*sin(2.0*M_PI*u_ini_frequency*((double)i + 0.5)/(double)NX);
            v[INDEX_2D(i,j)] = v_ini_amplitude*sin(2.0*M_PI*v_ini_frequency*((double)j + 0.5)/(double)NY);

#endif

#ifdef VELOCITY_CHANNEL

            y = (double)j + 0.5;

            if ((y > offset) && (y < offset + channel_width)) {
                y_shift = (double)j + 0.5 - d - offset;
                u[INDEX_2D(i,j)] = dp_dx/(2.0*nu)*(d*d - y_shift*y_shift);
            } else {
                u[INDEX_2D(i,j)] = 0.0;
            }
            
            v[INDEX_2D(i,j)] = 0.0;

#endif


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
            
        }
    }
    
}


void initialize_distribution_functions() {

#ifdef FLOW

    double rho_i, uhat, vhat, u2, uc;

    #ifdef DUALCOMPONENT
        double eq, rho_lava_i, rho_air_i;
    #else
        double feq;
    #endif

#endif

#ifdef TEMPERATURE

    double T_i;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

#ifdef FLOW

            rho_i  = rho[INDEX_2D(i,j)];

    #ifdef DUALCOMPONENT
            rho_lava_i = rho_lava[INDEX_2D(i,j)];
            rho_air_i = rho_air[INDEX_2D(i,j)];

            uhat = u[INDEX_2D(i,j)] - 1.0/(2.0*rho_i)*(Fx_lava[INDEX_2D(i,j)] + Fx_air[INDEX_2D(i,j)]);
            vhat = v[INDEX_2D(i,j)] - 1.0/(2.0*rho_i)*(Fy_lava[INDEX_2D(i,j)] + Fy_air[INDEX_2D(i,j)]);
    #else
            uhat = u[INDEX_2D(i,j)] - 1.0/(2.0*rho_i)*Fx[INDEX_2D(i,j)];
            vhat = v[INDEX_2D(i,j)] - 1.0/(2.0*rho_i)*Fy[INDEX_2D(i,j)];
    #endif

#endif

#ifdef TEMPERATURE

            T_i = T[INDEX_2D(i,j)];

#endif

#ifdef FLOW

            u2 = uhat*uhat + vhat*vhat;

#endif

            for (int p = 0; p < NP; p++) {

#ifdef FLOW

                uc = uhat*cx[p] + vhat*cy[p];

    #ifdef DUALCOMPONENT
                eq = w[p]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));

                f1_lava[INDEX_3D(i,j,p)] = rho_lava_i*eq;
                f1_air[INDEX_3D(i,j,p)] = rho_air_i*eq;
    #else
                f1[INDEX_3D(i,j,p)] = w[p]*rho_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
    #endif

#endif

#ifdef TEMPERATURE

    #ifdef FLOW
                g1[INDEX_3D(i,j,p)] = w[p]*T_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
    #else
                g1[INDEX_3D(i,j,p)] = w[p]*T_i;
    #endif

#endif

            }

        }
    }

}