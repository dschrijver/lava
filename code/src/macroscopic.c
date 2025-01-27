#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
#include "../include/macroscopic.h"


#ifdef FLOW

void calculate_hydrodynamic_macroscopic_quantities() {

    double rho_i, u_i, v_i;

#ifdef DUALCOMPONENT

    double rho_lava_i, rho_air_i;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            rho_i = 0.0;

#ifdef DUALCOMPONENT

            rho_lava_i = 0.0;
            rho_air_i = 0.0;       

#endif

            u_i = 0.0;
            v_i = 0.0;

            for (int p = 0; p < NP; p++) {

                

#ifndef DUALCOMPONENT

                rho_i += f1[INDEX_3D(i,j,p)];
                u_i += f1[INDEX_3D(i,j,p)]*cx[p];
                v_i += f1[INDEX_3D(i,j,p)]*cy[p];

#endif

#ifdef DUALCOMPONENT

                rho_i += f1_lava[INDEX_3D(i,j,p)] + f1_air[INDEX_3D(i,j,p)];
                rho_lava_i += f1_lava[INDEX_3D(i,j,p)];
                rho_air_i += f1_air[INDEX_3D(i,j,p)];
                u_i += (f1_lava[INDEX_3D(i,j,p)] + f1_air[INDEX_3D(i,j,p)])*cx[p];
                v_i += (f1_lava[INDEX_3D(i,j,p)] + f1_air[INDEX_3D(i,j,p)])*cy[p];

#endif

            }

            rho[INDEX_2D(i,j)] = rho_i;

            u[INDEX_2D(i,j)] = u_i / rho_i;
            v[INDEX_2D(i,j)] = v_i / rho_i;

        }
    }
}

#endif

#ifdef FLOW

void shift_velocity() {

    double rho_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            rho_i = rho[INDEX_2D(i,j)];

#ifndef DUALCOMPONENT

            u[INDEX_2D(i,j)] += Fx[INDEX_2D(i,j)]/(2.0*rho_i);
            v[INDEX_2D(i,j)] += Fy[INDEX_2D(i,j)]/(2.0*rho_i);

#endif

#ifdef DUALCOMPONENT

            u[INDEX_2D(i,j)] += (Fx_lava[INDEX_2D(i,j)] + Fx_air[INDEX_2D(i,j)])/(2.0*rho_i);
            v[INDEX_2D(i,j)] += (Fy_lava[INDEX_2D(i,j)] + Fy_air[INDEX_2D(i,j)])/(2.0*rho_i);

#endif

        }
    }

}

#endif

#ifdef TEMPERATURE

void calculate_thermal_macroscopic_quantities() {

    double T_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            T_i = 0.0;
            
            for (int p = 0; p < NP; p++) {
                T_i += g1[INDEX_3D(i,j,p)];
            
            }
            
            T[INDEX_2D(i,j)] = T_i;
        
        }
    }
    
}

#endif

#ifdef PHASECHANGE

void calculate_enthalpy_and_liquid_fraction() {

    double phi_i, T_i, H_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            phi_i = phi[INDEX_2D(i,j)];
            phi_old[INDEX_2D(i,j)] = phi_i;
            T_i = T[INDEX_2D(i,j)];
            H_i = (1.0-phi_i)*c_solid*T_i + phi_i*(c_liquid*(T_i-T_melt) + c_solid*T_melt) + L_f*phi_i;

            if (H_i < c_solid*T_melt) {
                phi[INDEX_2D(i,j)] = 0.0;
            }
            else if ((c_solid*T_melt <= H_i) && (H_i <= c_solid*T_melt + L_f)) {
                phi[INDEX_2D(i,j)] = (H_i-c_solid*T_melt)/L_f;
            }
            else {
                phi[INDEX_2D(i,j)] = 1.0;
            }

        }
    }

}

#endif