#include "../constants.h"
#include "../include/globals.h"
#include "../include/initialize.h"


void initialize() {

    double u2, uc;

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

            Fx[INDEX_2D(i,j)] = Fx_body;
            Fy[INDEX_2D(i,j)] = Fy_body;

#endif

            u2 = u[INDEX_2D(i,j)]*u[INDEX_2D(i,j)] + v[INDEX_2D(i,j)]*v[INDEX_2D(i,j)];

            for (int p = 0; p < NP; p++) {
                
                uc = u[INDEX_2D(i,j)]*cx[p] + v[INDEX_2D(i,j)]*cy[p];

#ifdef FLOW

                f1[INDEX_3D(i,j,p)] = w[p]*rho[INDEX_2D(i,j)]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));;

#endif
                
            }
            
        }
    }
    
}


void shift_velocity() {

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            u[INDEX_2D(i,j)] = -Fx[INDEX_2D(i,j)]/(2.0*rho[INDEX_2D(i,j)]);
            v[INDEX_2D(i,j)] = -Fy[INDEX_2D(i,j)]/(2.0*rho[INDEX_2D(i,j)]);

        }
    }

}