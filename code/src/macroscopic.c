#include "../constants.h"
#include "../include/globals.h"
#include "../include/macroscopic.h"


void calculate_hydrodynamic_macroscopic_quantities() {

    double rho_i, u_i, v_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            rho_i = 0.0;
            u_i = 0.0;
            v_i = 0.0;
            for (int p = 0; p < NP; p++) {
                rho_i += f1[INDEX_3D(i,j,p)];
                u_i += f1[INDEX_3D(i,j,p)]*cx[p];
                v_i += f1[INDEX_3D(i,j,p)]*cy[p];
            }
            rho[INDEX_2D(i,j)] = rho_i;
            u[INDEX_2D(i,j)] = (u_i + 0.5*Fx[INDEX_2D(i,j)]) / rho_i;
            v[INDEX_2D(i,j)] = (v_i + 0.5*Fy[INDEX_2D(i,j)]) / rho_i;

        }
    }
}