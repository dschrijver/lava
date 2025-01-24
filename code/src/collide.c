#include "../constants.h"
#include "../include/globals.h"
#include "../include/collide.h"


void collide_hydrodynamic_populations() {

    double rho_i, u_i, v_i, Fx_i, Fy_i, u2, uc, feq, S;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            rho_i = rho[INDEX_2D(i,j)];
            u_i = u[INDEX_2D(i,j)];
            v_i = v[INDEX_2D(i,j)];
            Fx_i = Fx[INDEX_2D(i,j)];
            Fy_i = Fy[INDEX_2D(i,j)];
            u2 = u_i*u_i + v_i*v_i;

            for (int p = 0; p < NP; p++) {

                uc = u_i*cx[p] + v_i*cy[p];
                feq = w[p]*rho_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                S = (1.0 - 1.0/(2.0*tau))*w[p]*(((cx[p]-u_i)/cs2 + uc/(cs2*cs2)*cx[p])*Fx_i + 
                                                ((cy[p]-u_i)/cs2 + uc/(cs2*cs2)*cy[p])*Fy_i);
                f2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau)*f1[INDEX_3D(i,j,p)] + 1.0/tau*feq + S;

            }

        }
    }

}