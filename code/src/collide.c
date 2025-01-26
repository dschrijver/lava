#include "../constants.h"
#include "../include/globals.h"
#include "../include/stencil.h"
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


void collide_thermal_populations() {

    double T_i, geq;

#ifdef FLOW

    double u_i, v_i, u2, uc;

#endif

#ifdef PHASECHANGE

    double phi_i, Delta_phi, tau_g_effective;

#endif

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            T_i = T[INDEX_2D(i,j)];

#ifdef FLOW

            u_i = u[INDEX_2D(i,j)];
            v_i = v[INDEX_2D(i,j)];
            u2 = u_i*u_i + v_i*v_i;

#endif

#ifdef PHASECHANGE

            phi_i = phi[INDEX_2D(i,j)];
            Delta_phi = phi_i - phi_old[INDEX_2D(i,j)];
            tau_g_effective = phi_i*tau_g_liquid + (1.0-phi_i)*tau_g_solid;

#endif
            
            for (int p = 0; p < NP; p++) {

#ifdef FLOW
                uc = u_i*cx[p] + v_i*cy[p];

#endif

#ifdef FLOW

                geq = w[p]*T_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));

#else

                geq = w[p]*T_i;

#endif

#ifdef PHASECHANGE

                g2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_g_effective)*g1[INDEX_3D(i,j,p)] + 1.0/tau_g_effective*geq - w[p]*L_f/c_solid*Delta_phi;

#else

                g2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_g)*g1[INDEX_3D(i,j,p)] + 1.0/tau_g*geq;

#endif

            }
        }

    }

}