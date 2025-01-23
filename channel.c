#include <stdio.h>  // printf, fflush, stdout
#include <stdlib.h> // malloc
#include <hdf5.h>   // Output data_*.h5 files
#include <time.h>   // time
#include <raylib.h>
#include <math.h>   // fmin, fmax


// --- SETTINGS ---
#define NTIME   10          // Number of timesteps
#define NSTORE  1            // Store macroscopic quantities after NSTORE timesteps
#define NLOG    1             // Print progress percentage after NLOG timesteps

#define NX 32                  // Number of cells in the x-direction
#define NY 128                  // Number of cells in the y-direction
#define NP 9                    // Number of velocity directions, DON'T CHANGE!

static double tau = 0.899;
static double tau_g_liquid = 0.53;
static double tau_g_solid = 0.53;

static double T_m = 0.0;
static double T_top = -1.0;
static double T_bottom = -1.0;
static double T_inlet = 0.0526;
static double T_initial = 0.0526;

static double c_liquid = 0.95;
static double c_solid = 0.95;
static double L_f = 1.0;

static double Delta_p = 1e-5;


// --- DISPLAY ---
#undef ANIM
#define  STORE
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
const int cell_size = 2;
const int steps_per_frame = 100;


// --- SIMULATION ---
#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)

static double cs2 = 1.0/3.0;
static int cx_i[NP] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
static int cy_i[NP] = {0, 0, 0, 1, -1, 1, -1, -1, 1};
static double cx[NP] = {0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
static double cy[NP] = {0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
static int p_bounceback[NP] = {0, 2, 1, 4, 3, 6, 5, 8, 7};

double *rho, *u, *v, *Fx, *Fy, *f1, *f2;
double *T, *g1, *g2;
double *phi, *phi_old;


// Output data_*.h5 files
void output_data(int t, int n_output) {
    char filename[32];
    sprintf(filename, "data_%d.h5", n_output);
    hid_t hdf5_fp = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dim_1d[1] = {1};
    hid_t dataspace_id_1d = H5Screate_simple(1, dim_1d, NULL);
    hsize_t dims[2] = {NX, NY};
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

    hid_t dataset_id;
    dataset_id = H5Dcreate2(hdf5_fp, "time", H5T_NATIVE_INT32, dataspace_id_1d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
    dataset_id = H5Dcreate2(hdf5_fp, "rho", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);
    dataset_id = H5Dcreate2(hdf5_fp, "u", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
    dataset_id = H5Dcreate2(hdf5_fp, "v", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
    dataset_id = H5Dcreate2(hdf5_fp, "Fx", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fx);
    dataset_id = H5Dcreate2(hdf5_fp, "Fy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fy);
    dataset_id = H5Dcreate2(hdf5_fp, "T", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, T);
    dataset_id = H5Dcreate2(hdf5_fp, "phi", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);

    H5Fflush(hdf5_fp, H5F_SCOPE_GLOBAL);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(dataspace_id_1d);
    H5Fclose(hdf5_fp);
}


int render_frame() {
    char phi_i;
    Color cell_color = {0,0,0,255};

    PollInputEvents(); // Poll input events (SUPPORT_CUSTOM_FRAME_CONTROL)
    if (WindowShouldClose()) {
        return 0;
    }

    BeginDrawing();

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                phi_i = (char)(fmin(fmax(phi[INDEX_2D(i,j)], 0.0), 1.0) * 255.0);
                cell_color.r = phi_i;
                cell_color.b = 255 - phi_i;
                DrawRectangle(i*cell_size, (NY-1-j)*cell_size, cell_size, cell_size, cell_color);
            }
        }

    EndDrawing();
    SwapScreenBuffer();

    return 1;
}


int main(void) {
#ifdef STORE
    int noutput = 0;
#endif

#ifdef ANIM
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "");
#endif

    // Variables used in functions
    double rho_i, u_i, v_i;
    double Fx_i, Fy_i;
    double T_i;
    double u2, uc, eq, S;
    int x_i, y_i, p_i;
    double tau_g, phi_i, phi_old_i, Delta_phi, phi2, H_i;

    // Allocate memory for macroscopic quantities
    rho = (double *) malloc(NX*NY*sizeof(double));
    u = (double *) malloc(NX*NY*sizeof(double));
    v = (double *) malloc(NX*NY*sizeof(double));
    Fx = (double *) malloc(NX*NY*sizeof(double));
    Fy = (double *) malloc(NX*NY*sizeof(double));
    T = (double *) malloc(NX*NY*sizeof(double));
    phi = (double *) malloc(NX*NY*sizeof(double));
    phi_old = (double *) malloc(NX*NY*sizeof(double));

    // Allocate memory for populations
    f1 = (double *) malloc(NX*NY*NP*sizeof(double));
    f2 = (double *) malloc(NX*NY*NP*sizeof(double));
    g1 = (double *) malloc(NX*NY*NP*sizeof(double));
    g2 = (double *) malloc(NX*NY*NP*sizeof(double));

    // Choose initial density and velocity fields and compute initial populations
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            rho[INDEX_2D(i,j)] = 1.0 - (Delta_p/cs2)/(double)(NX-1) * (double)i;
            u[INDEX_2D(i,j)] = 0.0;
            v[INDEX_2D(i,j)] = 0.0;
            T[INDEX_2D(i,j)] = T_initial;
            phi[INDEX_2D(i,j)] = 1.0;
            phi_old[INDEX_2D(i,j)] = 1.0;
            for (int p = 0; p < NP; p++) {
                f1[INDEX_3D(i,j,p)] = w[p]*rho[INDEX_2D(i,j)];
                g1[INDEX_3D(i,j,p)] = w[p]*T_initial;
            }
        }
    }

    // Initialize the force fields
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            Fx[INDEX_2D(i,j)] = 0.0;
            Fy[INDEX_2D(i,j)] = 0.0;
        }
    }

#ifdef STORE
    // Output initial macroscopic quantities
    output_data(0, noutput);
    noutput++;
#endif

    // MAIN LOOP
    for (int t = 0; t < NTIME; t++) {
#ifdef ANIM
        if ((t % steps_per_frame) == 0) {
            if (!render_frame()) break;
        }
#endif

#ifdef STORE
        // Output macroscopic quantities
        if (((t % NSTORE) == 0) && (t > 0)) {
            output_data(t, noutput);
            noutput++;
        }
#endif

        // Collide thermal populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                u_i = u[INDEX_2D(i,j)];
                v_i = v[INDEX_2D(i,j)];
                T_i = T[INDEX_2D(i,j)];
                u2 = u_i*u_i + v_i*v_i;
                phi_i = phi[INDEX_2D(i,j)];
                phi_old_i = phi_old[INDEX_2D(i,j)];
                Delta_phi = phi_i - phi_old_i;
                tau_g = phi_i*tau_g_liquid + (1.0-phi_i)*tau_g_solid;
                for (int p = 0; p < NP; p++) {
                    uc = u_i*cx[p] + v_i*cy[p];
                    eq = w[p]*T_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    g2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_g)*g1[INDEX_3D(i,j,p)] + 1.0/tau_g*eq - w[p]*L_f/c_solid*Delta_phi;
                }
            }
        }

        // Stream thermal populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int p = 0; p < NP; p++) {
                    x_i = i-cx_i[p];
                    y_i = j-cy_i[p];
                    p_i = p_bounceback[p];
                    if ((x_i < 0) || (x_i == NX)) {
                        continue;
                    }
                    else if (y_i < 0) {
                        g1[INDEX_3D(i,j,p)] = -g2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*T_bottom;
                    }
                    else if (y_i == NY) {
                        g1[INDEX_3D(i,j,p)] = -g2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*T_top;
                    }
                    else {
                        g1[INDEX_3D(i,j,p)] = g2[INDEX_3D(x_i, y_i, p)];
                    }
                }
            }
        }

        for (int j = 0; j < NY; j++) {
            if (phi[INDEX_2D(0,j)] == 0.0) {
                T_i = T[INDEX_2D(0,j)];
            }
            else if (phi[INDEX_2D(0,j)] == 1.0) {
                T_i = T_inlet;
            }
            else {
                T_i = phi[INDEX_2D(0,j)]*T_inlet + (1.0-phi[INDEX_2D(0,j)])*T[INDEX_2D(0,j)];
            }

            u_i = 1.0 - 1.0/T_i*(g1[INDEX_3D(0,j,0)] + g1[INDEX_3D(0,j,3)] + g1[INDEX_3D(0,j,4)] + 2.0*(g1[INDEX_3D(0,j,2)] + g1[INDEX_3D(0,j,6)] + g1[INDEX_3D(0,j,8)]));

            g1[INDEX_3D(0,j,1)] = g1[INDEX_3D(0,j,2)] + 2.0/3.0*T_i*u_i;

            g1[INDEX_3D(0,j,5)] = g1[INDEX_3D(0,j,6)] - 0.5*(g1[INDEX_3D(0,j,3)] - g1[INDEX_3D(0,j,4)]) + 1.0/6.0*T_i*u_i;

            g1[INDEX_3D(0,j,7)] = g1[INDEX_3D(0,j,8)] + 0.5*(g1[INDEX_3D(0,j,3)] - g1[INDEX_3D(0,j,4)]) + 1.0/6.0*T_i*u_i;

            T_i = T[INDEX_2D(NX-1, j)];

            u_i = -1.0 + 1.0/T_i*(g1[INDEX_3D(NX-1,j,0)] + g1[INDEX_3D(NX-1,j,3)] + g1[INDEX_3D(NX-1,j,4)] + 2.0*(g1[INDEX_3D(NX-1,j,1)] + g1[INDEX_3D(NX-1,j,5)] + g1[INDEX_3D(NX-1,j,7)]));

            g1[INDEX_3D(NX-1,j,2)] = g1[INDEX_3D(NX-1,j,1)] - 2.0/3.0*T_i*u_i;

            g1[INDEX_3D(NX-1,j,6)] = g1[INDEX_3D(NX-1,j,5)] + 0.5*(g1[INDEX_3D(NX-1,j,3)] - g1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*T_i*u_i;

            g1[INDEX_3D(NX-1,j,8)] = g1[INDEX_3D(NX-1,j,7)] - 0.5*(g1[INDEX_3D(NX-1,j,3)] - g1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*T_i*u_i;
        }

        // Calculate temperature
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                T_i = 0.0;
                for (int p = 0; p < NP; p++) {
                    T_i += g1[INDEX_3D(i,j,p)];
                }
                T[INDEX_2D(i,j)] = T_i;
            }
        }

        // Determine penalization force
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                phi_i = phi[INDEX_2D(i,j)];
                phi2 = phi_i*phi_i;
                rho_i = rho[INDEX_2D(i,j)];
                
                Fx[INDEX_2D(i,j)] = -(1.0-phi2)*rho_i*u[INDEX_2D(i,j)];
                Fy[INDEX_2D(i,j)] = -(1.0-phi2)*rho_i*v[INDEX_2D(i,j)];
            }
        }

        // Collide populations
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
                    eq = w[p]*rho_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    S = (1.0 - 1.0/(2.0*tau))*w[p]*(((cx[p]-u_i)/cs2 + uc/(cs2*cs2)*cx[p])*Fx_i + 
                                                    ((cy[p]-v_i)/cs2 + uc/(cs2*cs2)*cy[p])*Fy_i);
                    f2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau)*f1[INDEX_3D(i,j,p)] + 1.0/tau*eq + S;
                }
            }
        }

        // Stream populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int p = 0; p < NP; p++) {
                    x_i = i-cx_i[p];
                    y_i = j-cy_i[p];
                    p_i = p_bounceback[p];
                    if ((x_i < 0) || (x_i == NX)) {
                        continue;
                    }
                    else if ((y_i < 0) || (y_i == NY)) {
                        f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(i, j, p_i)];
                    }
                    else {
                        f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(x_i, y_i, p)];
                    }
                }
            }
        }

        for (int j = 0; j < NY; j++) {
            u_i = 1.0 - (f1[INDEX_3D(0,j,0)] + f1[INDEX_3D(0,j,3)] + f1[INDEX_3D(0,j,4)] + 2.0*(f1[INDEX_3D(0,j,2)] + f1[INDEX_3D(0,j,6)] + f1[INDEX_3D(0,j,8)]));

            f1[INDEX_3D(0,j,1)] = f1[INDEX_3D(0,j,2)] + 2.0/3.0*u_i;

            f1[INDEX_3D(0,j,5)] = f1[INDEX_3D(0,j,6)] - 0.5*(f1[INDEX_3D(0,j,3)] - f1[INDEX_3D(0,j,4)]) + 1.0/6.0*u_i;

            f1[INDEX_3D(0,j,7)] = f1[INDEX_3D(0,j,8)] + 0.5*(f1[INDEX_3D(0,j,3)] - f1[INDEX_3D(0,j,4)]) + 1.0/6.0*u_i;

            u_i = -1.0 + 1.0/(1.0-Delta_p/cs2)*(f1[INDEX_3D(NX-1,j,0)] + f1[INDEX_3D(NX-1,j,3)] + f1[INDEX_3D(NX-1,j,4)] + 2.0*(f1[INDEX_3D(NX-1,j,1)] + f1[INDEX_3D(NX-1,j,5)] + f1[INDEX_3D(NX-1,j,7)]));

            f1[INDEX_3D(NX-1,j,2)] = f1[INDEX_3D(NX-1,j,1)] - 2.0/3.0*(1.0-Delta_p/cs2)*u_i;

            f1[INDEX_3D(NX-1,j,6)] = f1[INDEX_3D(NX-1,j,5)] + 0.5*(f1[INDEX_3D(NX-1,j,3)] - f1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*(1.0-Delta_p/cs2)*u_i;

            f1[INDEX_3D(NX-1,j,8)] = f1[INDEX_3D(NX-1,j,7)] - 0.5*(f1[INDEX_3D(NX-1,j,3)] - f1[INDEX_3D(NX-1,j,4)]) - 1.0/6.0*(1.0-Delta_p/cs2)*u_i;
        }

        // Compute macroscopic quantities
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

        // Compute the enthalpy and liquid fraction
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                phi_i = phi[INDEX_2D(i,j)];
                phi_old[INDEX_2D(i,j)] = phi_i;
                T_i = T[INDEX_2D(i,j)];
                H_i = (1.0-phi_i)*c_solid*T_i + phi_i*(c_liquid*(T_i-T_m) + c_solid*T_m) + L_f*phi_i;

                if (H_i < c_solid*T_m) {
                    phi[INDEX_2D(i,j)] = 0.0;
                }
                else if ((c_solid*T_m <= H_i) && (H_i <= c_solid*T_m + L_f)) {
                    phi[INDEX_2D(i,j)] = (H_i-c_solid*T_m)/L_f;
                }
                else {
                    phi[INDEX_2D(i,j)] = 1.0;
                }
            }
        }

        // Log the progress
        if (((t+1) % NLOG) == 0) {
            printf("\rProgress: %.1f%%", ((double)(t+1))/((double)NTIME)*100.0);
            fflush(stdout);
        }
    }

#ifdef STORE
    // Output the final macroscopic quantities
    output_data(NTIME, noutput);
#endif

#ifdef ANIM
    CloseWindow();
#endif

    printf("\nDone!\n");

    return 0;
}
