#include <stdio.h>  // printf, fflush, stdout
#include <stdlib.h> // malloc
#include <hdf5.h>   // Output data_*.h5 files
#include <time.h>   // time


// --- SETTINGS ---
#define NTIME   1000000          // Number of timesteps
#define NSTORE  10000            // Store macroscopic quantities after NSTORE timesteps
#define NLOG    100             // Print progress percentage after NLOG timesteps

#define NX 2                  // Number of cells in the x-direction
#define NY 2048                  // Number of cells in the y-direction
#define NP 9                    // Number of velocity directions, DON'T CHANGE!

// --- PHYSICAL CONSTANTS ---
static double tau_g_liquid = 0.50498;
static double tau_g_solid = 0.50498;

static double T_m = 0.0;
static double T_top = 0.0;
static double T_bottom = -1.0/0.95;

static double c_liquid = 0.95;
static double c_solid = 0.95;
static double L_f = 1.0;

// --- SIMULATION ---
#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)

static int cx_i[NP] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
static int cy_i[NP] = {0, 0, -1, -1, -1, 0, 1, 1, 1};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};

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


// Used for periodic boundary conditions
int mod(int x, int n) {
    if (x < 0) return n-1;
    if (x == n) return 0;
    return x;
}


int main(void) {
    // Variables used in functions
    double T_i;
    double eq;
    int x_i, y_i, p_i;
    double tau_g, phi_i, phi_old_i, Delta_phi, H_i;

    int noutput = 0;

    // Allocate memory for macroscopic quantities
    T = (double *) malloc(NX*NY*sizeof(double));
    phi = (double *) malloc(NX*NY*sizeof(double));
    phi_old = (double *) malloc(NX*NY*sizeof(double));

    // Allocate memory for populations
    g1 = (double *) malloc(NX*NY*NP*sizeof(double));
    g2 = (double *) malloc(NX*NY*NP*sizeof(double));

    // Step (i) and (ii)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            phi[INDEX_2D(i,j)] = 1.0;
            phi_old[INDEX_2D(i,j)] = 1.0;
            T[INDEX_2D(i,j)] = T_m;
            for (int p = 0; p < NP; p++) {
                g1[INDEX_3D(i,j,p)] = w[p]*T[INDEX_2D(i,j)];
            }
        }
    }

    output_data(0, noutput);
    noutput++;

    // MAIN LOOP
    for (int t = 0; t < NTIME; t++) {
        // Output macroscopic quantities
        if (((t % NSTORE) == 0) && (t > 0)) {
            output_data(t, noutput);
            noutput++;
        }

        // Collide thermal populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                T_i = T[INDEX_2D(i,j)];
                phi_i = phi[INDEX_2D(i,j)];
                phi_old_i = phi_old[INDEX_2D(i,j)];
                Delta_phi = phi_i - phi_old_i;
                tau_g = phi_i*tau_g_liquid + (1.0-phi_i)*tau_g_solid;
                for (int p = 0; p < NP; p++) {
                    eq = w[p]*T_i;
                    g2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_g)*g1[INDEX_3D(i,j,p)] + 1.0/tau_g*eq - w[p]*L_f/c_solid*Delta_phi;
                }
            }
        }

        // Stream thermal populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int p = 0; p < NP; p++) {
                    y_i = j-cy_i[p];
                    if (y_i < 0) {
                        p_i = ((p+3)%8)+1;
                        g1[INDEX_3D(i,j,p)] = -g2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*T_bottom;
                    }
                    else if (y_i == NY) {
                        p_i = ((p+3)%8)+1;
                        g1[INDEX_3D(i,j,p)] = -g2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*T_top;
                    }
                    else {
                        x_i = mod(i-cx_i[p], NX);
                        g1[INDEX_3D(i,j,p)] = g2[INDEX_3D(x_i, y_i, p)];
                    }
                }
            }
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

    output_data(NTIME, noutput);

    printf("\nDone!\n");

    return 0;
}
