#include <stdio.h>  // printf, fflush, stdout
#include <stdlib.h> // malloc
#include <hdf5.h>   // Output data_*.h5 files


#define NTIME   100000          // Number of timesteps
#define NSTORE  1000            // Store macroscopic quantities after NSTORE timesteps
#define NLOG    100             // Print progress percentage after NLOG timesteps

#define NX 32                   // Number of cells in the x-direction
#define NY 128                  // Number of cells in the y-direction
#define NP 9                    // Number of velocity directions, DON'T CHANGE!

static double tau_lava = 1.0;   // Relaxation time of lava
static double tau_air = 0.625;  // Relaxation time of air
static double F_p = 1e-5;       // Body force from the external pressure
static double G = 4.0;          // Interaction strength between different fluid components


#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)

static double cs2 = 1.0/3.0;
static int cx_i[NP] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
static int cy_i[NP] = {0, 0, -1, -1, -1, 0, 1, 1, 1};
static double cx[NP] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
static double cy[NP] = {0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};

double *rho, *rho_lava, *rho_air, *u, *v;
double *Fx_lava, *Fy_lava, *Fx_air, *Fy_air;
double *f1_lava, *f2_lava, *f1_air, *f2_air;

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
    dataset_id = H5Dcreate2(hdf5_fp, "rho_lava", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_lava);
    dataset_id = H5Dcreate2(hdf5_fp, "rho_air", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_air);
    dataset_id = H5Dcreate2(hdf5_fp, "u", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
    dataset_id = H5Dcreate2(hdf5_fp, "v", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
    dataset_id = H5Dcreate2(hdf5_fp, "Fx_lava", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fx_lava);
    dataset_id = H5Dcreate2(hdf5_fp, "Fy_lava", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fy_lava);
    dataset_id = H5Dcreate2(hdf5_fp, "Fx_air", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fx_air);
    dataset_id = H5Dcreate2(hdf5_fp, "Fy_air", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fy_air);

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
    double rho_i, rho_air_i, rho_lava_i, u_i, v_i;
    double Fx_lava_i, Fy_lava_i, Fx_air_i, Fy_air_i; 
    double sx, sy, u2, uc, feq, S;
    int x_i, y_i, p_i;
    int noutput = 0;

    // Allocate memory for macroscopic quantities
    rho = (double *) malloc(NX*NY*sizeof(double));
    rho_lava = (double *) malloc(NX*NY*sizeof(double));
    rho_air = (double *) malloc(NX*NY*sizeof(double));
    u = (double *) malloc(NX*NY*sizeof(double));
    v = (double *) malloc(NX*NY*sizeof(double));
    Fx_lava = (double *) malloc(NX*NY*sizeof(double));
    Fy_lava = (double *) malloc(NX*NY*sizeof(double));
    Fx_air = (double *) malloc(NX*NY*sizeof(double));
    Fy_air = (double *) malloc(NX*NY*sizeof(double));

    // Allocate memory for populations
    f1_lava = (double *) malloc(NX*NY*NP*sizeof(double));
    f2_lava = (double *) malloc(NX*NY*NP*sizeof(double));
    f1_air = (double *) malloc(NX*NY*NP*sizeof(double));
    f2_air = (double *) malloc(NX*NY*NP*sizeof(double));

    // Choose initial density and velocity fields and compute initial populations
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            rho[INDEX_2D(i,j)] = 1.0;
            if ((j < NY/4) || (j >= 3*NY/4)) {
                rho_lava[INDEX_2D(i,j)] = 1.0;
                rho_air[INDEX_2D(i,j)] = 0.0;
            }
            else {
                rho_lava[INDEX_2D(i,j)] = 0.0;
                rho_air[INDEX_2D(i,j)] = 1.0;
            }
            u[INDEX_2D(i,j)] = 0.0;
            v[INDEX_2D(i,j)] = 0.0;
            for (int p = 0; p < NP; p++) {
                f1_lava[INDEX_3D(i,j,p)] = w[p]*rho_lava[INDEX_2D(i,j)];
                f1_air[INDEX_3D(i,j,p)] = w[p]*rho_air[INDEX_2D(i,j)];
            }
        }
    }

    // Initialize the force fields
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            Fx_lava[INDEX_2D(i,j)] = 0.0;
            Fy_lava[INDEX_2D(i,j)] = 0.0;
            Fx_air[INDEX_2D(i,j)] = 0.0;
            Fy_air[INDEX_2D(i,j)] = 0.0;
            for (int p = 0; p < NP; p++) {
                x_i = mod(i + cx[p], NX);
                y_i = j + cy[p];
                if ((y_i < 0) || (y_i == NY)) {
                    x_i = i;
                    y_i = j;
                }
                Fx_lava[INDEX_2D(i,j)] += w[p]*rho_air[INDEX_2D(x_i,y_i)]*cx[p];
                Fy_lava[INDEX_2D(i,j)] += w[p]*rho_air[INDEX_2D(x_i,y_i)]*cy[p];
                Fx_air[INDEX_2D(i,j)] += w[p]*rho_lava[INDEX_2D(x_i,y_i)]*cx[p];
                Fy_air[INDEX_2D(i,j)] += w[p]*rho_lava[INDEX_2D(x_i,y_i)]*cy[p];
            }
            Fx_lava[INDEX_2D(i,j)] = -G*rho_lava[INDEX_2D(i,j)]*Fx_lava[INDEX_2D(i,j)] + rho_lava[INDEX_2D(i,j)]/rho[INDEX_2D(i,j)]*F_p;
            Fy_lava[INDEX_2D(i,j)] = -G*rho_lava[INDEX_2D(i,j)]*Fy_lava[INDEX_2D(i,j)];
            Fx_air[INDEX_2D(i,j)] = -G*rho_air[INDEX_2D(i,j)]*Fx_air[INDEX_2D(i,j)] + rho_air[INDEX_2D(i,j)]/rho[INDEX_2D(i,j)]*F_p;
            Fy_air[INDEX_2D(i,j)] = -G*rho_air[INDEX_2D(i,j)]*Fy_air[INDEX_2D(i,j)];
        }
    }

    // Output initial macroscopic quantities
    output_data(0, noutput);
    noutput++;

    // Shift the velocity fields
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            u[INDEX_2D(i,j)] -= 1.0/(2.0*rho[INDEX_2D(i,j)])*(Fx_lava[INDEX_2D(i,j)] + Fx_air[INDEX_2D(i,j)]);
            v[INDEX_2D(i,j)] -= 1.0/(2.0*rho[INDEX_2D(i,j)])*(Fy_lava[INDEX_2D(i,j)] + Fy_air[INDEX_2D(i,j)]);
        }
    }

    // MAIN LOOP
    for (int t = 0; t < NTIME; t++) {

        // Output macroscopic quantities
        if (((t % NSTORE) == 0) && (t > 0)) {
            output_data(t, noutput);
            noutput++;
        }

        // Collide populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                rho_lava_i = rho_lava[INDEX_2D(i,j)];
                rho_air_i = rho_air[INDEX_2D(i,j)];
                u_i = u[INDEX_2D(i,j)];
                v_i = v[INDEX_2D(i,j)];
                Fx_lava_i = Fx_lava[INDEX_2D(i,j)];
                Fy_lava_i = Fy_lava[INDEX_2D(i,j)];
                Fx_air_i = Fx_air[INDEX_2D(i,j)];
                Fy_air_i = Fy_air[INDEX_2D(i,j)];
                u2 = u_i*u_i + v_i*v_i;
                for (int p = 0; p < NP; p++) {
                    uc = u_i*cx[p] + v_i*cy[p];
                    feq = w[p]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    sx = w[p]*((cx[p]-u_i)/cs2 + uc/(cs2*cs2)*cx[p]);
                    sy = w[p]*((cy[p]-v_i)/cs2 + uc/(cs2*cs2)*cy[p]);
                    S = (1.0 - 1.0/(2.0*tau_lava))*(sx*Fx_lava_i + sy*Fy_lava_i);
                    f2_lava[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_lava)*f1_lava[INDEX_3D(i,j,p)] + 1.0/tau_lava*rho_lava_i*feq + S;
                    S = (1.0 - 1.0/(2.0*tau_air))*(sx*Fx_air_i + sy*Fy_air_i);
                    f2_air[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_air)*f1_air[INDEX_3D(i,j,p)] + 1.0/tau_air*rho_air_i*feq + S;
                }
            }
        }

        // Stream populations while taking boundary conditions in mind
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int p = 0; p < NP; p++) {
                    x_i = mod(i-cx_i[p], NX);
                    y_i = j-cy_i[p];
                    p_i = p;
                    if ((y_i < 0) || (y_i == NY)) {
                        x_i = i;
                        y_i = j;
                        p_i = ((p+3)%8)+1;
                    }
                    f1_lava[INDEX_3D(i,j,p)] = f2_lava[INDEX_3D(x_i, y_i, p_i)];
                    f1_air[INDEX_3D(i,j,p)] = f2_air[INDEX_3D(x_i, y_i, p_i)];
                }
            }
        }

        // Compute macroscopic quantities
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                rho_i = 0.0;
                rho_lava_i = 0.0;
                rho_air_i = 0.0;
                u_i = 0.0;
                v_i = 0.0;
                for (int p = 0; p < NP; p++) {
                    rho_i += f1_lava[INDEX_3D(i,j,p)] + f1_air[INDEX_3D(i,j,p)];
                    rho_lava_i += f1_lava[INDEX_3D(i,j,p)];
                    rho_air_i += f1_air[INDEX_3D(i,j,p)];
                    u_i += (f1_lava[INDEX_3D(i,j,p)] + f1_air[INDEX_3D(i,j,p)])*cx[p];
                    v_i += (f1_lava[INDEX_3D(i,j,p)] + f1_air[INDEX_3D(i,j,p)])*cy[p];
                }
                rho[INDEX_2D(i,j)] = rho_i;
                rho_lava[INDEX_2D(i,j)] = rho_lava_i;
                rho_air[INDEX_2D(i,j)] = rho_air_i;

                u[INDEX_2D(i,j)] = u_i / rho_i;
                v[INDEX_2D(i,j)] = v_i / rho_i;
            }
        }

        // Calculate the force fields
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                Fx_lava_i = 0.0;
                Fy_lava_i = 0.0;
                Fx_air_i = 0.0;
                Fy_air_i = 0.0;
                for (int p = 0; p < NP; p++) {
                    x_i = mod(i + cx[p], NX);
                    y_i = j + cy[p];
                    if ((y_i < 0) || (y_i == NY)) {
                        x_i = i;
                        y_i = j;
                    }
                    Fx_lava_i += w[p]*rho_air[INDEX_2D(x_i,y_i)]*cx[p];
                    Fy_lava_i += w[p]*rho_air[INDEX_2D(x_i,y_i)]*cy[p];
                    Fx_air_i += w[p]*rho_lava[INDEX_2D(x_i,y_i)]*cx[p];
                    Fy_air_i += w[p]*rho_lava[INDEX_2D(x_i,y_i)]*cy[p];
                }
                Fx_lava[INDEX_2D(i,j)] = -G*rho_lava[INDEX_2D(i,j)]*Fx_lava_i + rho_lava[INDEX_2D(i,j)]/rho[INDEX_2D(i,j)]*F_p;
                Fy_lava[INDEX_2D(i,j)] = -G*rho_lava[INDEX_2D(i,j)]*Fy_lava_i;
                Fx_air[INDEX_2D(i,j)] = -G*rho_air[INDEX_2D(i,j)]*Fx_air_i + rho_air[INDEX_2D(i,j)]/rho[INDEX_2D(i,j)]*F_p;
                Fy_air[INDEX_2D(i,j)] = -G*rho_air[INDEX_2D(i,j)]*Fy_air_i;
            }
        }

        // Shift the velocity fields
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                u[INDEX_2D(i,j)] += 1.0/(2.0*rho[INDEX_2D(i,j)])*(Fx_lava[INDEX_2D(i,j)] + Fx_air[INDEX_2D(i,j)]);
                v[INDEX_2D(i,j)] += 1.0/(2.0*rho[INDEX_2D(i,j)])*(Fy_lava[INDEX_2D(i,j)] + Fy_air[INDEX_2D(i,j)]);
            }
        }

        // Log the progress
        if (((t+1) % NLOG) == 0) {
            printf("\rProgress: %.1f%%", ((double)(t+1))/((double)NTIME)*100.0);
            fflush(stdout);
        }
    }

    // Output the final macroscopic quantities
    output_data(NTIME, noutput);

    printf("\nDone!\n");
    return 0;
}