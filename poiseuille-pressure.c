#include <stdio.h>  // printf, fflush, stdout
#include <stdlib.h> // malloc
#include <hdf5.h>


#define NTIME   100000
#define NSTORE  1000
#define NLOG    100

#define NX 200
#define NY 200
#define NP 9

#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)

static int cx_i[NP] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
static int cy_i[NP] = {0, 0, -1, -1, -1, 0, 1, 1, 1};
static double cx[NP] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
static double cy[NP] = {0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};

static double cs2 = 1.0/3.0;
static double tau = 1.0;
static double Delta_p = 1e-5;

double *rho, *u, *v;
double *f1, *f2;

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

    H5Fflush(hdf5_fp, H5F_SCOPE_GLOBAL);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(dataspace_id_1d);
    H5Fclose(hdf5_fp);
}

int main(void) {
    double feq;
    double rho_i, u_i, v_i, u2, uc;
    int x_i, y_i, p_i;
    
    int noutput = 0;

    rho = (double *) malloc(NX*NY*sizeof(double));
    u = (double *) malloc(NX*NY*sizeof(double));
    v = (double *) malloc(NX*NY*sizeof(double));

    f1 = (double *) malloc(NX*NY*NP*sizeof(double));
    f2 = (double *) malloc(NX*NY*NP*sizeof(double));

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            rho[INDEX_2D(i,j)] = 1.0 - (Delta_p/cs2)/(double)NX * (0.5+(double)i);
            u[INDEX_2D(i,j)] = 0.0;
            v[INDEX_2D(i,j)] = 0.0;
            for (int p = 0; p < NP; p++) {
                f1[INDEX_3D(i,j,p)] = w[p]*rho[INDEX_2D(i,j)];
            }
        }
    }

    output_data(0, noutput);
    noutput++;

    for (int t = 0; t < NTIME; t++) {

        if (((t % NSTORE) == 0) && (t > 0)) {
            output_data(t, noutput);
            noutput++;
        }

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                rho_i = rho[INDEX_2D(i,j)];
                u_i = u[INDEX_2D(i,j)];
                v_i = v[INDEX_2D(i,j)];
                u2 = u_i*u_i + v_i*v_i;
                for (int p = 0; p < NP; p++) {
                    uc = u_i*cx[p] + v_i*cy[p];
                    feq = w[p]*rho_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    f2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau)*f1[INDEX_3D(i,j,p)] + 1.0/tau*feq;
                }
            }
        }

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                u_i = u[INDEX_2D(i,j)];
                v_i = v[INDEX_2D(i,j)];
                u2 = u_i*u_i + v_i*v_i;
                for (int p = 0; p < NP; p++) {
                    x_i = i-cx_i[p];
                    y_i = j-cy_i[p];
                    p_i = ((p+3)%8)+1;
                    if (x_i < 0) {
                        uc = u_i*cx[p] + v_i*cy[p];
                        f1[INDEX_3D(i,j,p)] = -f2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*(1.0 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    }
                    else if (x_i == NX) {
                        uc = u_i*cx[p] + v_i*cy[p];
                        f1[INDEX_3D(i,j,p)] = -f2[INDEX_3D(i, j, p_i)] + 2.0*w[p_i]*(1.0-Delta_p/cs2)*(1.0 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
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
                u[INDEX_2D(i,j)] = u_i / rho_i;
                v[INDEX_2D(i,j)] = v_i / rho_i;
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