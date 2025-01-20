#include <stdio.h>  // printf, fflush, stdout
#include <stdlib.h> // malloc
#include <hdf5.h>   // Output data_*.h5 files
#include <time.h>   // time
#include <math.h>   // sin
#include <raylib.h> 


// --- SETTINGS ---
#define NTIME   1000000          // Number of timesteps
#define NSTORE  10000            // Store macroscopic quantities after NSTORE timesteps
#define NLOG    100             // Print progress percentage after NLOG timesteps

#define NX 100                  // Number of cells in the x-direction
#define NY 50                  // Number of cells in the y-direction
#define NP 9                    // Number of velocity directions, DON'T CHANGE!

static double tau = 0.5012639224659765;
static double tau_g = 0.5017801724872908;
static double alpha = 207e-6;
static double T_top = 293.0;
static double T_bottom = 294.0;
static double g = 0.009661835748792274;


// --- DISPLAY ---
#define ANIM
#undef  STORE
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 400;
const int cell_size = 8;
const int steps_per_frame = 100;


// --- SIMULATION ---
#define M_PI 3.14159265358979323846
#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)

static double cs2 = 1.0/3.0;
static int cx_i[NP] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
static int cy_i[NP] = {0, 0, -1, -1, -1, 0, 1, 1, 1};
static double cx[NP] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
static double cy[NP] = {0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};

double *rho, *u, *v, *Fy, *f1, *f2;
double *T, *g1, *g2;


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
    dataset_id = H5Dcreate2(hdf5_fp, "Fy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fy);
    dataset_id = H5Dcreate2(hdf5_fp, "T", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, T);

    H5Fflush(hdf5_fp, H5F_SCOPE_GLOBAL);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(dataspace_id_1d);
    H5Fclose(hdf5_fp);
}

int render_frame() {
    char T_i;
    Color cell_color = {0,0,0,255};

    PollInputEvents(); // Poll input events (SUPPORT_CUSTOM_FRAME_CONTROL)
    if (WindowShouldClose()) {
        return 0;
    }

    BeginDrawing();

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                T_i = (char)(fmin(fmax(T[INDEX_2D(i,j)]-293.0, 0.0), 1.0) * 255.0);
                cell_color.r = T_i;
                cell_color.b = 255 - T_i;
                DrawRectangle(i*cell_size, (NY-1-j)*cell_size, cell_size, cell_size, cell_color);
            }
        }

    EndDrawing();
    SwapScreenBuffer();

    return 1;
}

// Used for periodic boundary conditions
int mod(int x, int n) {
    if (x < 0) return n-1;
    if (x == n) return 0;
    return x;
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
    double Fy_i;
    double T_i;
    double u2, uc, eq, S;
    int x_i, y_i, p_i;

    // Allocate memory for macroscopic quantities
    rho = (double *) malloc(NX*NY*sizeof(double));
    u = (double *) malloc(NX*NY*sizeof(double));
    v = (double *) malloc(NX*NY*sizeof(double));
    Fy = (double *) malloc(NX*NY*sizeof(double));
    T = (double *) malloc(NX*NY*sizeof(double));

    // Allocate memory for populations
    f1 = (double *) malloc(NX*NY*NP*sizeof(double));
    f2 = (double *) malloc(NX*NY*NP*sizeof(double));
    g1 = (double *) malloc(NX*NY*NP*sizeof(double));
    g2 = (double *) malloc(NX*NY*NP*sizeof(double));

    // Choose initial density and velocity fields and compute initial populations
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            rho[INDEX_2D(i,j)] = 1.0;
            y_i = 0.5+j;
            x_i = 0.5+i;
            T[INDEX_2D(i,j)] = (T_top-T_bottom)/NY*y_i + T_bottom;
            u[INDEX_2D(i,j)] = 0.0001*sin(2.0*M_PI*((double)x_i)/((double)NX));
            v[INDEX_2D(i,j)] = 0.0001*sin(2.0*M_PI*((double)y_i)/((double)NY));
            u2 = u[INDEX_2D(i,j)]*u[INDEX_2D(i,j)] + v[INDEX_2D(i,j)]*v[INDEX_2D(i,j)];
            for (int p = 0; p < NP; p++) {
                uc = u[INDEX_2D(i,j)]*cx[p] + v[INDEX_2D(i,j)]*cy[p];
                f1[INDEX_3D(i,j,p)] = w[p]*rho[INDEX_2D(i,j)]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                g1[INDEX_3D(i,j,p)] = w[p]*T[INDEX_2D(i,j)]*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
            }
        }
    }

    // Initialize the force fields
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            Fy[INDEX_2D(i,j)] = alpha * (T[INDEX_2D(i,j)] - T_top) * g;
        }
    }

#ifdef STORE
    // Output initial macroscopic quantities
    output_data(0, noutput);
    noutput++;
#endif

    // Shift the velocity fields
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            v[INDEX_2D(i,j)] -= 1.0/(2.0*rho[INDEX_2D(i,j)])*Fy[INDEX_2D(i,j)];
        }
    }

    // MAIN LOOP
    for (int t = 0; t < NTIME; t++) {

#ifdef ANIM
        if ((t % steps_per_frame) == 0) {
            if (!render_frame()) {
                break;
            }
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
                for (int p = 0; p < NP; p++) {
                    uc = u_i*cx[p] + v_i*cy[p];
                    eq = w[p]*T_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    g2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau_g)*g1[INDEX_3D(i,j,p)] + 1.0/tau_g*eq;
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

        // determine buoyant force
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                Fy[INDEX_2D(i,j)] = alpha * (T[INDEX_2D(i,j)] - T_top) * g;
            }
        }

        // Collide populations
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                rho_i = rho[INDEX_2D(i,j)];
                u_i = u[INDEX_2D(i,j)];
                v_i = v[INDEX_2D(i,j)];
                Fy_i = Fy[INDEX_2D(i,j)];
                u2 = u_i*u_i + v_i*v_i;
                for (int p = 0; p < NP; p++) {
                    uc = u_i*cx[p] + v_i*cy[p];
                    eq = w[p]*rho_i*(1.0 + uc/cs2 + uc*uc/(2.0*cs2*cs2) - u2/(2.0*cs2));
                    S = (1.0 - 1.0/(2.0*tau))*w[p]*((cy[p]-v_i)/cs2 + uc/(cs2*cs2)*cy[p])*Fy_i;
                    f2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau)*f1[INDEX_3D(i,j,p)] + 1.0/tau*eq + S;
                }
            }
        }

        // Stream populations
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
                    f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(x_i, y_i, p_i)];
                }
            }
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
                u[INDEX_2D(i,j)] = u_i / rho_i;
                v[INDEX_2D(i,j)] = (v_i + 0.5*Fy[INDEX_2D(i,j)]) / rho_i;
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
