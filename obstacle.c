#include <stdio.h>  // printf, fflush, stdout
#include <stdlib.h> // malloc


#define NTIME   17500
#define NLOG    100
#define NSTORE  17500

#define NX 240
#define NY 80
#define NP 9

#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)

static int cx_i[NP] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
static int cy_i[NP] = {0, 0, -1, -1, -1, 0, 1, 1, 1};
static double cx[NP] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
static double cy[NP] = {0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};

static double cs2 = 1.0/3.0;
static double tau = 0.52;
static double F = 1e-5;

static double obstacle_width = 20;
static double obstacle_height = 40;
static double obstacle_shift = 1;

int *obstacle_grid;
double *rho, *u, *v;
double *f1, *f2;

int mod(int x, int n) {
    if (x < 0) return n-1;
    if (x == n) return 0;
    return x;
}

int main(void) {
    obstacle_grid = (int *) malloc(NX*NY*sizeof(int));

    rho = (double *) malloc(NX*NY*sizeof(double));
    u = (double *) malloc(NX*NY*sizeof(double));
    v = (double *) malloc(NX*NY*sizeof(double));

    f1 = (double *) malloc(NX*NY*NP*sizeof(double));
    f2 = (double *) malloc(NX*NY*NP*sizeof(double));

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            if ((i > NX/2 - obstacle_width/2) && (i < NX/2 + obstacle_width/2) &&
                (j > NY/2 - obstacle_height/2 + obstacle_shift) && (j < NY/2 + obstacle_height/2 + obstacle_shift)) {
                obstacle_grid[INDEX_2D(i,j)] = 1;
            } 
            else {
                obstacle_grid[INDEX_2D(i,j)] = 0;
            }
        }
    }

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            if (obstacle_grid[INDEX_2D(i,j)]) {
                rho[INDEX_2D(i,j)] = 0.0;
            }
            else {
                rho[INDEX_2D(i,j)] = 1.0;
            }
            u[INDEX_2D(i,j)] = 0.0;
            v[INDEX_2D(i,j)] = 0.0;
            for (int p = 0; p < NP; p++) {
                f1[INDEX_3D(i,j,p)] = w[p]*rho[INDEX_2D(i,j)];
            }
        }
    }

    FILE *fp;
    FILE *fp_t = fopen("time.dat", "w");
    fclose(fp_t);
    char filename[32];
    int s = 0;
    double feq, S;
    double rho_i, u_i, v_i, u2, uc;
    int x_i, y_i, p_i;
    for (int t = 0; t < NTIME; t++) {

        if ((t % NSTORE) == 0) {
            sprintf(filename, "u_%d.dat", s);
            fp = fopen(filename, "w");
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    fprintf(fp, "%-23.15e", u[INDEX_2D(i,j)]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            s++;
            fp_t = fopen("time.dat", "a");
            fprintf(fp_t, "%d\n", t);
            fclose(fp_t);
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
                    S = (1.0 - 1.0/(2.0*tau))*w[p]*((cx[p]-u_i)/cs2 + uc/(cs2*cs2)*cx[p])*F;
                    f2[INDEX_3D(i,j,p)] = (1.0 - 1.0/tau)*f1[INDEX_3D(i,j,p)] + 1.0/tau*feq + S;
                }
            }
        }

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int p = 0; p < NP; p++) {
                    if (obstacle_grid[INDEX_2D(i,j)]) {
                        f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(i,j,p)];
                    }
                    else {
                        x_i = i-cx_i[p];
                        y_i = j-cy_i[p];
                        p_i = p;
                        x_i = mod(x_i, NX);
                        if ((y_i < 0) || (y_i == NY)) {
                            x_i = i;
                            y_i = j;
                            p_i = ((p+3)%8)+1;
                        }
                        else if (obstacle_grid[INDEX_2D(x_i, y_i)]) {
                            x_i = i;
                            y_i = j;
                            p_i = ((p+3)%8)+1;
                        }
                        f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(x_i, y_i, p_i)];
                    }
                }
            }
        }

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                rho_i = 0.0;
                u_i = 0.0;
                v_i = 0.0;
                if (!obstacle_grid[INDEX_2D(i,j)]) {
                    for (int p = 0; p < NP; p++) {
                        rho_i += f1[INDEX_3D(i,j,p)];
                        u_i += f1[INDEX_3D(i,j,p)]*cx[p];
                        v_i += f1[INDEX_3D(i,j,p)]*cy[p];
                    }
                    rho[INDEX_2D(i,j)] = rho_i;
                    u[INDEX_2D(i,j)] = (u_i + 0.5*F) / rho_i;
                    v[INDEX_2D(i,j)] = v_i / rho_i;
                }
            }
        }

        if ((t % NLOG) == 0) {
            printf("\rProgress: %.1f%%", ((double)(t+1))/((double)NTIME)*100.0);
            fflush(stdout);
        }
    }

    sprintf(filename, "u_%d.dat", s);
    fp = fopen(filename, "w");
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            fprintf(fp, "%-23.15e ", u[INDEX_2D(i,j)]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    fp_t = fopen("time.dat", "a");
    fprintf(fp_t, "%d\n", NTIME);
    fclose(fp_t);

    printf("\nDone!\n");

    return 0;
}