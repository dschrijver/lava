#ifndef STENCIL_H
#define STENCIL_H


#include "../constants.h"

// --- D2Q9 STENCIL ---
static int cx_i[NP] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
static int cy_i[NP] = {0, 0, 0, 1, -1, 1, -1, -1, 1};
static double cx[NP] = {0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
static double cy[NP] = {0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
static int p_bounceback[NP] = {0, 2, 1, 4, 3, 6, 5, 8, 7};

#endif