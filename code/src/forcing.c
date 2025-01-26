#include "../constants.h"
#include "../include/globals.h"
#include "../include/forcing.h"



void calculate_body_forces() {

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

#ifdef FLOW

            Fx[INDEX_2D(i,j)] = Fx_body;

    #ifdef TEMPERATURE
            Fy[INDEX_2D(i,j)] = Fy_body + alpha * (T[INDEX_2D(i,j)] - T_top) * g;
    #else
            Fy[INDEX_2D(i,j)] = Fy_body;
    #endif

#endif

        }
    }

}