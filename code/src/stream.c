#include "../constants.h"
#include "../include/globals.h"
#include "../include/stream.h"


void stream_hydrodynamic_populations() {

    int x_i, y_i;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int p = 0; p < NP; p++) {

                x_i = i-cx_i[p];
                y_i = j-cy_i[p];

#ifdef XPERIODIC
                if (x_i < 0) x_i = NX-1;
                else if (x_i == NX) x_i = 0;
#endif

#ifdef SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i < 0) {
                    f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(i, j, p_bounceback[p])];
                    continue;
                }
#endif

#ifdef NORTH_NOSLIP_HALFWAY_BOUNCEBACK
                if (y_i < 0) {
                    f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(i, j, p_bounceback[p])];
                    continue;
                }
#endif

                f1[INDEX_3D(i,j,p)] = f2[INDEX_3D(x_i, y_i, p)];
                
            }
        }
    }

}