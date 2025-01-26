#include <stdio.h>  // printf(), fflush(), stdout()

#include "constants.h"
#include "include/animate.h"
#include "include/memory.h"
#include "include/initialize.h"
#include "include/output.h"
#include "include/collide.h"
#include "include/stream.h"
#include "include/macroscopic.h"
#include "include/forcing.h"
#include "include/main.h"


int main(void) {

#ifdef ANIMATE
    initialize_animation();
#endif

    allocate_memory();

    initialize();

#ifdef FLOW
    calculate_body_forces();
#endif

#ifdef OUTPUT
    output_data(0);
#endif

#ifdef FLOW
    shift_velocity();
#endif

    // --- MAIN LOOP ---
    for (int t = 0; t < NTIME; t++) {
        
#ifdef OUTPUT
        if (((t % NSTORE) == 0) && (t > 0)) {
            output_data(t);
        }
#endif

#ifdef ANIMATE
        if ((t % NFRAME) == 0) {
            if (!render_frame()) break;
        }
#endif

#ifdef TEMPERATURE
        collide_thermal_populations();
#endif

#ifdef TEMPERATURE
        stream_thermal_populations();
#endif

#ifdef TEMPERATURE
        calculate_thermal_macroscopic_quantities();
#endif

#ifdef PHASECHANGE
        calculate_enthalpy_and_liquid_fraction();
#endif

#ifdef FLOW
        calculate_body_forces();
#endif

#ifdef FLOW
        collide_hydrodynamic_populations();
#endif

#ifdef FLOW
        stream_hydrodynamic_populations();
#endif

#ifdef FLOW
        calculate_hydrodynamic_macroscopic_quantities();
#endif

        if (((t+1) % NLOG) == 0) {
            printf("\rProgress: %.1f%%", ((double)(t+1))/((double)NTIME)*100.0);
            fflush(stdout);
        }
    }

#ifdef OUTPUT
    output_data(NTIME);
#endif

    free_memory();

#ifdef ANIMATE
    close_animation();
#endif

    printf("\nDone!\n");

    return 0;
}